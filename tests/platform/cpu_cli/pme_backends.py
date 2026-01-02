"""
PME energy backend for gcmc_cpu ("Mode D": energy_method=pme).

These tests are end-to-end and do not depend on stdout/stderr log text.
They assert behavior via --dump-accept (JSONL) and --dump-params (JSON).
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from cpu_cli.inp_units_compat_helpers import _first_accept_record, _run_gcmc_cpu, _write_inp


def _write_text(path: Path, content: str) -> None:
    path.write_text(content.strip() + "\n")


def _read_dump_params(path: Path) -> dict:
    return json.loads(path.read_text())


def test_pme_adds_long_range_coulomb_beyond_cutoff_for_deletion(gcmc_cpu, temp_dir):
    """
    For a charged pair separated beyond the DIRECT cutoff:
    - DIRECT+cutoff deletion sees ~0 interaction energy -> deltaU ~ 0
    - PME deletion includes long-range electrostatics -> deltaU is significantly positive
    """
    work = Path(temp_dir) / "pme_backend" / "long_range_deletion"
    work.mkdir(parents=True, exist_ok=True)

    def _write_pdb(path: Path, *, frg_x_ang: float) -> None:
        _write_text(
            path,
            f"""
CRYST1   40.000   40.000   40.000  90.00  90.00  90.00 P 1           1
ATOM      1  Q   HST A   1      10.000  20.000  20.000  1.00  0.00           C
ATOM      2  Q   FRG B   1      {frg_x_ang:6.3f}  20.000  20.000  1.00  0.00           C
END
""",
        )

    pdb_r16 = work / "sys_r16.pdb"  # min-image distance 16 Å -> 1.6 nm
    pdb_r20 = work / "sys_r20.pdb"  # min-image distance 20 Å -> 2.0 nm
    _write_pdb(pdb_r16, frg_x_ang=26.0)
    _write_pdb(pdb_r20, frg_x_ang=30.0)

    top = work / "sys.top"
    _write_text(
        top,
        """
[ defaults ]
1 2 yes 1.0 1.0

[ moleculetype ]
HST  1

[ atoms ]
; nr  type  resnr  residue  atom  cgnr  charge  mass
1   Q     1      HST      Q     1     1.000   12.011

[ moleculetype ]
FRG  1

[ atoms ]
; nr  type  resnr  residue  atom  cgnr  charge  mass
1   Q     1      FRG      Q     1    -1.000   12.011

[ system ]
Minimal

[ molecules ]
HST  1
FRG  1
""",
    )

    par = work / "par.itp"
    _write_text(
        par,
        """
[ defaults ]
1 2 yes 1.0 1.0

[ atomtypes ]
; name  at.num  mass   charge  ptype  sigma   epsilon
Q       0       12.011 0.000   A      0.300   0.000
""",
    )

    frag = work / "frag.itp"
    _write_text(
        frag,
        """
[ moleculetype ]
FRG  1

[ atoms ]
; nr  type  resnr  residue  atom  cgnr  charge  mass
1   Q     1      FRG      Q     1    -1.000  12.011
""",
    )

    def _run(energy_method: str, *, pdb_path: Path) -> tuple[float, dict]:
        run_dir = work / f"{energy_method}_{pdb_path.stem}"
        run_dir.mkdir(parents=True, exist_ok=True)

        inp = run_dir / "test.inp"
        _write_inp(
            inp,
            f"""
random_seed:123
energy_method:{energy_method}
par:{par}
fragitp:{frag}

fragname:FRG
fragconc:0.01
fragmuex:0.0

pdb:{pdb_path}
top:{top}
box_size:40.0 40.0 40.0
cutoff:6.0
temperature:300.0

moves_per_step:1
mcsteps:1
nprint:1
mc_move_prob:0 1 0 0

use_cavity_bias:no
use_conf_bias:no
""",
        )

        out_prefix = run_dir / "out" / "gcmc"
        out_prefix.parent.mkdir(parents=True, exist_ok=True)
        accept_log = run_dir / "out" / "acceptance.jsonl"
        dump_params = run_dir / "out" / "params.json"

        result = _run_gcmc_cpu(
            gcmc_cpu,
            workdir=run_dir,
            inp=inp,
            out_prefix=out_prefix,
            extra_args=["--dump-accept", str(accept_log), "--dump-params", str(dump_params)],
            timeout=60,
        )
        assert result.returncode == 0, result.stdout + result.stderr

        params = _read_dump_params(dump_params)
        assert params["basic"]["energy_method"] == energy_method

        rec = _first_accept_record(accept_log, move="deletion", species="FRG")
        assert int(rec.get("nBefore", 0)) == 1
        return float(rec["deltaU"]), rec

    delta_direct, _ = _run("direct", pdb_path=pdb_r16)
    delta_pme_r16, _ = _run("pme", pdb_path=pdb_r16)
    delta_pme_r20, _ = _run("pme", pdb_path=pdb_r20)

    assert abs(delta_direct) < 1e-3
    assert abs(delta_pme_r16) > 20.0
    assert abs(delta_pme_r16 - delta_direct) > 20.0
    assert abs(delta_pme_r16 - delta_pme_r20) > 5.0

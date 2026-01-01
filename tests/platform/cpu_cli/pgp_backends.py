"""
PGP energy backends for gcmc_cpu (Mode B: pgp_host, Mode C: pgp_full).

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


def test_dump_params_includes_energy_method(gcmc_cpu, temp_dir):
    work = Path(temp_dir) / "pgp_dump_params_energy_method"
    work.mkdir(parents=True, exist_ok=True)

    pdb = work / "sys.pdb"
    _write_text(
        pdb,
        """
CRYST1   50.000   50.000   50.000  90.00  90.00  90.00 P 1           1
ATOM      1  Q   MOL A   1      10.000  10.000  10.000  1.00  0.00           C
END
""",
    )

    top = work / "sys.top"
    _write_text(
        top,
        """
[ defaults ]
1 2 yes 1.0 1.0

[ moleculetype ]
MOL  1

[ atoms ]
; nr  type  resnr  residue  atom  cgnr  charge  mass
1   Q     1      MOL      Q     1     1.000   12.011

[ system ]
Minimal

[ molecules ]
MOL  1
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
1   Q     1      FRG      Q     1     -1.000  12.011
""",
    )

    inp = work / "test.inp"
    _write_inp(
        inp,
        f"""
random_seed:123
energy_method:pgp_host
par:{par}
fragitp:{frag}
fragname:FRG
fragconc:1.0
fragmuex:0.0

pdb:{pdb}
top:{top}
box_size:50.0 50.0 50.0
cutoff:6.0
temperature:300.0
moves_per_step:1
mcsteps:1
nprint:1
mc_move_prob:1 0 0 0

use_cavity_bias:no
use_conf_bias:no
""",
    )

    out_prefix = work / "out" / "gcmc"
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    dump_params = work / "out" / "params.json"

    result = _run_gcmc_cpu(
        gcmc_cpu,
        workdir=work,
        inp=inp,
        out_prefix=out_prefix,
        extra_args=["--dump-params", str(dump_params)],
        timeout=40,
    )
    assert result.returncode == 0, result.stdout + result.stderr

    params = _read_dump_params(dump_params)
    assert params["basic"]["energy_method"] == "pgp_host"


def test_pgp_host_fails_without_fixed_residue(gcmc_cpu, temp_dir):
    work = Path(temp_dir) / "pgp_host_requires_fixed"
    work.mkdir(parents=True, exist_ok=True)

    pdb = work / "sys.pdb"
    _write_text(
        pdb,
        """
CRYST1   50.000   50.000   50.000  90.00  90.00  90.00 P 1           1
ATOM      1  Q   FRG A   1      10.000  10.000  10.000  1.00  0.00           C
END
""",
    )

    top = work / "sys.top"
    _write_text(
        top,
        """
[ defaults ]
1 2 yes 1.0 1.0

[ moleculetype ]
FRG  1

[ atoms ]
; nr  type  resnr  residue  atom  cgnr  charge  mass
1   Q     1      FRG      Q     1     1.000   12.011

[ system ]
Minimal

[ molecules ]
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
1   Q     1      FRG      Q     1     1.000  12.011
""",
    )

    def _write_common_inp(path: Path, *, energy_method: str) -> None:
        _write_inp(
            path,
            f"""
random_seed:123
energy_method:{energy_method}
par:{par}
fragitp:{frag}
fragname:FRG
fragconc:0.0
fragmuex:0.0

pdb:{pdb}
top:{top}
box_size:50.0 50.0 50.0
cutoff:6.0
temperature:300.0
moves_per_step:1
mcsteps:1
nprint:1
mc_move_prob:0 0 1 0

use_cavity_bias:no
use_conf_bias:no
""",
        )

    out_prefix = work / "out" / "gcmc"
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    direct_inp = work / "direct.inp"
    _write_common_inp(direct_inp, energy_method="direct")
    direct = _run_gcmc_cpu(gcmc_cpu, workdir=work, inp=direct_inp, out_prefix=out_prefix, timeout=40)
    assert direct.returncode == 0, direct.stdout + direct.stderr

    pgp_inp = work / "pgp_host.inp"
    _write_common_inp(pgp_inp, energy_method="pgp_host")
    pgp = _run_gcmc_cpu(gcmc_cpu, workdir=work, inp=pgp_inp, out_prefix=out_prefix, timeout=40)
    assert pgp.returncode != 0


def test_pgp_host_adds_long_range_coulomb_beyond_cutoff(gcmc_cpu, temp_dir):
    work = Path(temp_dir) / "pgp_host_long_range_coulomb"
    work.mkdir(parents=True, exist_ok=True)

    # Fixed host residue (resname != fragname) at (10,10,10) Å.
    pdb = work / "host.pdb"
    _write_text(
        pdb,
        """
CRYST1   50.000   50.000   50.000  90.00  90.00  90.00 P 1           1
ATOM      1  Q   HST A   1      10.000  10.000  10.000  1.00  0.00           C
END
""",
    )

    top = work / "host.top"
    _write_text(
        top,
        """
[ defaults ]
1 2 yes 1.0 1.0

[ moleculetype ]
HST  1

[ atoms ]
; nr  type  resnr  residue  atom  cgnr  charge  mass
1   QP    1      HST      Q     1     1.000   12.011

[ system ]
Host

[ molecules ]
HST  1
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
QP      0       12.011 0.000   A      0.300   0.000
QM      0       12.011 0.000   A      0.300   0.000
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
1   QM    1      FRG      Q     1     -1.000  12.011
""",
    )

    def _run(method: str) -> float:
        out_prefix = work / method / "out" / "gcmc"
        out_prefix.parent.mkdir(parents=True, exist_ok=True)
        accept_log = work / method / "out" / "accept.jsonl"

        inp = work / method / "test.inp"
        inp.parent.mkdir(parents=True, exist_ok=True)
        _write_inp(
            inp,
            f"""
random_seed:12345
energy_method:{method}
par:{par}
fragitp:{frag}
fragname:FRG
fragconc:1.0
fragmuex:0.0

pdb:{pdb}
top:{top}
box_size:50.0 50.0 50.0
gcmc_region:box 39.5 9.5 9.5 40.5 10.5 10.5
cutoff:6.0
temperature:300.0
moves_per_step:1
mcsteps:1
nprint:1
mc_move_prob:1 0 0 0

use_cavity_bias:no
use_conf_bias:no
""",
        )

        result = _run_gcmc_cpu(
            gcmc_cpu,
            workdir=work / method,
            inp=inp,
            out_prefix=out_prefix,
            extra_args=["--dump-accept", str(accept_log)],
            timeout=60,
        )
        assert result.returncode == 0, result.stdout + result.stderr
        rec = _first_accept_record(accept_log, move="insertion", species="FRG")
        return float(rec["deltaU"])

    delta_direct = _run("direct")
    assert abs(delta_direct) < 1e-6

    delta_pgp = _run("pgp_host")
    assert delta_pgp < -1.0


def test_pgp_full_adds_long_range_coulomb_beyond_cutoff_for_nonfixed_background(gcmc_cpu, temp_dir):
    work = Path(temp_dir) / "pgp_full_long_range_coulomb"
    work.mkdir(parents=True, exist_ok=True)

    # One existing FRG in the system (non-fixed), insertion region far away so DIRECT+cutoff sees 0 interaction.
    pdb = work / "sys.pdb"
    _write_text(
        pdb,
        """
CRYST1   50.000   50.000   50.000  90.00  90.00  90.00 P 1           1
ATOM      1  Q   FRG A   1      10.000  10.000  10.000  1.00  0.00           C
END
""",
    )

    top = work / "sys.top"
    _write_text(
        top,
        """
[ defaults ]
1 2 yes 1.0 1.0

[ moleculetype ]
FRG  1

[ atoms ]
; nr  type  resnr  residue  atom  cgnr  charge  mass
1   QP    1      FRG      Q     1     1.000   12.011

[ system ]
FRG

[ molecules ]
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
QP      0       12.011 0.000   A      0.300   0.000
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
1   QP    1      FRG      Q     1     1.000  12.011
""",
    )

    def _run(method: str) -> float:
        out_prefix = work / method / "out" / "gcmc"
        out_prefix.parent.mkdir(parents=True, exist_ok=True)
        accept_log = work / method / "out" / "accept.jsonl"

        inp = work / method / "test.inp"
        inp.parent.mkdir(parents=True, exist_ok=True)
        _write_inp(
            inp,
            f"""
random_seed:12345
energy_method:{method}
par:{par}
fragitp:{frag}
fragname:FRG
fragconc:1.0
fragmuex:0.0

pdb:{pdb}
top:{top}
box_size:50.0 50.0 50.0
gcmc_region:box 39.5 9.5 9.5 40.5 10.5 10.5
cutoff:6.0
temperature:300.0
moves_per_step:1
mcsteps:1
nprint:1
mc_move_prob:1 0 0 0

use_cavity_bias:no
use_conf_bias:no
""",
        )

        result = _run_gcmc_cpu(
            gcmc_cpu,
            workdir=work / method,
            inp=inp,
            out_prefix=out_prefix,
            extra_args=["--dump-accept", str(accept_log)],
            timeout=60,
        )
        assert result.returncode == 0, result.stdout + result.stderr
        rec = _first_accept_record(accept_log, move="insertion", species="FRG")
        return float(rec["deltaU"])

    delta_direct = _run("direct")
    assert abs(delta_direct) < 1e-6

    delta_pgp = _run("pgp_full")
    assert delta_pgp > 1.0


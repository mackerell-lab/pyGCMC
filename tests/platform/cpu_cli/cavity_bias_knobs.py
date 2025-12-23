"""
Cavity-bias knob regressions (probe_radius / exclude_hydrogens / use_vdw_radius).
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from .inp_units_compat_helpers import _run_gcmc_cpu, _write_inp


def _write_text(path: Path, content: str) -> None:
    path.write_text(content.strip() + "\n")


def _mean_cavity_fraction(accept_log: Path) -> float:
    records = [json.loads(line) for line in accept_log.read_text().splitlines() if line.strip()]
    insertions = [r for r in records if str(r.get("move", "")).lower() == "insertion"]
    assert insertions, "No insertion records found in acceptance log"
    fractions = [float(r.get("cavityFraction", 0.0)) for r in insertions]
    assert all(0.0 <= val <= 1.0 for val in fractions), "Invalid cavityFraction values"
    return sum(fractions) / len(fractions)


def _run_cavity_case(
    *,
    gcmc_cpu: str,
    work: Path,
    inp: Path,
    out_prefix: Path,
    accept_log: Path,
    timeout: int = 60,
) -> float:
    result = _run_gcmc_cpu(
        gcmc_cpu,
        workdir=work,
        inp=inp,
        out_prefix=out_prefix,
        extra_args=["--dump-accept", str(accept_log), "--store-probabilities"],
        timeout=timeout,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert accept_log.exists()
    return _mean_cavity_fraction(accept_log)


def test_cavity_bias_probe_radius_changes_cavity_fraction(gcmc_cpu, temp_dir):
    work = Path(temp_dir) / "cavity_bias_probe_radius"
    work.mkdir(parents=True, exist_ok=True)

    pdb = work / "sys.pdb"
    _write_text(
        pdb,
        """
CRYST1   20.000   20.000   20.000  90.00  90.00  90.00 P 1           1
ATOM      1  X   FRG A   1      10.000  10.000  10.000  1.00  0.00           C
END
""",
    )

    top = work / "sys.top"
    _write_text(
        top,
        """
[ defaults ]
1 2 yes 1.0 1.0

[ atomtypes ]
; name  at.num  mass   charge  ptype  sigma   epsilon
X       0       1.000  0.000   A      0.600   0.000

[ moleculetype ]
FRG  1

[ atoms ]
; nr  type  resnr  residue  atom  cgnr  charge  mass
1   X     1      FRG      X     1     0.000   1.000

[ system ]
CavityProbe

[ molecules ]
FRG 1
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
X       0       1.000  0.000   A      0.600   0.000
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
1   X     1      FRG      X     1     0.000   1.000
""",
    )

    def run_case(tag: str, probe_radius: float) -> float:
        case = work / tag
        case.mkdir(parents=True, exist_ok=True)
        inp = case / "run.inp"
        _write_inp(
            inp,
            f"""
version:gcmc_2.0
par:{par}
fragitp:{frag}
fragname:FRG
fragconc:55.0
fragmuex:0.0

pdb:{pdb}
top:{top}
box_size:20.0 20.0 20.0
cutoff:12.0
grid_dx:1.0
probe_radius:{probe_radius}
use_cavity_bias:yes
use_vdw_radius_for_grid:no
exclude_hydrogens_from_grid:no

temperature:300.0
mcsteps:20
nprint:20
moves_per_step:1
mc_move_prob:1.0 0.0 0.0 0.0
random_seed:99
""",
        )
        out_prefix = case / "out" / "gcmc"
        out_prefix.parent.mkdir(parents=True, exist_ok=True)
        accept_log = case / "out" / "acceptance.jsonl"
        return _run_cavity_case(
            gcmc_cpu=gcmc_cpu,
            work=case,
            inp=inp,
            out_prefix=out_prefix,
            accept_log=accept_log,
            timeout=90,
        )

    frac_small = run_case("probe_small", probe_radius=0.5)
    frac_large = run_case("probe_large", probe_radius=3.0)

    assert frac_large < frac_small - 0.05


def test_cavity_bias_exclude_hydrogens_changes_cavity_fraction(gcmc_cpu, temp_dir):
    work = Path(temp_dir) / "cavity_bias_exclude_h"
    work.mkdir(parents=True, exist_ok=True)

    pdb = work / "sys.pdb"
    _write_text(
        pdb,
        """
CRYST1   20.000   20.000   20.000  90.00  90.00  90.00 P 1           1
ATOM      1  O   WAT A   1      10.000  10.000  10.000  1.00  0.00           O
ATOM      2  H1  WAT A   1      10.957  10.000  10.000  1.00  0.00           H
ATOM      3  H2  WAT A   1       9.043  10.000  10.000  1.00  0.00           H
ATOM      4  O   WAT A   2      14.000  10.000  10.000  1.00  0.00           O
ATOM      5  H1  WAT A   2      14.957  10.000  10.000  1.00  0.00           H
ATOM      6  H2  WAT A   2      13.043  10.000  10.000  1.00  0.00           H
ATOM      7  O   WAT A   3       6.000  10.000  10.000  1.00  0.00           O
ATOM      8  H1  WAT A   3       6.957  10.000  10.000  1.00  0.00           H
ATOM      9  H2  WAT A   3       5.043  10.000  10.000  1.00  0.00           H
END
""",
    )

    top = work / "sys.top"
    _write_text(
        top,
        """
[ defaults ]
1 2 yes 1.0 1.0

[ atomtypes ]
O   8   15.9994  0.000  A  0.315  0.636
H   1   1.008    0.000  A  0.000  0.000

[ moleculetype ]
WAT  3

[ atoms ]
1  O   1  WAT  O   1  0.000  15.9994
2  H   1  WAT  H1  2  0.000  1.008
3  H   1  WAT  H2  3  0.000  1.008

[ system ]
CavityExcludeH

[ molecules ]
WAT 3
""",
    )

    par = work / "par.itp"
    _write_text(
        par,
        """
[ defaults ]
1 2 yes 1.0 1.0

[ atomtypes ]
O   8   15.9994  0.000  A  0.315  0.636
H   1   1.008    0.000  A  0.000  0.000
""",
    )

    frag = work / "frag.itp"
    _write_text(
        frag,
        """
[ moleculetype ]
WAT  3

[ atoms ]
1  O   1  WAT  O   1  0.000  15.9994
2  H   1  WAT  H1  2  0.000  1.008
3  H   1  WAT  H2  3  0.000  1.008
""",
    )

    def run_case(tag: str, exclude_h: str) -> float:
        case = work / tag
        case.mkdir(parents=True, exist_ok=True)
        inp = case / "run.inp"
        _write_inp(
        inp,
        f"""
version:gcmc_2.0
par:{par}
fragitp:{frag}
fragname:WAT
fragconc:55.0
fragmuex:0.0

pdb:{pdb}
top:{top}
box_size:20.0 20.0 20.0
cutoff:12.0
grid_dx:1.0
probe_radius:1.4
use_cavity_bias:yes
use_vdw_radius_for_grid:yes
exclude_hydrogens_from_grid:{exclude_h}

temperature:300.0
mcsteps:20
nprint:20
moves_per_step:1
mc_move_prob:1.0 0.0 0.0 0.0
random_seed:101
""",
        )
        out_prefix = case / "out" / "gcmc"
        out_prefix.parent.mkdir(parents=True, exist_ok=True)
        accept_log = case / "out" / "acceptance.jsonl"
        return _run_cavity_case(
            gcmc_cpu=gcmc_cpu,
            work=case,
            inp=inp,
            out_prefix=out_prefix,
            accept_log=accept_log,
            timeout=90,
        )

    frac_keep = run_case("keep_h", exclude_h="no")
    frac_exclude = run_case("exclude_h", exclude_h="yes")

    assert frac_exclude > frac_keep + 0.004

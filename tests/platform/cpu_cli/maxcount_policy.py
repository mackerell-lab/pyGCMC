"""
Regression tests for the hard capacity (`maxCount`) policy.

We use a gcmc_gpu-aligned safety bound:
  maxCount(type) = N_initial(type) + mcsteps + 200

This is treated as a hard memory/capacity limit, not a physical constraint.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from .inp_units_compat_helpers import _run_gcmc_cpu, _write_inp


def _write_text(path: Path, content: str) -> None:
    path.write_text(content.strip() + "\n")


def _read_fragment_max_count(top_path: Path, fragname: str) -> int:
    want = fragname.strip().upper()
    in_frag = False
    for line in top_path.read_text().splitlines():
        if line.strip() == "[ fragments ]":
            in_frag = True
            continue
        if not in_frag:
            continue
        if line.startswith("["):
            break
        if not line.strip() or line.lstrip().startswith(";"):
            continue
        cols = line.split()
        if cols and cols[0].strip().upper() == want:
            assert len(cols) >= 3, f"Unexpected [ fragments ] row: {line}"
            return int(cols[2])
    raise AssertionError(f"Fragment {fragname} not found in {top_path}")


def test_maxcount_is_initial_plus_mcsteps_plus_buffer(gcmc_cpu, temp_dir):
    work = Path(temp_dir) / "maxcount_policy"
    work.mkdir(parents=True, exist_ok=True)

    pdb = work / "sys.pdb"
    _write_text(
        pdb,
        """
CRYST1   10.000   10.000   10.000  90.00  90.00  90.00 P 1           1
ATOM      1  X   FRG A   1       5.000   5.000   5.000  1.00  0.00           C
END
""",
    )

    top = work / "sys.top"
    _write_text(
        top,
        """
[ defaults ]
1 2 yes 0.5 0.8333

[ atomtypes ]
; name  at.num  mass   charge  ptype  sigma   epsilon
X       0       1.000  0.000   A      0.000   0.000

[ moleculetype ]
FRG  1

[ atoms ]
; nr  type  resnr  residue  atom  cgnr  charge  mass
1   X     1      FRG      X     1     0.000   1.000

[ system ]
MaxCountPolicy

[ molecules ]
FRG 1
""",
    )

    par = work / "par.itp"
    _write_text(
        par,
        """
[ defaults ]
1 2 yes 0.5 0.8333

[ atomtypes ]
; name  at.num  mass   charge  ptype  sigma   epsilon
X       0       1.000  0.000   A      0.000   0.000
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

    mcsteps = 10
    out_prefix = work / "out" / "gcmc"
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    inp = work / "run.inp"
    _write_inp(
        inp,
        f"""
random_seed:123
par:{par}
fragitp:{frag}
fragname:FRG
fragconc:55.0
fragmuex:0.0

pdb:{pdb}
top:{top}
box_size:10.0 10.0 10.0
cutoff:4.0
temperature:300.0
moves_per_step:1
mcsteps:{mcsteps}
nprint:1
mc_move_prob:0 0 1 0
""",
    )

    result = _run_gcmc_cpu(
        gcmc_cpu,
        workdir=work,
        inp=inp,
        out_prefix=out_prefix,
        timeout=30,
    )
    assert result.returncode == 0, result.stdout + result.stderr

    out_top = Path(f"{out_prefix}_final.top")
    assert out_top.exists()
    out_pdb = Path(f"{out_prefix}_final.pdb")
    assert out_pdb.exists()
    resnames = set()
    for line in out_pdb.read_text().splitlines():
        if line.startswith(("ATOM", "HETATM")):
            resnames.add(line[17:20].strip().upper())
    assert "FRG" in resnames, f"Initial system did not survive into output PDB: {sorted(resnames)}"
    assert _read_fragment_max_count(out_top, "FRG") == (1 + mcsteps + 200)

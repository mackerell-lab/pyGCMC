"""
End-to-end validation that fragment templates load real geometry from monomerdir/<frag>.pdb.

Historically FragmentLibrary defaulted all coordinates to 0, which degenerates insertions.
This test ensures inserted molecules have non-degenerate internal distances.
"""

from __future__ import annotations

import math
import subprocess
from pathlib import Path

import pytest


def _parse_first_residue_coords(pdb_path: Path) -> list[tuple[str, float, float, float]]:
    atoms: list[tuple[str, float, float, float]] = []
    first_resid = None

    for line in pdb_path.read_text().splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        atom_name = line[12:16].strip()
        resid = int(line[22:26])
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])

        if first_resid is None:
            first_resid = resid
        if resid != first_resid:
            break

        atoms.append((atom_name, x, y, z))

    return atoms


def _dist(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    dx = a[0] - b[0]
    dy = a[1] - b[1]
    dz = a[2] - b[2]
    return math.sqrt(dx * dx + dy * dy + dz * dz)


def test_fragment_template_geometry_from_monomerdir(gcmc_cpu, test_data_dir, temp_dir):
    work = Path(temp_dir) / "frag_geom"
    work.mkdir(parents=True, exist_ok=True)

    mol_dir = test_data_dir / "charmm36.ff" / "mol"
    itp = mol_dir / "sol.itp"
    assert itp.exists()

    out_prefix = work / "out" / "gcmc"
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    inp = work / "test.inp"
    inp.write_text(
        f"""
inp_units:gcmc_gpu
monomerdir:{mol_dir}
fragitp:{itp}
box_size:30.0 30.0 30.0
cutoff:12.0
temperature:300.0
moves_per_step:1
mcsteps:1
nprint:1
fragname:SOL
fragmuex:5.0
attempt_prob_ins:1.0
attempt_prob_del:0.0
attempt_prob_trn:0.0
attempt_prob_rot:0.0
""".strip()
        + "\n"
    )

    result = subprocess.run(
        [gcmc_cpu, "--inp", str(inp), "--prefix", str(out_prefix)],
        cwd=str(work),
        capture_output=True,
        text=True,
        timeout=20,
    )
    assert result.returncode == 0, result.stdout + result.stderr

    out_pdb = Path(f"{out_prefix}_final.pdb")
    assert out_pdb.exists()

    atoms = _parse_first_residue_coords(out_pdb)
    assert len(atoms) >= 3, f"Expected >=3 atoms in first residue, got {len(atoms)}"

    # Use first 3 atoms (TIP3P-like water): OW, HW1, HW2. Distances are in Angstrom in PDB.
    a0 = (atoms[0][1], atoms[0][2], atoms[0][3])
    a1 = (atoms[1][1], atoms[1][2], atoms[1][3])
    a2 = (atoms[2][1], atoms[2][2], atoms[2][3])

    d01 = _dist(a0, a1)
    d02 = _dist(a0, a2)
    d12 = _dist(a1, a2)

    # Expected geometry from the template PDB (≈0.957 Å bonds, ≈1.514 Å H-H)
    assert d01 == pytest.approx(0.957, abs=0.05)
    assert d02 == pytest.approx(0.957, abs=0.05)
    assert d12 == pytest.approx(1.514, abs=0.08)


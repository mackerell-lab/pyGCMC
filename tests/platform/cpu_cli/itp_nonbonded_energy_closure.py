"""
End-to-end regression: insertion -> deletion symmetry for ITP nonbonded LJ parsing.

This is a "non-cheating" check because the expected deltaU is reconstructed from:
  - the written intermediate PDB coordinates (Å), and
  - known sigma/epsilon values from a minimal [atomtypes] PAR file (nm / kJ/mol),
without relying on stdout/stderr log text.
"""

from __future__ import annotations

import json
import math
import subprocess
from pathlib import Path

import pytest


def _write_text(path: Path, content: str) -> None:
    path.write_text(content.strip() + "\n")


def _first_accept_record(path: Path, *, move: str, species: str) -> dict:
    records = [json.loads(line) for line in path.read_text().splitlines() if line.strip()]
    want_move = move.strip().lower()
    want_species = species.strip().upper()
    for rec in records:
        if str(rec.get("move", "")).strip().lower() != want_move:
            continue
        if str(rec.get("species", "")).strip().upper() != want_species:
            continue
        return rec
    raise AssertionError(f"No {move}/{species} record found in {path}; first records: {records[:3]}")


def _first_atom_xyz_angstrom(pdb_path: Path, *, resname: str) -> tuple[float, float, float]:
    want = resname.strip().upper()
    for line in pdb_path.read_text().splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        if line[17:20].strip().upper() != want:
            continue
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        return x, y, z
    raise AssertionError(f"Residue {resname} not found in {pdb_path}")


def _count_residues_in_pdb(pdb_path: Path, *, resname: str) -> int:
    want = resname.strip().upper()
    return sum(
        1
        for line in pdb_path.read_text().splitlines()
        if line.startswith(("ATOM", "HETATM")) and line[17:20].strip().upper() == want
    )


def _lj_energy_kj_mol(*, r_nm: float, sigma_nm: float, eps_kj_mol: float) -> float:
    r2 = r_nm * r_nm
    sigma_r2 = (sigma_nm * sigma_nm) / r2
    sigma_r6 = sigma_r2 * sigma_r2 * sigma_r2
    sigma_r12 = sigma_r6 * sigma_r6
    return 4.0 * eps_kj_mol * (sigma_r12 - sigma_r6)


def test_itp_atomtypes_insertion_then_deletion_deltaU_matches_analytic_lj_pair_energy(
    gcmc_cpu, temp_dir
):
    """
    Scenario:
    - Initial system: one neutral C atom (MOL) at box center.
    - Fragment: one neutral NA atom.
    - PAR (ITP): defines sigma/epsilon for both types.

    Steps:
    1) Run 1 step of insertion-only (accepted by construction).
       - Check deltaU_ins matches analytic LJ(C-NA) from intermediate PDB coords.
    2) Use the produced *_final.{pdb,top} as initial system and run 1 step of deletion-only.
       - Check deltaU_del == -deltaU_ins == -E_LJ(C-NA), and NA disappears from final PDB.
    """
    work = Path(temp_dir) / "itp_nonbonded_energy" / "ins_del_closure"
    work.mkdir(parents=True, exist_ok=True)

    pdb0 = work / "sys.pdb"
    _write_text(
        pdb0,
        """
CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1           1
ATOM      1  C   MOL A   1      15.000  15.000  15.000  1.00  0.00           C
END
""",
    )

    top0 = work / "sys.top"
    _write_text(
        top0,
        """
[ defaults ]
1 2 yes 0.5 0.8333

[ moleculetype ]
MOL  2

[ atoms ]
; nr  type  resnr  residue  atom  cgnr  charge    mass
1   C     1      MOL      C     1     0.000   12.011

[ system ]
Minimal

[ molecules ]
MOL  1
""",
    )

    frag_itp = work / "na.itp"
    _write_text(
        frag_itp,
        """
[ moleculetype ]
NA  1

[ atoms ]
; nr  type  resnr  residue  atom  cgnr  charge  mass
1   NA    1      NA       NA    1     0.000   22.9898
""",
    )

    # Sigma/epsilon are nm / kJ/mol (GROMACS convention)
    sigma_c = 0.300
    eps_c = 1.000
    sigma_na = 0.400
    eps_na = 2.000

    par = work / "ffnonbonded.itp"
    _write_text(
        par,
        f"""
[ atomtypes ]
; name  at.num  mass     charge   ptype    sigma      epsilon
C       6       12.011   0.000    A        {sigma_c:.6f}   {eps_c:.6f}
NA      11      22.990   0.000    A        {sigma_na:.6f}  {eps_na:.6f}
""",
    )

    out_prefix_ins = work / "out_ins" / "gcmc"
    out_prefix_ins.parent.mkdir(parents=True, exist_ok=True)
    accept_log_ins = work / "out_ins" / "acceptance.jsonl"

    inp_ins = work / "insert.inp"
    _write_text(
        inp_ins,
        f"""
random_seed:123
par:{par}
top:{top0}
pdb:{pdb0}
fragitp:{frag_itp}

box_size:30.0 30.0 30.0
cutoff:12.0
gcmc_region:box 20.0 12.5 12.5 25.0 17.5 17.5

temperature:300.0
moves_per_step:1
mcsteps:1
nprint:1

fragname:NA
fragconc:0.01
fragmuex:10.0
mc_move_prob:1 0 0 0
""",
    )

    result_ins = subprocess.run(
        [gcmc_cpu, "--inp", str(inp_ins), "--prefix", str(out_prefix_ins), "--dump-accept", str(accept_log_ins)],
        cwd=str(work),
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert result_ins.returncode == 0, result_ins.stdout + result_ins.stderr

    rec_ins = _first_accept_record(accept_log_ins, move="insertion", species="NA")
    assert bool(rec_ins.get("accepted")) is True

    final_pdb_ins = Path(f"{out_prefix_ins}_final.pdb")
    final_top_ins = Path(f"{out_prefix_ins}_final.top")
    assert final_pdb_ins.exists()
    assert final_top_ins.exists()

    # Compute expected LJ(C-NA) from intermediate geometry
    x0, y0, z0 = _first_atom_xyz_angstrom(final_pdb_ins, resname="MOL")
    x1, y1, z1 = _first_atom_xyz_angstrom(final_pdb_ins, resname="NA")

    dx = (x1 - x0) / 10.0
    dy = (y1 - y0) / 10.0
    dz = (z1 - z0) / 10.0
    r_nm = math.sqrt(dx * dx + dy * dy + dz * dz)

    sigma_mixed = 0.5 * (sigma_c + sigma_na)
    eps_mixed = math.sqrt(eps_c * eps_na)
    expected_pair = _lj_energy_kj_mol(r_nm=r_nm, sigma_nm=sigma_mixed, eps_kj_mol=eps_mixed)

    assert float(rec_ins["deltaU"]) == pytest.approx(expected_pair, rel=5e-4, abs=5e-4)

    # Build a minimal *valid* topology for the deletion run.
    # Note: gcmc_cpu's generated *_final.top is a summary and is not parseable as a full GROMACS topology.
    top_with_na = work / "sys_with_na.top"
    _write_text(
        top_with_na,
        """
[ defaults ]
1 2 yes 0.5 0.8333

[ moleculetype ]
MOL  2

[ atoms ]
; nr  type  resnr  residue  atom  cgnr  charge    mass
1   C     1      MOL      C     1     0.000   12.011

[ moleculetype ]
NA  1

[ atoms ]
; nr  type  resnr  residue  atom  cgnr  charge  mass
1   NA    1      NA       NA    1     0.000   22.9898

[ system ]
Minimal

[ molecules ]
MOL  1
NA   1
""",
    )

    out_prefix_del = work / "out_del" / "gcmc"
    out_prefix_del.parent.mkdir(parents=True, exist_ok=True)
    accept_log_del = work / "out_del" / "acceptance.jsonl"

    inp_del = work / "delete.inp"
    _write_text(
        inp_del,
        f"""
random_seed:456
par:{par}
top:{top_with_na}
pdb:{final_pdb_ins}
fragitp:{frag_itp}

box_size:30.0 30.0 30.0
cutoff:12.0
gcmc_region:box 20.0 12.5 12.5 25.0 17.5 17.5

temperature:300.0
moves_per_step:1
mcsteps:1
nprint:1

fragname:NA
fragconc:0.01
fragmuex:-10.0
mc_move_prob:0 1 0 0
""",
    )

    result_del = subprocess.run(
        [gcmc_cpu, "--inp", str(inp_del), "--prefix", str(out_prefix_del), "--dump-accept", str(accept_log_del)],
        cwd=str(work),
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert result_del.returncode == 0, result_del.stdout + result_del.stderr

    rec_del = _first_accept_record(accept_log_del, move="deletion", species="NA")
    assert bool(rec_del.get("accepted")) is True
    assert float(rec_del["deltaU"]) == pytest.approx(-expected_pair, rel=5e-4, abs=5e-4)

    final_pdb_del = Path(f"{out_prefix_del}_final.pdb")
    assert final_pdb_del.exists()
    assert _count_residues_in_pdb(final_pdb_del, resname="NA") == 0
    assert _count_residues_in_pdb(final_pdb_del, resname="MOL") == 1

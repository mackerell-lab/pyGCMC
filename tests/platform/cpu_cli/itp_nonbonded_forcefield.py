"""
End-to-end tests for GROMACS ITP nonbonded parsing -> MCState LJ matrix wiring.

These tests provide an independent ("truth") check by reconstructing the expected
Lennard-Jones energy from:
  - the written final PDB coordinates (Å), and
  - known sigma/epsilon values from a minimal .itp PAR file (nm / kJ/mol),
then comparing against the structured acceptance log deltaU (kJ/mol).

This catches regressions where the simulation silently falls back to a placeholder/zero LJ matrix.
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


def _lj_energy_kj_mol(*, r_nm: float, sigma_nm: float, eps_kj_mol: float) -> float:
    r2 = r_nm * r_nm
    sigma_r2 = (sigma_nm * sigma_nm) / r2
    sigma_r6 = sigma_r2 * sigma_r2 * sigma_r2
    sigma_r12 = sigma_r6 * sigma_r6
    return 4.0 * eps_kj_mol * (sigma_r12 - sigma_r6)


def test_itp_atomtypes_builds_lj_matrix_and_deltaU_matches_analytic_lj(
    gcmc_cpu, temp_dir
):
    """
    Independent energy check:
    - Initial system: one neutral atom type C at box center
    - Insert: one neutral atom type NA within a region near the initial atom
    - PAR: defines sigma/epsilon for both types (nm/kJ)
    Expect: acceptance log deltaU equals analytic LJ(C-NA) computed from final PDB coords.
    """
    work = Path(temp_dir) / "itp_nonbonded_energy" / "atomtypes_only"
    work.mkdir(parents=True, exist_ok=True)

    pdb = work / "sys.pdb"
    _write_text(
        pdb,
        """
CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1           1
ATOM      1  C   MOL A   1      15.000  15.000  15.000  1.00  0.00           C
END
""",
    )

    top = work / "sys.top"
    _write_text(
        top,
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

    out_prefix = work / "out" / "gcmc"
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    accept_log = work / "out" / "acceptance.jsonl"

    inp = work / "run.inp"
    _write_text(
        inp,
        f"""
inp_units:gcmc_gpu
random_seed:123
par:{par}
top:{top}
pdb:{pdb}
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

    result = subprocess.run(
        [gcmc_cpu, "--inp", str(inp), "--prefix", str(out_prefix), "--dump-accept", str(accept_log)],
        cwd=str(work),
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert result.returncode == 0, result.stdout + result.stderr

    rec = _first_accept_record(accept_log, move="insertion", species="NA")
    assert bool(rec.get("accepted")) is True

    final_pdb = Path(f"{out_prefix}_final.pdb")
    assert final_pdb.exists()

    x0, y0, z0 = _first_atom_xyz_angstrom(final_pdb, resname="MOL")
    x1, y1, z1 = _first_atom_xyz_angstrom(final_pdb, resname="NA")

    dx = (x1 - x0) / 10.0
    dy = (y1 - y0) / 10.0
    dz = (z1 - z0) / 10.0
    r_nm = math.sqrt(dx * dx + dy * dy + dz * dz)

    sigma_mixed = 0.5 * (sigma_c + sigma_na)
    eps_mixed = math.sqrt(eps_c * eps_na)
    expected = _lj_energy_kj_mol(r_nm=r_nm, sigma_nm=sigma_mixed, eps_kj_mol=eps_mixed)

    assert float(rec["deltaU"]) == pytest.approx(expected, rel=5e-4, abs=5e-4)


def test_itp_nonbond_params_override_applies_nbfix_in_lj_matrix(
    gcmc_cpu, temp_dir
):
    """
    End-to-end check that [ nonbond_params ] entries are applied as NBFIX overrides.

    We set an explicit C-NA override (sigma/epsilon) and verify deltaU matches the
    override parameters (not the Lorentz-Berthelot mixed parameters).
    """
    work = Path(temp_dir) / "itp_nonbonded_energy" / "nonbond_params_override"
    work.mkdir(parents=True, exist_ok=True)

    pdb = work / "sys.pdb"
    _write_text(
        pdb,
        """
CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1           1
ATOM      1  C   MOL A   1      15.000  15.000  15.000  1.00  0.00           C
END
""",
    )

    top = work / "sys.top"
    _write_text(
        top,
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

    sigma_c = 0.300
    eps_c = 1.000
    sigma_na = 0.400
    eps_na = 2.000

    sigma_override = 0.400
    eps_override = 20.000

    par = work / "ffnonbonded.itp"
    _write_text(
        par,
        f"""
[ atomtypes ]
; name  at.num  mass     charge   ptype    sigma      epsilon
C       6       12.011   0.000    A        {sigma_c:.6f}   {eps_c:.6f}
NA      11      22.990   0.000    A        {sigma_na:.6f}  {eps_na:.6f}

[ nonbond_params ]
; type1  type2  func  sigma  epsilon
C   NA   1   {sigma_override:.6f}  {eps_override:.6f}
""",
    )

    out_prefix = work / "out" / "gcmc"
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    accept_log = work / "out" / "acceptance.jsonl"

    inp = work / "run.inp"
    _write_text(
        inp,
        f"""
inp_units:gcmc_gpu
random_seed:123
par:{par}
top:{top}
pdb:{pdb}
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

    result = subprocess.run(
        [gcmc_cpu, "--inp", str(inp), "--prefix", str(out_prefix), "--dump-accept", str(accept_log)],
        cwd=str(work),
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert result.returncode == 0, result.stdout + result.stderr

    rec = _first_accept_record(accept_log, move="insertion", species="NA")
    assert bool(rec.get("accepted")) is True

    final_pdb = Path(f"{out_prefix}_final.pdb")
    assert final_pdb.exists()

    x0, y0, z0 = _first_atom_xyz_angstrom(final_pdb, resname="MOL")
    x1, y1, z1 = _first_atom_xyz_angstrom(final_pdb, resname="NA")

    dx = (x1 - x0) / 10.0
    dy = (y1 - y0) / 10.0
    dz = (z1 - z0) / 10.0
    r_nm = math.sqrt(dx * dx + dy * dy + dz * dz)

    expected = _lj_energy_kj_mol(r_nm=r_nm, sigma_nm=sigma_override, eps_kj_mol=eps_override)
    assert float(rec["deltaU"]) == pytest.approx(expected, rel=5e-4, abs=5e-4)


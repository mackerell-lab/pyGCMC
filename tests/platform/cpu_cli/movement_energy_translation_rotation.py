"""
End-to-end (CLI) movement+energy closure tests for TRANSLATE/ROTATE.

Goal: prevent "silently runs but wrong" regressions by locking down:
- translation/rotation moves are actually attempted and exported in --dump-accept JSONL
- logged deltaU matches an independent analytic LJ energy difference computed from PDB geometry
- logged pAcc matches Metropolis: min(1, exp(-beta * deltaU) * bias) for these moves

Notes:
- We deliberately set an astronomically high temperature to make Metropolis acceptance ~1.0
  (and numerically round to 1.0 in double precision), so the final PDB reflects the proposed move.
  This keeps the test deterministic while still exercising real energy calculations.
- All length-like inputs/outputs use Å (gcmc_gpu/opencl convention); internal is nm.
"""

from __future__ import annotations

import json
import math
import subprocess
from pathlib import Path

import pytest


def _write_text(path: Path, content: str) -> None:
    path.write_text(content.strip() + "\n")


def _read_cryst1_box_angstrom(pdb_path: Path) -> tuple[float, float, float]:
    for line in pdb_path.read_text().splitlines():
        if line.startswith("CRYST1"):
            parts = line.split()
            assert len(parts) >= 4, f"Unexpected CRYST1 format: {line}"
            return float(parts[1]), float(parts[2]), float(parts[3])
    raise AssertionError(f"CRYST1 not found in {pdb_path}")


def _first_accept_record(path: Path, *, move: str, species: str) -> dict:
    want_move = move.strip().lower()
    want_species = species.strip().upper()
    records = [json.loads(line) for line in path.read_text().splitlines() if line.strip()]
    for rec in records:
        if str(rec.get("move", "")).strip().lower() != want_move:
            continue
        if str(rec.get("species", "")).strip().upper() != want_species:
            continue
        return rec
    raise AssertionError(f"No {move}/{species} record found in {path}; first records: {records[:3]}")


def _atom_xyz_angstrom(pdb_path: Path, *, resname: str, atom: str) -> tuple[float, float, float]:
    want_res = resname.strip().upper()
    want_atom = atom.strip().upper()
    for line in pdb_path.read_text().splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        if line[17:20].strip().upper() != want_res:
            continue
        if line[12:16].strip().upper() != want_atom:
            continue
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        return x, y, z
    raise AssertionError(f"Atom {want_atom}/{want_res} not found in {pdb_path}")


def _residue_atoms_xyz_angstrom(pdb_path: Path, *, resname: str) -> dict[str, tuple[float, float, float]]:
    want_res = resname.strip().upper()
    out: dict[str, tuple[float, float, float]] = {}
    for line in pdb_path.read_text().splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        if line[17:20].strip().upper() != want_res:
            continue
        atom = line[12:16].strip().upper()
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        out[atom] = (x, y, z)
    if not out:
        raise AssertionError(f"Residue {want_res} not found in {pdb_path}")
    return out


def _min_image_delta_nm(
    a_ang: tuple[float, float, float],
    b_ang: tuple[float, float, float],
    *,
    box_ang: tuple[float, float, float],
) -> tuple[float, float, float]:
    ax, ay, az = a_ang
    bx, by, bz = b_ang
    dx = (bx - ax) / 10.0
    dy = (by - ay) / 10.0
    dz = (bz - az) / 10.0
    box_nm = (box_ang[0] / 10.0, box_ang[1] / 10.0, box_ang[2] / 10.0)
    dx -= box_nm[0] * round(dx / box_nm[0])
    dy -= box_nm[1] * round(dy / box_nm[1])
    dz -= box_nm[2] * round(dz / box_nm[2])
    return dx, dy, dz


def _lj_energy_kj_mol(*, r_nm: float, sigma_nm: float, eps_kj_mol: float) -> float:
    r2 = r_nm * r_nm
    sigma_r2 = (sigma_nm * sigma_nm) / r2
    sigma_r6 = sigma_r2 * sigma_r2 * sigma_r2
    sigma_r12 = sigma_r6 * sigma_r6
    return 4.0 * eps_kj_mol * (sigma_r12 - sigma_r6)


def _metropolis_pacc(*, beta: float, delta_u: float, bias: float = 1.0) -> float:
    return min(1.0, math.exp(-beta * delta_u) * bias)


def test_translation_move_deltaU_matches_analytic_lj_and_is_logged_in_dump_accept(gcmc_cpu, temp_dir):
    work = Path(temp_dir) / "move_energy" / "translation"
    work.mkdir(parents=True, exist_ok=True)

    # Simple 2-residue system:
    # - MOL: one atom type C (static receptor)
    # - SOD: one atom type NA (active fragment, will be translated)
    pdb = work / "sys.pdb"
    _write_text(
        pdb,
        """
CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1           1
ATOM      1  C   MOL A   1      15.000  15.000  15.000  1.00  0.00           C
ATOM      2  NA  SOD B   1      22.000  15.000  15.000  1.00  0.00          NA
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

    [ moleculetype ]
    SOD  1
    
    [ atoms ]
    ; nr  type  resnr  residue  atom  cgnr  charge  mass
    1   NA    1      SOD      NA    1     0.000   22.9898

[ system ]
Minimal

    [ molecules ]
    MOL  1
    SOD  1
""",
    )

    frag_itp = work / "na.itp"
    _write_text(
        frag_itp,
        """
[ moleculetype ]
SOD  1

[ atoms ]
; nr  type  resnr  residue  atom  cgnr  charge  mass
1   NA    1      SOD      NA    1     0.000   22.9898
""",
    )

    # Template coordinates (next to ITP, as SimulationInputBuilder fallback expects)
    frag_pdb = work / "na.pdb"
    _write_text(
        frag_pdb,
        """
ATOM      1  NA  SOD B   1       0.000   0.000   0.000  1.00  0.00          NA
END
""",
    )

    sigma_c = 0.300
    eps_c = 1.000
    sigma_na = 0.400
    eps_na = 2.000
    sigma_mixed = 0.5 * (sigma_c + sigma_na)
    eps_mixed = math.sqrt(eps_c * eps_na)

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

    # Make acceptance deterministic: beta ~ 0 => pAcc rounds to 1.0
    huge_temperature = 1.0e30

    inp = work / "run.inp"
    _write_text(
        inp,
        f"""
random_seed:123
par:{par}
top:{top}
pdb:{pdb}
fragitp:{frag_itp}

box_size:30.0 30.0 30.0
cutoff:12.0

temperature:{huge_temperature}
moves_per_step:1
mcsteps:1
nprint:1

fragname:SOD
fragconc:0.01
fragmuex:0.0
mc_move_prob:0 0 1 0
""",
    )

    result = subprocess.run(
        [gcmc_cpu, "--inp", str(inp), "--prefix", str(out_prefix), "--dump-accept", str(accept_log)],
        cwd=str(work),
        capture_output=True,
        text=True,
        timeout=60,
    )
    assert result.returncode == 0, result.stdout + result.stderr

    rec = _first_accept_record(accept_log, move="translation", species="SOD")
    assert bool(rec.get("accepted")) is True

    final_pdb = Path(f"{out_prefix}_final.pdb")
    assert final_pdb.exists()

    box_ang = _read_cryst1_box_angstrom(final_pdb)
    c0 = _atom_xyz_angstrom(pdb, resname="MOL", atom="C")
    na0 = _atom_xyz_angstrom(pdb, resname="SOD", atom="NA")
    c1 = _atom_xyz_angstrom(final_pdb, resname="MOL", atom="C")
    na1 = _atom_xyz_angstrom(final_pdb, resname="SOD", atom="NA")
    assert c1 == pytest.approx(c0, abs=1e-6)

    dx0, dy0, dz0 = _min_image_delta_nm(c0, na0, box_ang=box_ang)
    dx1, dy1, dz1 = _min_image_delta_nm(c1, na1, box_ang=box_ang)
    r0 = math.sqrt(dx0 * dx0 + dy0 * dy0 + dz0 * dz0)
    r1 = math.sqrt(dx1 * dx1 + dy1 * dy1 + dz1 * dz1)

    e0 = _lj_energy_kj_mol(r_nm=r0, sigma_nm=sigma_mixed, eps_kj_mol=eps_mixed)
    e1 = _lj_energy_kj_mol(r_nm=r1, sigma_nm=sigma_mixed, eps_kj_mol=eps_mixed)
    expected_delta_u = e1 - e0

    assert float(rec["deltaU"]) == pytest.approx(expected_delta_u, rel=5e-4, abs=5e-4)

    beta = float(rec["beta"])
    expected_pacc = _metropolis_pacc(beta=beta, delta_u=float(rec["deltaU"]), bias=float(rec["bias"]))
    assert float(rec["pAcc"]) == pytest.approx(expected_pacc, rel=1e-12, abs=1e-12)


def test_rotation_move_deltaU_matches_analytic_lj_for_multi_atom_fragment(gcmc_cpu, temp_dir):
    work = Path(temp_dir) / "move_energy" / "rotation"
    work.mkdir(parents=True, exist_ok=True)

    # 2-residue system:
    # - MOL: one atom type C (static receptor)
    # - FRG: two atoms A/B; rotation changes their distances to C -> deltaU != 0
    pdb = work / "sys.pdb"
    _write_text(
        pdb,
        """
CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1           1
ATOM      1  C   MOL A   1      15.000  15.000  15.000  1.00  0.00           C
ATOM      2  A   FRG B   1      22.000  15.000  15.000  1.00  0.00           A
ATOM      3  B   FRG B   1      23.000  15.000  15.000  1.00  0.00           B
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

[ moleculetype ]
FRG  1

[ atoms ]
; nr  type  resnr  residue  atom  cgnr  charge  mass
1   A     1      FRG      A     1     0.000   10.000
2   B     1      FRG      B     1     0.000   10.000

[ system ]
Minimal

[ molecules ]
MOL  1
FRG  1
""",
    )

    frag_itp = work / "frg.itp"
    _write_text(
        frag_itp,
        """
[ moleculetype ]
FRG  1

[ atoms ]
; nr  type  resnr  residue  atom  cgnr  charge  mass
1   A     1      FRG      A     1     0.000   10.000
2   B     1      FRG      B     1     0.000   10.000
""",
    )

    # Template coordinates (next to ITP). Use a non-collinear geometry to make rotation measurable.
    frag_pdb = work / "frg.pdb"
    _write_text(
        frag_pdb,
        """
ATOM      1  A   FRG B   1       0.000   0.000   0.000  1.00  0.00           A
ATOM      2  B   FRG B   1       1.000   0.000   0.000  1.00  0.00           B
END
""",
    )

    sigma_c = 0.300
    eps_c = 1.000
    sigma_a = 0.350
    eps_a = 1.500
    sigma_b = 0.450
    eps_b = 0.800

    sigma_ca = 0.5 * (sigma_c + sigma_a)
    eps_ca = math.sqrt(eps_c * eps_a)
    sigma_cb = 0.5 * (sigma_c + sigma_b)
    eps_cb = math.sqrt(eps_c * eps_b)

    par = work / "ffnonbonded.itp"
    _write_text(
        par,
        f"""
[ atomtypes ]
; name  at.num  mass     charge   ptype    sigma      epsilon
C       6       12.011   0.000    A        {sigma_c:.6f}   {eps_c:.6f}
A       1       10.000   0.000    A        {sigma_a:.6f}   {eps_a:.6f}
B       2       10.000   0.000    A        {sigma_b:.6f}   {eps_b:.6f}
""",
    )

    out_prefix = work / "out" / "gcmc"
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    accept_log = work / "out" / "acceptance.jsonl"

    huge_temperature = 1.0e30

    inp = work / "run.inp"
    _write_text(
        inp,
        f"""
random_seed:321
par:{par}
top:{top}
pdb:{pdb}
fragitp:{frag_itp}

box_size:30.0 30.0 30.0
cutoff:12.0

temperature:{huge_temperature}
moves_per_step:1
mcsteps:1
nprint:1

fragname:FRG
fragconc:0.01
fragmuex:0.0
mc_move_prob:0 0 0 1
""",
    )

    result = subprocess.run(
        [gcmc_cpu, "--inp", str(inp), "--prefix", str(out_prefix), "--dump-accept", str(accept_log)],
        cwd=str(work),
        capture_output=True,
        text=True,
        timeout=60,
    )
    assert result.returncode == 0, result.stdout + result.stderr

    rec = _first_accept_record(accept_log, move="rotation", species="FRG")
    assert bool(rec.get("accepted")) is True

    final_pdb = Path(f"{out_prefix}_final.pdb")
    assert final_pdb.exists()

    box_ang = _read_cryst1_box_angstrom(final_pdb)
    c0 = _atom_xyz_angstrom(pdb, resname="MOL", atom="C")
    frg0 = _residue_atoms_xyz_angstrom(pdb, resname="FRG")
    c1 = _atom_xyz_angstrom(final_pdb, resname="MOL", atom="C")
    frg1 = _residue_atoms_xyz_angstrom(final_pdb, resname="FRG")
    assert c1 == pytest.approx(c0, abs=1e-6)
    assert set(frg0.keys()) == {"A", "B"}
    assert set(frg1.keys()) == {"A", "B"}

    def energy_from_geom(frg_atoms: dict[str, tuple[float, float, float]]) -> float:
        dx_a, dy_a, dz_a = _min_image_delta_nm(c0, frg_atoms["A"], box_ang=box_ang)
        dx_b, dy_b, dz_b = _min_image_delta_nm(c0, frg_atoms["B"], box_ang=box_ang)
        r_a = math.sqrt(dx_a * dx_a + dy_a * dy_a + dz_a * dz_a)
        r_b = math.sqrt(dx_b * dx_b + dy_b * dy_b + dz_b * dz_b)
        return _lj_energy_kj_mol(r_nm=r_a, sigma_nm=sigma_ca, eps_kj_mol=eps_ca) + _lj_energy_kj_mol(
            r_nm=r_b, sigma_nm=sigma_cb, eps_kj_mol=eps_cb
        )

    e0 = energy_from_geom(frg0)
    e1 = energy_from_geom(frg1)
    expected_delta_u = e1 - e0

    assert float(rec["deltaU"]) == pytest.approx(expected_delta_u, rel=5e-4, abs=5e-4)

    beta = float(rec["beta"])
    expected_pacc = _metropolis_pacc(beta=beta, delta_u=float(rec["deltaU"]), bias=float(rec["bias"]))
    assert float(rec["pAcc"]) == pytest.approx(expected_pacc, rel=1e-12, abs=1e-12)

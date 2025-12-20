"""
Physical-contract CLI tests for DIRECT+cutoff GCMC.

These are end-to-end checks that:
- acceptance probabilities follow the exact μVT formulas (ideal gas limit)
- number fluctuations follow a Poisson law in the ideal gas limit
- cutoff enforces zero interaction beyond the cutoff distance

All assertions are file-driven via --dump-accept JSONL and output PDBs.
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import pytest

from .inp_units_compat_helpers import _run_gcmc_cpu, _write_inp

COULOMB = 138.935458  # kJ·nm/mol/e^2


def _write_text(path: Path, content: str) -> None:
    path.write_text(content.strip() + "\n")


def _ensure_pygcmc_on_path() -> None:
    repo_root = Path(__file__).resolve().parents[3]
    bindings = repo_root / "build" / "modules" / "bindings"
    if bindings.exists():
        bindings_str = str(bindings)
        if bindings_str not in sys.path:
            sys.path.insert(0, bindings_str)


def _load_accept_records(path: Path) -> list[dict]:
    return [json.loads(line) for line in path.read_text().splitlines() if line.strip()]


def _first_accept_record(path: Path, *, move: str, species: str) -> dict:
    want_move = move.strip().lower()
    want_species = species.strip().upper()
    for rec in _load_accept_records(path):
        if str(rec.get("move", "")).strip().lower() != want_move:
            continue
        if str(rec.get("species", "")).strip().upper() != want_species:
            continue
        return rec
    raise AssertionError(f"No {move}/{species} record found in {path}")


def _read_cryst1_box_angstrom(pdb_path: Path) -> tuple[float, float, float]:
    for line in pdb_path.read_text().splitlines():
        if line.startswith("CRYST1"):
            parts = line.split()
            assert len(parts) >= 4, f"Unexpected CRYST1 format: {line}"
            return float(parts[1]), float(parts[2]), float(parts[3])
    raise AssertionError(f"CRYST1 not found in {pdb_path}")


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


def _coulomb_energy_kj_mol(*, r_nm: float, q1: float, q2: float) -> float:
    return COULOMB * q1 * q2 / r_nm


def _write_ideal_gas_itp(work: Path) -> tuple[Path, Path]:
    par = work / "ideal_par.itp"
    frag = work / "ideal_frag.itp"
    _write_text(
        par,
        """
[ defaults ]
1 2 yes 0.5 0.8333

[ atomtypes ]
; name  at.num  mass   charge  ptype  sigma   epsilon
X       0       1.000  0.000   A      0.300   0.000
""",
    )
    _write_text(
        frag,
        """
[ moleculetype ]
X   1

[ atoms ]
; nr  type  resnr  residue  atom  cgnr  charge  mass
1   X     1      X        X     1     0.000   1.000
""",
    )
    return par, frag


def test_acceptance_formula_insertion_and_deletion_ideal_gas(gcmc_cpu, temp_dir):
    work = Path(temp_dir) / "physical_contracts" / "acceptance"
    work.mkdir(parents=True, exist_ok=True)

    par, frag = _write_ideal_gas_itp(work)

    out_prefix = work / "ins" / "gcmc"
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    accept_log = work / "ins" / "acceptance.jsonl"

    inp = work / "ins" / "test.inp"
    _write_inp(
        inp,
        f"""
random_seed:123
par:{par}
fragitp:{frag}
fragname:X
fragconc:0.5
fragmuex:0.0

box_size:10.0 10.0 10.0
gcmc_region:box 0 0 0 10 10 10
cutoff:4.0
temperature:300.0
moves_per_step:1
mcsteps:1
nprint:1
mc_move_prob:1 0 0 0
""",
    )

    result = _run_gcmc_cpu(
        gcmc_cpu,
        workdir=work / "ins",
        inp=inp,
        out_prefix=out_prefix,
        extra_args=["--dump-accept", str(accept_log)],
        timeout=30,
    )
    assert result.returncode == 0, result.stdout + result.stderr

    rec_ins = _first_accept_record(accept_log, move="insertion", species="X")
    assert float(rec_ins["deltaU"]) == pytest.approx(0.0, abs=1e-12)
    expected_ins = min(
        1.0,
        float(rec_ins["z"]) * float(rec_ins["vBox"]) / (int(rec_ins["nBefore"]) + 1),
    )
    assert float(rec_ins["pAcc"]) == pytest.approx(expected_ins, rel=1e-12, abs=1e-12)

    del_prefix = work / "del" / "gcmc"
    del_prefix.parent.mkdir(parents=True, exist_ok=True)
    del_accept = work / "del" / "acceptance.jsonl"
    del_inp = work / "del" / "test.inp"
    _write_inp(
        del_inp,
        f"""
random_seed:456
par:{par}
fragitp:{frag}
fragname:X
fragconc:3.321
fragmuex:0.0

box_size:10.0 10.0 10.0
gcmc_region:box 0 0 0 10 10 10
cutoff:4.0
temperature:300.0
moves_per_step:1
mcsteps:50
nprint:50
mc_move_prob:1 1 0 0
""",
    )

    result_del = _run_gcmc_cpu(
        gcmc_cpu,
        workdir=work / "del",
        inp=del_inp,
        out_prefix=del_prefix,
        extra_args=["--dump-accept", str(del_accept)],
        timeout=30,
    )
    assert result_del.returncode == 0, result_del.stdout + result_del.stderr

    del_records = _load_accept_records(del_accept)
    rec_del = next(
        (r for r in del_records
         if r.get("move") == "deletion"
         and str(r.get("species", "")).strip().upper() == "X"
         and int(r.get("nBefore", 0)) > 0),
        None,
    )
    assert rec_del is not None, "Expected at least one deletion record with nBefore > 0"
    assert float(rec_del["deltaU"]) == pytest.approx(0.0, abs=1e-12)
    expected_del = min(
        1.0,
        int(rec_del["nBefore"]) / (float(rec_del["z"]) * float(rec_del["vBox"])),
    )
    assert float(rec_del["pAcc"]) == pytest.approx(expected_del, rel=1e-12, abs=1e-12)


def test_cutoff_excludes_far_lj_interactions_in_insertion_energy(gcmc_cpu, temp_dir):
    work = Path(temp_dir) / "physical_contracts" / "cutoff"
    work.mkdir(parents=True, exist_ok=True)

    pdb = work / "sys.pdb"
    _write_text(
        pdb,
        """
CRYST1   40.000   40.000   40.000  90.00  90.00  90.00 P 1           1
ATOM      1  C   MOL A   1       5.000   5.000   5.000  1.00  0.00           C
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
; nr  type  resnr  residue  atom  cgnr  charge  mass
1   C     1      MOL      C     1     0.000   12.011

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
1 2 yes 0.5 0.8333

[ atomtypes ]
; name  at.num  mass   charge  ptype  sigma   epsilon
C       0       12.011 0.000   A      0.340   0.300
X       0       1.000  0.000   A      0.250   0.200
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

    out_prefix = work / "out" / "gcmc"
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    accept_log = work / "out" / "acceptance.jsonl"

    inp = work / "test.inp"
    _write_inp(
        inp,
        f"""
random_seed:1357
par:{par}
fragitp:{frag}
fragname:FRG
fragconc:3.5
fragmuex:0.0

pdb:{pdb}
top:{top}
box_size:40.0 40.0 40.0
    gcmc_region:box 20 20 20 30 30 30
    cutoff:6.0
    temperature:1000000.0
moves_per_step:1
mcsteps:1
nprint:1
mc_move_prob:1 0 0 0
""",
    )

    result = _run_gcmc_cpu(
        gcmc_cpu,
        workdir=work,
        inp=inp,
        out_prefix=out_prefix,
        extra_args=["--dump-accept", str(accept_log)],
        timeout=30,
    )
    assert result.returncode == 0, result.stdout + result.stderr

    rec = _first_accept_record(accept_log, move="insertion", species="FRG")
    assert bool(rec.get("accepted")) is True

    final_pdb = Path(f"{out_prefix}_final.pdb")
    assert final_pdb.exists()

    box_ang = _read_cryst1_box_angstrom(final_pdb)
    c_pos = _first_atom_xyz_angstrom(final_pdb, resname="MOL")
    x_pos = _first_atom_xyz_angstrom(final_pdb, resname="FRG")
    dx, dy, dz = _min_image_delta_nm(c_pos, x_pos, box_ang=box_ang)
    r_nm = math.sqrt(dx * dx + dy * dy + dz * dz)

    cutoff_nm = 6.0 / 10.0
    assert r_nm > cutoff_nm + 0.4

    sigma_mix = 0.5 * (0.340 + 0.250)
    eps_mix = math.sqrt(0.300 * 0.200)
    expected = 0.0 if r_nm > cutoff_nm else _lj_energy_kj_mol(
        r_nm=r_nm, sigma_nm=sigma_mix, eps_kj_mol=eps_mix
    )
    assert float(rec["deltaU"]) == pytest.approx(expected, abs=1e-6)


def test_insertion_deltaU_matches_analytic_coulomb_energy(gcmc_cpu, temp_dir):
    work = Path(temp_dir) / "physical_contracts" / "coulomb_insertion"
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
; nr  type  resnr  residue  atom  cgnr  charge  mass
1   C     1      MOL      C     1     0.500   12.011

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
1 2 yes 0.5 0.8333

[ atomtypes ]
; name  at.num  mass   charge  ptype  sigma   epsilon
C       0       12.011 0.000   A      0.320   0.000
X       0       1.000  0.000   A      0.280   0.000
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
1   X     1      FRG      X     1     -0.250  1.000
""",
    )

    out_prefix = work / "out" / "gcmc"
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    accept_log = work / "out" / "acceptance.jsonl"

    inp = work / "test.inp"
    _write_inp(
        inp,
        f"""
random_seed:2468
par:{par}
fragitp:{frag}
fragname:FRG
fragconc:0.0
fragmuex:50.0

pdb:{pdb}
top:{top}
box_size:30.0 30.0 30.0
gcmc_region:box 14.0 14.0 14.0 16.0 16.0 16.0
cutoff:6.0
temperature:300.0
moves_per_step:1
mcsteps:1
nprint:1
mc_move_prob:1 0 0 0
""",
    )

    result = _run_gcmc_cpu(
        gcmc_cpu,
        workdir=work,
        inp=inp,
        out_prefix=out_prefix,
        extra_args=["--dump-accept", str(accept_log)],
        timeout=30,
    )
    assert result.returncode == 0, result.stdout + result.stderr

    rec = _first_accept_record(accept_log, move="insertion", species="FRG")
    assert bool(rec.get("accepted")) is True

    final_pdb = Path(f"{out_prefix}_final.pdb")
    assert final_pdb.exists()

    box_ang = _read_cryst1_box_angstrom(final_pdb)
    mol_pos = _first_atom_xyz_angstrom(final_pdb, resname="MOL")
    frg_pos = _first_atom_xyz_angstrom(final_pdb, resname="FRG")
    dx, dy, dz = _min_image_delta_nm(mol_pos, frg_pos, box_ang=box_ang)
    r_nm = math.sqrt(dx * dx + dy * dy + dz * dz)

    cutoff_nm = 6.0 / 10.0
    assert r_nm < cutoff_nm

    expected = _coulomb_energy_kj_mol(r_nm=r_nm, q1=0.500, q2=-0.250)
    assert float(rec["deltaU"]) == pytest.approx(expected, rel=1e-3, abs=1e-2)


def test_poisson_number_distribution_ideal_gas():
    """
    Delegate to the engine-level Poisson distribution test.

    This keeps the physical-contract coverage while avoiding CLI move-selection
    side effects that can bias the stationary distribution.
    """
    _ensure_pygcmc_on_path()
    try:
        import pygcmc  # noqa: F401
    except ModuleNotFoundError as exc:
        raise AssertionError(
            "pygcmc bindings not found; build with `cmake --build build --target pygcmc`"
        ) from exc

    from movementGCMC.poisson_distribution import test_ideal_gas_poisson_distribution

    test_ideal_gas_poisson_distribution()

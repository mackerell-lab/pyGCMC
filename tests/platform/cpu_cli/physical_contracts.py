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
import statistics
import sys
from pathlib import Path

import pytest

from .inp_units_compat_helpers import _run_gcmc_cpu, _write_inp

COULOMB = 138.935458  # kJ·nm/mol/e^2
MOLAR_TO_NM3 = 0.6022140857  # mol/L -> particles/nm^3


def _write_text(path: Path, content: str) -> None:
    path.write_text(content.strip() + "\n")


def _ensure_pygcmc_on_path() -> None:
    repo_root = Path(__file__).resolve().parents[3]
    bindings = repo_root / "build" / "modules" / "bindings"
    if bindings.exists():
        bindings_str = str(bindings)
        if bindings_str not in sys.path:
            sys.path.insert(0, bindings_str)


def _parse_itp_atomtypes(path: Path) -> dict[str, tuple[float, float]]:
    atomtypes: dict[str, tuple[float, float]] = {}
    in_section = False
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(("[", ";", "#")):
            if line.startswith("["):
                in_section = line.lower().startswith("[ atomtypes")
            continue
        if not in_section:
            continue
        if ";" in line:
            line = line.split(";", 1)[0].strip()
        if not line:
            continue
        parts = line.split()
        if len(parts) < 7:
            continue
        name = parts[0]
        sigma = float(parts[-2])
        eps = float(parts[-1])
        atomtypes[name] = (sigma, eps)
    if not atomtypes:
        raise AssertionError(f"No atomtypes parsed from {path}")
    return atomtypes


def _build_state_from_pdb_top(
    par_path: Path,
    pdb_path: Path,
    top_path: Path,
    *,
    cutoff_nm: float,
    box_nm: tuple[float, float, float],
):
    _ensure_pygcmc_on_path()
    try:
        import pygcmc
    except ModuleNotFoundError as exc:
        raise AssertionError(
            "pygcmc bindings not found; build with `cmake --build build --target pygcmc`"
        ) from exc

    structure = pygcmc.PDBParser.parse_file(str(pdb_path))
    topology = pygcmc.TOPParser.parse_file(str(top_path))
    molecular = pygcmc.MolecularSystem().combine(structure, topology)

    mc_system = pygcmc.MonteCarloSystem()
    info = pygcmc.MCInfo()
    info.max_residues = 50
    info.max_atoms = 200
    mc_system.initialize(info)
    mc_system.initialize_from_molecular(molecular)

    state = mc_system.get_state_mutable()
    state.info.cutoff = float(cutoff_nm)
    state.info.box = [float(x) for x in box_nm]

    atomtypes = _parse_itp_atomtypes(par_path)
    for name in sorted(atomtypes):
        state.atomTypes.get_or_add_type(name)

    num_types = len(state.atomTypes.atomTypes)
    state.forcefield.numTotalTypes = num_types
    state.forcefield.numMovementTypes = num_types
    state.forcefield.mixingRule = pygcmc.MCForceField.MixingRule.LorentzBerthelot
    state.forcefield.ljSigmaType = [0.0] * num_types
    state.forcefield.ljEpsType = [0.0] * num_types

    for name in sorted(atomtypes):
        sigma, eps = atomtypes[name]
        idx = state.atomTypes.get_or_add_type(name)
        if 0 <= idx < num_types:
            state.forcefield.ljSigmaType[idx] = float(sigma)
            state.forcefield.ljEpsType[idx] = float(eps)

    state.forcefield.rebuildLJMatrix()
    return mc_system, state


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


def _read_statistics_n_total(stats_path: Path) -> list[int]:
    n_total: list[int] = []
    for line in stats_path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 3:
            continue
        n_total.append(int(parts[2]))
    return n_total


def _lj_energy_kj_mol(*, r_nm: float, sigma_nm: float, eps_kj_mol: float) -> float:
    r2 = r_nm * r_nm
    sigma_r2 = (sigma_nm * sigma_nm) / r2
    sigma_r6 = sigma_r2 * sigma_r2 * sigma_r2
    sigma_r12 = sigma_r6 * sigma_r6
    return 4.0 * eps_kj_mol * (sigma_r12 - sigma_r6)


def _coulomb_energy_kj_mol(*, r_nm: float, q1: float, q2: float) -> float:
    return COULOMB * q1 * q2 / r_nm


def _expected_pacc_from_record(rec: dict) -> float | None:
    move = str(rec.get("move", "")).strip().lower()
    beta = float(rec["beta"])
    deltaU = float(rec["deltaU"])

    if move in ("insertion", "deletion"):
        n_before = int(rec["nBefore"])
        if move == "deletion" and n_before <= 0:
            return 0.0
        z = max(float(rec["z"]), 1e-30)
        v_eff = max(float(rec.get("vEff", rec["vBox"])), 1e-30)
        q_forward = max(float(rec["qForward"]), 1e-30)
        q_reverse = max(float(rec["qReverse"]), 1e-30)
        proposal_ratio = float(rec.get("proposalRatio", 1.0))

        if move == "insertion":
            ratio = (z * v_eff / (n_before + 1)) * math.exp(-beta * deltaU)
        else:
            ratio = (n_before / (z * v_eff)) * math.exp(-beta * deltaU)

        ratio *= (q_forward / q_reverse) * proposal_ratio
        return min(1.0, ratio)

    if move in ("translation", "rotation"):
        bias = float(rec.get("bias", 1.0))
        return min(1.0, math.exp(-beta * deltaU) * bias)

    return None


def _expected_pacc_with_cavity(rec: dict, move: str) -> float:
    expected = _expected_pacc_from_record(rec)
    if expected is None:
        raise AssertionError(f"Unsupported move for acceptance closure: {move}")
    return expected


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


def _write_charged_single_atom_system(
    work: Path,
    *,
    mol_charge: float,
    frag_charge: float,
    mol_sigma_nm: float,
    mol_eps_kj_mol: float,
    frag_sigma_nm: float,
    frag_eps_kj_mol: float,
) -> dict[str, Path]:
    pdb = work / "sys.pdb"
    _write_text(
        pdb,
        """
CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1           1
ATOM      1  C   MOL A   1      15.000  15.000  15.000  1.00  0.00           C
END
""",
    )

    par = work / "par.itp"
    _write_text(
        par,
        f"""
[ defaults ]
1 2 yes 0.5 0.8333

[ atomtypes ]
; name  at.num  mass   charge  ptype  sigma   epsilon
C       0       12.011 0.000   A      {mol_sigma_nm:.6f}   {mol_eps_kj_mol:.6f}
X       0       1.000  0.000   A      {frag_sigma_nm:.6f}   {frag_eps_kj_mol:.6f}
""",
    )

    frag = work / "frag.itp"
    _write_text(
        frag,
        f"""
[ moleculetype ]
FRG  1

[ atoms ]
; nr  type  resnr  residue  atom  cgnr  charge  mass
1   X     1      FRG      X     1     {frag_charge:.6f}   1.000
""",
    )

    top_ins = work / "sys_ins.top"
    _write_text(
        top_ins,
        f"""
[ defaults ]
1 2 yes 0.5 0.8333

[ moleculetype ]
MOL  2

[ atoms ]
; nr  type  resnr  residue  atom  cgnr  charge  mass
1   C     1      MOL      C     1     {mol_charge:.6f}   12.011

[ moleculetype ]
FRG  1

[ atoms ]
; nr  type  resnr  residue  atom  cgnr  charge  mass
1   X     1      FRG      X     1     {frag_charge:.6f}   1.000

[ system ]
Minimal

[ molecules ]
MOL  1
""",
    )

    top_del = work / "sys_del.top"
    _write_text(
        top_del,
        f"""
[ defaults ]
1 2 yes 0.5 0.8333

[ moleculetype ]
MOL  2

[ atoms ]
; nr  type  resnr  residue  atom  cgnr  charge  mass
1   C     1      MOL      C     1     {mol_charge:.6f}   12.011

[ moleculetype ]
FRG  1

[ atoms ]
; nr  type  resnr  residue  atom  cgnr  charge  mass
1   X     1      FRG      X     1     {frag_charge:.6f}   1.000

[ system ]
Minimal

[ molecules ]
MOL  1
FRG  1
""",
    )

    return {
        "pdb": pdb,
        "par": par,
        "frag": frag,
        "top_ins": top_ins,
        "top_del": top_del,
    }


def _write_cavity_bias_system(work: Path) -> dict[str, Path]:
    pdb = work / "cavity.pdb"
    _write_text(
        pdb,
        """
CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1           1
ATOM      1  C   MOL A   1      10.000  10.000  10.000  1.00  0.00           C
ATOM      2  X   FRG A   2      20.000  20.000  20.000  1.00  0.00           C
END
""",
    )

    par = work / "cavity_par.itp"
    _write_text(
        par,
        """
[ defaults ]
1 2 yes 0.5 0.8333

[ atomtypes ]
; name  at.num  mass   charge  ptype  sigma   epsilon
C       0       12.011 0.000   A      0.300   0.200
X       0       1.000  0.000   A      0.250   0.100
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

    top = work / "cavity.top"
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

[ moleculetype ]
FRG  1

[ atoms ]
; nr  type  resnr  residue  atom  cgnr  charge  mass
1   X     1      FRG      X     1     0.000   1.000

[ system ]
CavityBias

[ molecules ]
MOL  1
FRG  1
""",
    )

    return {"pdb": pdb, "par": par, "frag": frag, "top": top}


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


def test_acceptance_log_pacc_matches_formula_all_moves(gcmc_cpu, temp_dir):
    work = Path(temp_dir) / "physical_contracts" / "acceptance_all_moves"
    work.mkdir(parents=True, exist_ok=True)

    files = _write_charged_single_atom_system(
        work,
        mol_charge=0.0,
        frag_charge=0.0,
        mol_sigma_nm=0.34,
        mol_eps_kj_mol=0.20,
        frag_sigma_nm=0.28,
        frag_eps_kj_mol=0.15,
    )

    out_prefix = work / "run" / "gcmc"
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    accept_log = work / "run" / "acceptance.jsonl"

    inp = work / "run" / "test.inp"
    _write_inp(
        inp,
        f"""
random_seed:2468
par:{files["par"]}
fragitp:{files["frag"]}
fragname:FRG
fragconc:0.5
fragmuex:0.0

pdb:{files["pdb"]}
top:{files["top_ins"]}
box_size:30.0 30.0 30.0
cutoff:8.0
temperature:300.0
moves_per_step:2
mcsteps:150
nprint:150
mc_move_prob:1 1 1 1
""",
    )

    result = _run_gcmc_cpu(
        gcmc_cpu,
        workdir=work / "run",
        inp=inp,
        out_prefix=out_prefix,
        extra_args=["--dump-accept", str(accept_log)],
        timeout=60,
    )
    assert result.returncode == 0, result.stdout + result.stderr

    records = _load_accept_records(accept_log)
    assert records, "Acceptance log unexpectedly empty"

    seen = {"insertion": 0, "deletion": 0, "translation": 0, "rotation": 0}
    for rec in records:
        move = str(rec.get("move", "")).strip().lower()
        if move not in seen:
            continue
        if float(rec.get("pAcc", -1.0)) < 0.0:
            continue
        expected = _expected_pacc_from_record(rec)
        assert expected is not None
        assert float(rec["pAcc"]) == pytest.approx(expected, rel=1e-6, abs=1e-10)
        seen[move] += 1

    for move, count in seen.items():
        assert count >= 1, f"Expected at least one {move} record"


def test_cavity_bias_pacc_matches_acceptance_formula(gcmc_cpu, temp_dir):
    work = Path(temp_dir) / "physical_contracts" / "cavity_bias"
    work.mkdir(parents=True, exist_ok=True)

    files = _write_cavity_bias_system(work)

    out_prefix = work / "run" / "gcmc"
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    accept_log = work / "run" / "acceptance.jsonl"

    inp = work / "run" / "test.inp"
    _write_inp(
        inp,
        f"""
random_seed:777
par:{files["par"]}
fragitp:{files["frag"]}
fragname:FRG
fragconc:1.0
fragmuex:0.0

pdb:{files["pdb"]}
top:{files["top"]}
box_size:30.0 30.0 30.0
cutoff:8.0
grid_dx:2.0
probe_radius:1.4
use_cavity_bias:yes
use_vdw_radius_for_grid:yes
exclude_hydrogens_from_grid:no
temperature:300.0
moves_per_step:1
mcsteps:200
nprint:200
mc_move_prob:0.5 0.5 0 0
""",
    )

    result = _run_gcmc_cpu(
        gcmc_cpu,
        workdir=work / "run",
        inp=inp,
        out_prefix=out_prefix,
        extra_args=["--dump-accept", str(accept_log)],
        timeout=60,
    )
    assert result.returncode == 0, result.stdout + result.stderr

    records = _load_accept_records(accept_log)
    ins_rec = next(
        (r for r in records
         if r.get("move") == "insertion"
         and str(r.get("species", "")).strip().upper() == "FRG"),
        None,
    )
    del_rec = next(
        (r for r in records
         if r.get("move") == "deletion"
         and str(r.get("species", "")).strip().upper() == "FRG"
         and int(r.get("nBefore", 0)) > 0),
        None,
    )
    assert ins_rec is not None, "Expected at least one FRG insertion record"
    assert del_rec is not None, "Expected at least one FRG deletion record with nBefore > 0"

    expected_ins = _expected_pacc_from_record(ins_rec)
    assert float(ins_rec["pAcc"]) == pytest.approx(expected_ins, rel=1e-6, abs=1e-6)

    expected_del = _expected_pacc_from_record(del_rec)
    assert float(del_rec["pAcc"]) == pytest.approx(expected_del, rel=1e-6, abs=1e-6)


def test_cbmc_cavity_bias_acceptance_closure(gcmc_cpu, temp_dir):
    work = Path(temp_dir) / "physical_contracts" / "cbmc_cavity"
    work.mkdir(parents=True, exist_ok=True)

    files = _write_cavity_bias_system(work)

    out_prefix = work / "run" / "gcmc"
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    accept_log = work / "run" / "acceptance.jsonl"
    params_json = work / "run" / "params.json"

    inp = work / "run" / "test.inp"
    _write_inp(
        inp,
        f"""
random_seed:4242
par:{files["par"]}
fragitp:{files["frag"]}
fragname:FRG
fragconc:0.5
fragmuex:0.0

pdb:{files["pdb"]}
top:{files["top"]}
box_size:30.0 30.0 30.0
cutoff:8.0
grid_dx:2.0
probe_radius:2.8
use_cavity_bias:yes
use_conf_bias:yes
num_conf_bias_trial:5
use_vdw_radius_for_grid:yes
exclude_hydrogens_from_grid:no
temperature:300.0
moves_per_step:1
mcsteps:200
nprint:200
mc_move_prob:1 1 0 0
""",
    )

    result = _run_gcmc_cpu(
        gcmc_cpu,
        workdir=work / "run",
        inp=inp,
        out_prefix=out_prefix,
        extra_args=["--dump-accept", str(accept_log), "--dump-params", str(params_json)],
        timeout=90,
    )
    assert result.returncode == 0, result.stdout + result.stderr

    params = json.loads(params_json.read_text())
    assert params["bias"]["use_cavity_bias"] is True
    assert params["bias"]["use_conf_bias"] is True
    assert int(params["bias"]["num_conf_bias_trials"]) == 5

    records = _load_accept_records(accept_log)
    insertions = [
        r for r in records
        if r.get("move") == "insertion"
        and str(r.get("species", "")).strip().upper() == "FRG"
    ]
    deletions = [
        r for r in records
        if r.get("move") == "deletion"
        and str(r.get("species", "")).strip().upper() == "FRG"
        and int(r.get("nBefore", 0)) > 0
    ]
    assert insertions, "Expected insertion records for FRG"
    assert deletions, "Expected deletion records for FRG with nBefore > 0"

    max_trials = max(int(r.get("cbmcTrials", 0)) for r in insertions)
    assert max_trials >= 5
    min_cavity = min(float(r.get("cavityFraction", 1.0)) for r in insertions)
    assert 0.0 < min_cavity < 0.999

    for rec in insertions[:10] + deletions[:10]:
        expected = _expected_pacc_from_record(rec)
        assert expected is not None
        assert float(rec["pAcc"]) == pytest.approx(expected, rel=1e-6, abs=1e-10)


def test_mu_targets_stable_concentration_with_insertion_and_deletion(gcmc_cpu, temp_dir):
    work = Path(temp_dir) / "physical_contracts" / "mu_stability"
    work.mkdir(parents=True, exist_ok=True)

    par, frag = _write_ideal_gas_itp(work)

    out_prefix = work / "out" / "gcmc"
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    accept_log = work / "out" / "acceptance.jsonl"

    target_mean = 3.0
    fragconc_m = 10.0
    beta = 1.0 / (0.008314462618 * 300.0)
    target_z = target_mean  # V = 1 nm^3
    mu_kj = math.log(target_z / (fragconc_m * MOLAR_TO_NM3)) / beta
    fragmuex_kcal = mu_kj / 4.184

    inp = work / "test.inp"
    _write_inp(
        inp,
        f"""
random_seed:314159
par:{par}
fragitp:{frag}
fragname:X
fragconc:{fragconc_m}
fragmuex:{fragmuex_kcal}

box_size:10.0 10.0 10.0
gcmc_region:box 0 0 0 10 10 10
cutoff:4.0
temperature:300.0
moves_per_step:1
mcsteps:2000
nprint:2000
mc_move_prob:1 1 0 0
""",
    )

    result = _run_gcmc_cpu(
        gcmc_cpu,
        workdir=work,
        inp=inp,
        out_prefix=out_prefix,
        extra_args=["--dump-accept", str(accept_log)],
        timeout=60,
    )
    assert result.returncode == 0, result.stdout + result.stderr

    records = _load_accept_records(accept_log)
    assert records, "No acceptance records were produced"

    species_records = [
        r for r in records
        if str(r.get("species", "")).strip().upper() == "X"
    ]
    assert species_records, "No acceptance records for species X"

    expected_z = fragconc_m * MOLAR_TO_NM3 * math.exp(beta * mu_kj)
    log_z = float(species_records[0]["z"])
    assert log_z == pytest.approx(expected_z, rel=2e-3, abs=1e-4)

    burnin = 200
    sample_records = species_records[burnin:] if len(species_records) > burnin else species_records

    n_values = [int(r["nBefore"]) for r in sample_records if "nBefore" in r]
    assert n_values, "No nBefore values found in acceptance records"

    mean_n = statistics.mean(n_values)
    expected_mean = expected_z * float(species_records[0]["vBox"])
    assert mean_n == pytest.approx(expected_mean, rel=0.25, abs=0.5)

    ins_attempts = [r for r in sample_records if str(r.get("move", "")).strip().lower() == "insertion"]
    del_attempts = [r for r in sample_records if str(r.get("move", "")).strip().lower() == "deletion"]
    assert len(ins_attempts) > 0
    assert len(del_attempts) > 0

    ins_accepts = [r for r in ins_attempts if bool(r.get("accepted"))]
    del_accepts = [r for r in del_attempts if bool(r.get("accepted"))]
    assert len(ins_accepts) > 0
    assert len(del_accepts) > 0


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


def test_insertion_deltaU_matches_analytic_lj_plus_coulomb_energy(gcmc_cpu, temp_dir):
    work = Path(temp_dir) / "physical_contracts" / "lj_coulomb_insertion"
    work.mkdir(parents=True, exist_ok=True)

    files = _write_charged_single_atom_system(
        work,
        mol_charge=0.500,
        frag_charge=-0.250,
        mol_sigma_nm=0.300,
        mol_eps_kj_mol=0.200,
        frag_sigma_nm=0.280,
        frag_eps_kj_mol=0.100,
    )

    out_prefix = work / "out" / "gcmc"
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    accept_log = work / "out" / "acceptance.jsonl"

    inp = work / "test.inp"
    _write_inp(
        inp,
        f"""
random_seed:97531
par:{files["par"]}
fragitp:{files["frag"]}
fragname:FRG
fragconc:0.0
fragmuex:30.0

pdb:{files["pdb"]}
top:{files["top_ins"]}
box_size:30.0 30.0 30.0
gcmc_region:box 18.0 15.0 15.0 19.0 16.0 16.0
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

    sigma_mix = 0.5 * (0.300 + 0.280)
    eps_mix = math.sqrt(0.200 * 0.100)
    expected = _lj_energy_kj_mol(
        r_nm=r_nm, sigma_nm=sigma_mix, eps_kj_mol=eps_mix
    ) + _coulomb_energy_kj_mol(r_nm=r_nm, q1=0.500, q2=-0.250)
    assert float(rec["deltaU"]) == pytest.approx(expected, rel=1e-3, abs=1e-2)

    _ensure_pygcmc_on_path()
    try:
        import pygcmc
    except ModuleNotFoundError as exc:
        raise AssertionError(
            "pygcmc bindings not found; build with `cmake --build build --target pygcmc`"
        ) from exc

    cutoff_nm = 6.0 / 10.0
    box_nm = (3.0, 3.0, 3.0)
    mc_before, state_before = _build_state_from_pdb_top(
        files["par"],
        files["pdb"],
        files["top_ins"],
        cutoff_nm=cutoff_nm,
        box_nm=box_nm,
    )
    pygcmc.computeSystemEnergyPBCCutoff(state_before)
    elec_before, vdw_before = pygcmc.getTotalEnergyComponents(state_before)

    mc_after, state_after = _build_state_from_pdb_top(
        files["par"],
        final_pdb,
        files["top_del"],
        cutoff_nm=cutoff_nm,
        box_nm=box_nm,
    )
    pygcmc.computeSystemEnergyPBCCutoff(state_after)
    elec_after, vdw_after = pygcmc.getTotalEnergyComponents(state_after)

    delta_components = (elec_after - elec_before) + (vdw_after - vdw_before)
    assert float(rec["deltaU"]) == pytest.approx(delta_components, rel=5e-3, abs=0.15)


def test_deletion_deltaU_is_negative_of_insertion_for_charged_system(gcmc_cpu, temp_dir):
    work = Path(temp_dir) / "physical_contracts" / "charged_deltaU_symmetry"
    work.mkdir(parents=True, exist_ok=True)

    files = _write_charged_single_atom_system(
        work,
        mol_charge=0.400,
        frag_charge=-0.200,
        mol_sigma_nm=0.320,
        mol_eps_kj_mol=0.000,
        frag_sigma_nm=0.280,
        frag_eps_kj_mol=0.000,
    )

    ins_prefix = work / "ins" / "gcmc"
    ins_prefix.parent.mkdir(parents=True, exist_ok=True)
    ins_accept = work / "ins" / "acceptance.jsonl"
    ins_inp = work / "ins" / "test.inp"
    _write_inp(
        ins_inp,
        f"""
random_seed:8642
par:{files["par"]}
fragitp:{files["frag"]}
fragname:FRG
fragconc:0.0
fragmuex:30.0

pdb:{files["pdb"]}
top:{files["top_ins"]}
box_size:30.0 30.0 30.0
gcmc_region:box 18.0 15.0 15.0 19.0 16.0 16.0
cutoff:6.0
temperature:300.0
moves_per_step:1
mcsteps:1
nprint:1
mc_move_prob:1 0 0 0
""",
    )

    ins_result = _run_gcmc_cpu(
        gcmc_cpu,
        workdir=work / "ins",
        inp=ins_inp,
        out_prefix=ins_prefix,
        extra_args=["--dump-accept", str(ins_accept)],
        timeout=30,
    )
    assert ins_result.returncode == 0, ins_result.stdout + ins_result.stderr

    ins_rec = _first_accept_record(ins_accept, move="insertion", species="FRG")
    assert bool(ins_rec.get("accepted")) is True

    deltaU_ins = float(ins_rec["deltaU"])
    assert abs(deltaU_ins) > 1e-6

    final_pdb = Path(f"{ins_prefix}_final.pdb")
    assert final_pdb.exists()

    del_prefix = work / "del" / "gcmc"
    del_prefix.parent.mkdir(parents=True, exist_ok=True)
    del_accept = work / "del" / "acceptance.jsonl"
    del_inp = work / "del" / "test.inp"
    _write_inp(
        del_inp,
        f"""
random_seed:97531
par:{files["par"]}
fragitp:{files["frag"]}
fragname:FRG
fragconc:0.0
fragmuex:30.0

pdb:{final_pdb}
top:{files["top_del"]}
box_size:30.0 30.0 30.0
gcmc_region:box 18.0 15.0 15.0 19.0 16.0 16.0
cutoff:6.0
temperature:300.0
moves_per_step:1
mcsteps:1
nprint:1
mc_move_prob:0 1 0 0
""",
    )

    del_result = _run_gcmc_cpu(
        gcmc_cpu,
        workdir=work / "del",
        inp=del_inp,
        out_prefix=del_prefix,
        extra_args=["--dump-accept", str(del_accept)],
        timeout=30,
    )
    assert del_result.returncode == 0, del_result.stdout + del_result.stderr

    del_rec = _first_accept_record(del_accept, move="deletion", species="FRG")
    assert int(del_rec.get("nBefore", 0)) == 1

    deltaU_del = float(del_rec["deltaU"])
    assert deltaU_del == pytest.approx(-deltaU_ins, rel=1e-4, abs=1e-2)


def test_poisson_number_distribution_ideal_gas(gcmc_cpu, temp_dir):
    """
    CLI Poisson distribution check using --dump-accept records.

    Validates <N> = z * V and Var(N) ≈ <N> in the ideal-gas limit.
    """
    work = Path(temp_dir) / "physical_contracts" / "poisson_cli"
    work.mkdir(parents=True, exist_ok=True)

    par, frag = _write_ideal_gas_itp(work)

    target_mean = 5.0
    beta = 1.0 / (0.008314462618 * 300.0)
    target_z = target_mean  # V = 1 nm^3 for 10 Å box
    mu_kj = math.log(target_z) / beta
    fragmuex_kcal = mu_kj / 4.184

    out_prefix = work / "out" / "gcmc"
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    accept_log = work / "out" / "acceptance.jsonl"

    inp = work / "run.inp"
    _write_inp(
        inp,
        f"""
random_seed:424242
par:{par}
fragitp:{frag}
fragname:X
fragconc:0.0
fragmuex:{fragmuex_kcal}

box_size:10.0 10.0 10.0
gcmc_region:box 0 0 0 10 10 10
cutoff:4.0
temperature:300.0
moves_per_step:1
mcsteps:4000
nprint:1
mc_move_prob:1 1 0 0
""",
    )

    result = _run_gcmc_cpu(
        gcmc_cpu,
        workdir=work,
        inp=inp,
        out_prefix=out_prefix,
        extra_args=["--dump-accept", str(accept_log)],
        timeout=60,
    )
    assert result.returncode == 0, result.stdout + result.stderr

    records = _load_accept_records(accept_log)
    species_records = [
        r for r in records
        if str(r.get("species", "")).strip().upper() == "X"
    ]
    assert species_records, "No acceptance records for species X"

    expected_mean = float(species_records[0]["z"]) * float(species_records[0]["vBox"])

    stats_path = Path(f"{out_prefix}_statistics.dat")
    assert stats_path.exists(), "statistics.dat missing for Poisson check"
    n_values = _read_statistics_n_total(stats_path)
    assert len(n_values) >= 1000, "Insufficient statistics.dat samples for Poisson check"

    burnin = 200
    n_values = n_values[burnin:] if len(n_values) > burnin else n_values
    assert len(n_values) >= 500, "Insufficient post-burnin samples for Poisson statistics"

    sample_mean = statistics.mean(n_values)
    sample_var = statistics.pvariance(n_values)

    assert sample_mean == pytest.approx(expected_mean, rel=0.25, abs=0.5)
    assert abs(sample_var - expected_mean) / max(expected_mean, 1e-6) < 0.35

"""
INP compatibility tests for gcmc_gpu-style units and legacy keys.

These are end-to-end (CLI) tests that exercise InpParserGCMC via gcmc_cpu.
"""

from __future__ import annotations

import json
import math
import subprocess
from pathlib import Path

import pytest


def _write_inp(path: Path, content: str) -> None:
    path.write_text(content.strip() + "\n")


def _read_cryst1_box_angstrom(pdb_path: Path) -> tuple[float, float, float]:
    for line in pdb_path.read_text().splitlines():
        if line.startswith("CRYST1"):
            parts = line.split()
            assert len(parts) >= 4, f"Unexpected CRYST1 format: {line}"
            return float(parts[1]), float(parts[2]), float(parts[3])
    raise AssertionError(f"CRYST1 not found in {pdb_path}")

def _count_residues_by_resname(pdb_path: Path, resname: str) -> int:
    want = resname.strip().upper()
    resids: set[int] = set()
    for line in pdb_path.read_text().splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        if line[17:20].strip().upper() != want:
            continue
        resids.add(int(line[22:26]))
    return len(resids)

def _run_gcmc_cpu(
    gcmc_cpu: str,
    *,
    workdir: Path,
    inp: Path,
    out_prefix: Path,
    extra_args: list[str] | None = None,
    timeout: int = 20,
) -> subprocess.CompletedProcess[str]:
    args = [gcmc_cpu, "--inp", str(inp), "--prefix", str(out_prefix)]
    if extra_args:
        args.extend(extra_args)
    return subprocess.run(
        args,
        cwd=str(workdir),
        capture_output=True,
        text=True,
        timeout=timeout,
    )


def test_inp_units_gcmc_gpu_box_roundtrip_cryst1(gcmc_cpu, test_data_dir, temp_dir):
    """If INP is in Å, internal nm should roundtrip back to Å in output CRYST1."""
    work = Path(temp_dir) / "units_roundtrip"
    work.mkdir(parents=True, exist_ok=True)

    itp = test_data_dir / "charmm36.ff" / "mol" / "sol.itp"
    assert itp.exists()

    out_prefix = work / "out" / "gcmc"
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    inp = work / "test.inp"
    _write_inp(
        inp,
        f"""
inp_units:gcmc_gpu
fragitp:{itp}
box_size:30.0 31.0 32.0
cutoff:12.0
temperature:300.0
moves_per_step:1
mcsteps:5
nprint:1
fragname:SOL
fragmuex:0.0
attempt_prob_ins:1.0
attempt_prob_del:0.0
attempt_prob_trn:0.0
attempt_prob_rot:0.0
""",
    )

    result = subprocess.run(
        [gcmc_cpu, "--inp", str(inp), "--prefix", str(out_prefix), "--verbose"],
        cwd=str(work),
        capture_output=True,
        text=True,
        timeout=20,
    )
    assert result.returncode == 0, result.stdout + result.stderr

    out_pdb = Path(f"{out_prefix}_final.pdb")
    assert out_pdb.exists()

    lx, ly, lz = _read_cryst1_box_angstrom(out_pdb)
    assert lx == pytest.approx(30.0, abs=1e-3)
    assert ly == pytest.approx(31.0, abs=1e-3)
    assert lz == pytest.approx(32.0, abs=1e-3)


def test_inp_units_gcmc_gpu_fragmuex_kcal_to_kj_log(gcmc_cpu, test_data_dir, temp_dir):
    """fragmuex is kcal/mol in legacy mode and should be converted to kJ/mol internally."""
    work = Path(temp_dir) / "muex_units"
    work.mkdir(parents=True, exist_ok=True)

    itp = test_data_dir / "charmm36.ff" / "mol" / "sol.itp"
    assert itp.exists()

    out_prefix = work / "out" / "gcmc"
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    accept_log = work / "out" / "acceptance.jsonl"

    inp = work / "test.inp"
    _write_inp(
        inp,
        f"""
inp_units:gcmc_gpu
fragitp:{itp}
box_size:30.0 30.0 30.0
cutoff:12.0
temperature:300.0
moves_per_step:1
mcsteps:1
nprint:1
fragname:SOL
fragmuex:-1.0
attempt_prob_ins:1.0
attempt_prob_del:0.0
attempt_prob_trn:0.0
attempt_prob_rot:0.0
""",
    )

    result = _run_gcmc_cpu(
        gcmc_cpu,
        workdir=work,
        inp=inp,
        out_prefix=out_prefix,
        extra_args=["--dump-accept", str(accept_log)],
    )
    assert result.returncode == 0, result.stdout + result.stderr

    assert accept_log.exists()
    records = [
        json.loads(line)
        for line in accept_log.read_text().splitlines()
        if line.strip()
    ]
    assert records, "Acceptance log unexpectedly empty"

    mu_values = [
        float(r["mu"])
        for r in records
        if str(r.get("species", "")).strip().upper() == "SOL"
    ]
    assert mu_values, f"No SOL moves found in acceptance log: {records[:3]}"
    assert mu_values[0] == pytest.approx(-4.184, rel=1e-3, abs=1e-3)

def test_fragmuex_scales_activity_and_acceptance_in_ideal_gas_limit(
    gcmc_cpu, test_data_dir, temp_dir
):
    """
    Behavior-driven check that fragmuex influences acceptance probability via activity.

    Use an empty 1 nm^3 box with a single-atom fragment so deltaU ~= 0 and:
        pAcc = min(1, z * V / (nBefore+1))
    For small z (z<1), this reduces to pAcc == z, and changing μ should scale z by exp(beta*μ).
    """
    itp = test_data_dir / "charmm36.ff" / "mol" / "na.itp"
    assert itp.exists()

    conc_m = 0.1
    temperature_k = 300.0

    def run_once(muex_kcal_mol: float) -> dict:
        tag = "mu0" if muex_kcal_mol == 0.0 else "mu_nonzero"
        work = Path(temp_dir) / "muex_pacc_effect" / tag
        work.mkdir(parents=True, exist_ok=True)

        out_prefix = work / "out" / "gcmc"
        out_prefix.parent.mkdir(parents=True, exist_ok=True)
        accept_log = work / "out" / "acceptance.jsonl"

        inp = work / "test.inp"
        _write_inp(
            inp,
            f"""
inp_units:gcmc_gpu
random_seed:123
fragitp:{itp}
fragname:NA
fragconc:{conc_m}
fragmuex:{muex_kcal_mol}

box_size:10.0 10.0 10.0
cutoff:4.0
temperature:{temperature_k}
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

        rec = _first_accept_record(accept_log, move="insertion", species="NA")
        return rec

    rec_mu0 = run_once(0.0)
    rec_mu1 = run_once(-1.0)

    # Ideal-gas-ish: in an empty box, insertion should have no interaction energy.
    assert float(rec_mu0["deltaU"]) == pytest.approx(0.0, abs=1e-12)
    assert float(rec_mu1["deltaU"]) == pytest.approx(0.0, abs=1e-12)
    assert int(rec_mu0["nBefore"]) == 0
    assert int(rec_mu1["nBefore"]) == 0
    assert float(rec_mu0["vBox"]) == pytest.approx(1.0, abs=1e-12)
    assert float(rec_mu1["vBox"]) == pytest.approx(1.0, abs=1e-12)
    assert float(rec_mu0["cavityFraction"]) == pytest.approx(1.0, abs=1e-12)
    assert float(rec_mu1["cavityFraction"]) == pytest.approx(1.0, abs=1e-12)
    assert float(rec_mu0["rosenbluthWeight"]) == pytest.approx(1.0, abs=1e-12)
    assert float(rec_mu1["rosenbluthWeight"]) == pytest.approx(1.0, abs=1e-12)
    assert float(rec_mu0["proposalRatio"]) == pytest.approx(1.0, abs=1e-12)
    assert float(rec_mu1["proposalRatio"]) == pytest.approx(1.0, abs=1e-12)

    # Activity should follow z = conc * (M->nm^-3) * exp(beta*mu).
    # (The implementation uses a rounded conversion constant, so keep tolerances loose for absolute z.)
    beta = 1.0 / (8.314e-3 * temperature_k)  # mol/kJ
    expected_z0 = conc_m * 0.6022
    assert float(rec_mu0["mu"]) == pytest.approx(0.0, abs=1e-12)
    assert float(rec_mu0["z"]) == pytest.approx(expected_z0, rel=1e-4, abs=1e-12)

    expected_mu1_kj = -4.184  # -1.0 kcal/mol -> kJ/mol
    assert float(rec_mu1["mu"]) == pytest.approx(expected_mu1_kj, rel=1e-6, abs=1e-6)
    expected_z1 = float(rec_mu0["z"]) * math.exp(beta * expected_mu1_kj)
    assert float(rec_mu1["z"]) == pytest.approx(expected_z1, rel=1e-6, abs=1e-12)

    # In this limit, pAcc should reduce to z*V/(n+1). With V=1 and n=0, pAcc==z.
    assert float(rec_mu0["pAcc"]) == pytest.approx(float(rec_mu0["z"]), rel=1e-12, abs=1e-12)
    assert float(rec_mu1["pAcc"]) == pytest.approx(float(rec_mu1["z"]), rel=1e-12, abs=1e-12)

    # And changing μ should scale pAcc by exp(beta*Δμ) (and here Δμ == μ1 since μ0==0).
    pacc_ratio = float(rec_mu1["pAcc"]) / float(rec_mu0["pAcc"])
    assert float(rec_mu0["pAcc"]) < 0.2
    assert float(rec_mu1["pAcc"]) < float(rec_mu0["pAcc"])
    assert pacc_ratio == pytest.approx(math.exp(beta * float(rec_mu1["mu"])), rel=1e-6, abs=1e-12)

def test_dump_accept_log_consistent_with_final_pdb_counts(gcmc_cpu, test_data_dir, temp_dir):
    """
    No-cheating consistency check: the acceptance log 'accepted' flags must match the final system state.

    Run with a single-atom fragment (NA) and only insertion/deletion moves enabled.

    Verify: final molecule count (from *_final.pdb) == (#accepted insertions - #accepted deletions).
    """
    work = Path(temp_dir) / "accept_log_consistency"
    work.mkdir(parents=True, exist_ok=True)

    itp = test_data_dir / "charmm36.ff" / "mol" / "na.itp"
    assert itp.exists()

    out_prefix = work / "out" / "gcmc"
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    accept_log = work / "out" / "acceptance.jsonl"

    inp = work / "test.inp"
    _write_inp(
        inp,
        f"""
inp_units:gcmc_gpu
random_seed:999
fragitp:{itp}
fragname:NA
fragconc:55.0
fragmuex:0.0

box_size:10.0 10.0 10.0
cutoff:4.0
temperature:300.0
moves_per_step:1
mcsteps:20
nprint:1000
mc_move_prob:1 1 0 0
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

    records = [json.loads(line) for line in accept_log.read_text().splitlines() if line.strip()]
    na_records = [r for r in records if str(r.get("species", "")).strip().upper() == "NA"]
    assert na_records

    assert any(r.get("move") == "insertion" for r in na_records)
    assert any(r.get("move") == "deletion" for r in na_records)
    assert any(bool(r.get("accepted")) for r in na_records if r.get("move") == "insertion")

    accepted_insert = sum(1 for r in na_records if r.get("move") == "insertion" and bool(r.get("accepted")))
    accepted_delete = sum(1 for r in na_records if r.get("move") == "deletion" and bool(r.get("accepted")))
    expected_final = accepted_insert - accepted_delete
    assert expected_final >= 0

    out_pdb = Path(f"{out_prefix}_final.pdb")
    assert out_pdb.exists()
    final_count = _count_residues_by_resname(out_pdb, "NA")

    assert final_count == expected_final


def test_inp_random_seed_used_when_cli_missing(gcmc_cpu, test_data_dir, temp_dir):
    """random_seed/seed in INP should make runs reproducible when --seed is not provided."""
    itp = test_data_dir / "charmm36.ff" / "mol" / "sol.itp"
    assert itp.exists()

    def run_once(work: Path, seed: int) -> tuple[tuple[str, str, int, float, float, float], ...]:
        work.mkdir(parents=True, exist_ok=True)
        out_prefix = work / "out" / "gcmc"
        out_prefix.parent.mkdir(parents=True, exist_ok=True)

        inp = work / "test.inp"
        _write_inp(
            inp,
            f"""
inp_units:gcmc_gpu
random_seed:{seed}
fragitp:{itp}
box_size:30.0 30.0 30.0
cutoff:12.0
temperature:300.0
moves_per_step:1
mcsteps:10
nprint:5
fragname:SOL
fragconc:55.0
fragmuex:50.0
attempt_prob_ins:1.0
attempt_prob_del:0.0
attempt_prob_trn:0.0
attempt_prob_rot:0.0
""",
        )

        result = _run_gcmc_cpu(
            gcmc_cpu,
            workdir=work,
            inp=inp,
            out_prefix=out_prefix,
        )
        assert result.returncode == 0, result.stdout + result.stderr

        out_pdb = Path(f"{out_prefix}_final.pdb")
        assert out_pdb.exists()
        atoms: list[tuple[str, str, int, float, float, float]] = []
        for line in out_pdb.read_text().splitlines():
            if not line.startswith(("ATOM", "HETATM")):
                continue
            atom_name = line[12:16].strip().upper()
            resname = line[17:20].strip().upper()
            resid = int(line[22:26])
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            atoms.append((resname, atom_name, resid, round(x, 3), round(y, 3), round(z, 3)))
        return tuple(sorted(atoms))

    run1 = Path(temp_dir) / "seed_run1"
    run2 = Path(temp_dir) / "seed_run2"
    assert run_once(run1, 123) == run_once(run2, 123)

    # And with a different seed, we should almost certainly get a different trajectory/output.
    run3 = Path(temp_dir) / "seed_run3"
    run4 = Path(temp_dir) / "seed_run4"
    assert run_once(run3, 123) != run_once(run4, 124)


def test_inp_units_gcmc_gpu_converts_grid_dx_and_cutoffs_and_target_volume(
    gcmc_cpu, test_data_dir, temp_dir
):
    """
    Verify unit conversions for less-visible keys via a stable JSON dump (no log matching):
    - grid_dx (Å -> nm)
    - energy_cutoff_frag (Å -> nm) and mirrored cutoff
    - target_volume (Å^3 -> nm^3)
    """
    work = Path(temp_dir) / "units_dump_params"
    work.mkdir(parents=True, exist_ok=True)

    itp = test_data_dir / "charmm36.ff" / "mol" / "na.itp"
    assert itp.exists()

    out_prefix = work / "out" / "gcmc"
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    params_json = work / "out" / "params.json"

    inp = work / "test.inp"
    _write_inp(
        inp,
        f"""
inp_units:gcmc_gpu
fragitp:{itp}
fragname:NA
fragconc:1.0
fragmuex:0.0

box_size:10.0 10.0 10.0
cutoff:12.0
grid_dx:1.0
energy_cutoff_frag:12.0
target_volume:1000.0

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
        extra_args=["--dump-params", str(params_json)],
    )
    assert result.returncode == 0, result.stdout + result.stderr

    assert params_json.exists(), "Expected --dump-params to create params.json"
    params = json.loads(params_json.read_text())

    # Internal units are nm / kJ/mol.
    assert float(params["space"]["grid_spacing_nm"]) == pytest.approx(0.1, abs=1e-6)
    assert float(params["space"]["cutoff_nm"]) == pytest.approx(1.2, abs=1e-6)
    assert float(params["energy"]["fragment_cutoff_nm"]) == pytest.approx(1.2, abs=1e-6)
    assert float(params["space"]["target_volume_nm3"]) == pytest.approx(1.0, abs=1e-6)


@pytest.mark.parametrize(
    ("region_spec_angstrom", "expected_volume_nm3"),
    [
        # Use relatively large regions to avoid occasional rejection due to hard region constraint.
        ("sphere 15.0 15.0 15.0 15.0", (4.0 / 3.0) * 3.141592653589793 * (1.5**3)),  # r=15Å=1.5nm
        ("box 2.0 2.0 2.0 28.0 28.0 28.0", 2.6**3),  # 0.2..2.8 nm inside a 3 nm box
        ("cylinder 15.0 15.0 15.0 15.0 20.0 z", 3.141592653589793 * (1.5**2) * 2.0),
    ],
)
def test_inp_units_gcmc_gpu_gcmc_region_numeric_conversion_affects_volume(
    gcmc_cpu, test_data_dir, temp_dir, region_spec_angstrom, expected_volume_nm3
):
    """
    gcmc_region numeric values are in Å for gcmc_gpu inputs, but RegionConstraint expects nm.
    Validate end-to-end by asserting the effective acceptance volume (vEff) equals the converted region volume.
    """
    work = Path(temp_dir) / "units_region" / region_spec_angstrom.split()[0]
    work.mkdir(parents=True, exist_ok=True)

    itp = test_data_dir / "charmm36.ff" / "mol" / "sol.itp"
    assert itp.exists()

    out_prefix = work / "out" / "gcmc"
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    accept_log = work / "out" / "acceptance.jsonl"

    inp = work / "test.inp"
    _write_inp(
        inp,
        f"""
inp_units:gcmc_gpu
random_seed:123
fragitp:{itp}
fragname:SOL
fragmuex:0.0

box_size:30.0 30.0 30.0
cutoff:12.0
gcmc_region:{region_spec_angstrom}

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

    records = [
        json.loads(line)
        for line in accept_log.read_text().splitlines()
        if line.strip()
    ]
    assert records, "Acceptance log unexpectedly empty"

    ins = next(
        (r for r in records if r.get("move") == "insertion" and str(r.get("species", "")).upper() == "SOL"),
        None,
    )
    assert ins is not None, f"Expected a SOL insertion record, got: {records[:3]}"
    assert float(ins["cavityFraction"]) == pytest.approx(1.0, abs=1e-12)
    assert float(ins["vEff"]) == pytest.approx(expected_volume_nm3, rel=1e-6, abs=1e-6)


def _first_accept_record(path: Path, *, move: str, species: str) -> dict:
    records = [json.loads(line) for line in path.read_text().splitlines() if line.strip()]
    want_move = move.strip().lower()
    want_species = species.strip().upper()
    for r in records:
        if str(r.get("move", "")).strip().lower() != want_move:
            continue
        if str(r.get("species", "")).strip().upper() != want_species:
            continue
        return r
    raise AssertionError(f"No {move}/{species} record found in {path}; first records: {records[:3]}")


def test_nbar_volume_based_mode_sets_activity_from_concentration(gcmc_cpu, test_data_dir, temp_dir):
    """
    Default (volume-based) mode: activity should be conc(M) * 0.602214... (molecules/nm^3).
    Use a 1 nm^3 box so activity is directly observable from acceptance log 'z'.
    """
    work = Path(temp_dir) / "nbar_volume_based"
    work.mkdir(parents=True, exist_ok=True)

    itp = test_data_dir / "charmm36.ff" / "mol" / "na.itp"
    assert itp.exists()

    out_prefix = work / "out" / "gcmc"
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    accept_log = work / "out" / "acceptance.jsonl"

    inp = work / "test.inp"
    _write_inp(
        inp,
        f"""
inp_units:gcmc_gpu
random_seed:7
fragitp:{itp}
fragname:NA
fragconc:1.0
fragmuex:0.0

box_size:10.0 10.0 10.0
cutoff:12.0
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
    )
    assert result.returncode == 0, result.stdout + result.stderr

    rec = _first_accept_record(accept_log, move="insertion", species="NA")
    assert float(rec["z"]) == pytest.approx(0.6022, rel=1e-4, abs=1e-4)


def test_nbar_const_water_nbar_scales_all_fragments_by_fragconc(gcmc_cpu, test_data_dir, temp_dir):
    """
    const_water_nbar mode (gcmc_gpu semantics): for any fragment i,
      nbar_i = const_water_nbar / water_density * fragconc_i
      activity_i = nbar_i / V * exp(beta*muex_i)
    This must apply to non-water fragments too (not just SOL/WAT).
    """
    work = Path(temp_dir) / "nbar_const_water_all_species"
    work.mkdir(parents=True, exist_ok=True)

    itp = test_data_dir / "charmm36.ff" / "mol" / "na.itp"
    assert itp.exists()

    out_prefix = work / "out" / "gcmc"
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    accept_log = work / "out" / "acceptance.jsonl"

    inp = work / "test.inp"
    _write_inp(
        inp,
        f"""
inp_units:gcmc_gpu
random_seed:11
use_const_water_nbar:55
fragitp:{itp}
fragname:NA
fragconc:1.0
fragmuex:0.0

box_size:10.0 10.0 10.0
cutoff:12.0
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
    )
    assert result.returncode == 0, result.stdout + result.stderr

    rec = _first_accept_record(accept_log, move="insertion", species="NA")
    # V = 1 nm^3, water_density defaults to 55 M, muex=0 => activity = 55/55 * 1 / 1 = 1.0
    assert float(rec["z"]) == pytest.approx(1.0, abs=1e-12)


def test_nbar_number_water_nbar_updates_activity_after_first_insertion(gcmc_cpu, test_data_dir, temp_dir):
    """
    number_water_nbar mode: after the first accepted water insertion, waterCount==1 so
    activity for water should update to z = 1/V (with V=1 nm^3 here).
    """
    work = Path(temp_dir) / "nbar_number_water"
    work.mkdir(parents=True, exist_ok=True)

    itp = test_data_dir / "charmm36.ff" / "mol" / "sol.itp"
    assert itp.exists()

    out_prefix = work / "out" / "gcmc"
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    accept_log = work / "out" / "acceptance.jsonl"

    inp = work / "test.inp"
    _write_inp(
        inp,
        f"""
inp_units:gcmc_gpu
random_seed:123
use_number_water_nbar:yes
fragitp:{itp}
fragname:SOL
fragconc:55.0
fragmuex:0.0

box_size:10.0 10.0 10.0
cutoff:12.0
temperature:300.0
moves_per_step:1
mcsteps:2
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

    records = [json.loads(line) for line in accept_log.read_text().splitlines() if line.strip()]
    sol_ins = [r for r in records if r.get("move") == "insertion" and str(r.get("species", "")).upper() == "SOL"]
    assert len(sol_ins) >= 2, f"Expected >=2 SOL insertions, got {len(sol_ins)}"

    # First move: before any waters exist, leave base activity from concentration.
    assert float(sol_ins[0]["z"]) == pytest.approx(55.0 * 0.6022, rel=1e-4, abs=1e-4)

    # Second move: after 1 accepted insertion, number-water nbar sets activity to N/V (N=1, V=1 nm^3).
    assert float(sol_ins[1]["z"]) == pytest.approx(1.0, abs=1e-12)


def test_gcmc_gpu_active_and_muex_files_are_written(gcmc_cpu, test_data_dir, temp_dir):
    """
    gcmc_gpu-style outputs should be produced without relying on stdout:
    - active_<frag>.dat
    - muex_<frag>.dat (kcal/mol)
    Files should be written next to the chosen --prefix path.
    """
    work = Path(temp_dir) / "outputs_active_muex"
    work.mkdir(parents=True, exist_ok=True)

    itp = test_data_dir / "charmm36.ff" / "mol" / "sol.itp"
    assert itp.exists()

    out_prefix = work / "out" / "gcmc"
    out_dir = out_prefix.parent
    out_dir.mkdir(parents=True, exist_ok=True)

    inp = work / "test.inp"
    _write_inp(
        inp,
        f"""
inp_units:gcmc_gpu
random_seed:321
fragitp:{itp}
fragname:SOL
fragconc:55.0
fragmuex:-1.0

box_size:10.0 10.0 10.0
cutoff:12.0
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
        timeout=30,
    )
    assert result.returncode == 0, result.stdout + result.stderr

    active = out_dir / "active_SOL.dat"
    muex = out_dir / "muex_SOL.dat"
    assert active.exists()
    assert muex.exists()

    active_lines = [line.strip() for line in active.read_text().splitlines() if line.strip()]
    muex_lines = [line.strip() for line in muex.read_text().splitlines() if line.strip()]
    assert active_lines, "active_<frag>.dat unexpectedly empty"
    assert muex_lines, "muex_<frag>.dat unexpectedly empty"

    active_last = int(active_lines[-1])
    assert active_last >= 0

    out_pdb = Path(f"{out_prefix}_final.pdb")
    assert out_pdb.exists()
    sol_resids = set()
    for line in out_pdb.read_text().splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        resname = line[17:20].strip().upper()
        if resname != "SOL":
            continue
        sol_resids.add(int(line[22:26]))
    assert active_last == len(sol_resids)

    muex_last = float(muex_lines[-1])
    assert muex_last == pytest.approx(-1.0, abs=0.02)

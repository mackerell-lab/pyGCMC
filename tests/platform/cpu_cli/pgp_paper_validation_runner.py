"""
PGP paper validation runner for B1/B2/B3/B4 evidence package.

This script implements the experiment workflow described in docs/025:
- B1: Production-style T4 lysozyme occupancy validation (Mode C vs Mode E)
- B2: Multi-charge mesh-self diagnostic with a rigid TIP3P-like water probe
- B3: Mesh sensitivity diagnostic with a single-charge probe (power-of-two meshes)
- B4: Convergence trace extraction from B1 occupancy trajectories

The script intentionally writes machine-readable outputs (JSON/CSV) for paper tables.
"""

from __future__ import annotations

import argparse
import json
import math
import random
import statistics
import subprocess
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np

try:
    import pygcmc
except ImportError as import_error:
    raise RuntimeError(
        "pygcmc import failed. Build bindings first and set "
        "PYTHONPATH=pygcmc_dev/build/modules/bindings"
    ) from import_error


SCRIPT_PATH = Path(__file__).resolve()
PYGCMC_DEV_ROOT = SCRIPT_PATH.parents[3]
WORKSPACE_ROOT = SCRIPT_PATH.parents[4]
TEST_DATA_ROOT = PYGCMC_DEV_ROOT / "tests" / "data"
DEFAULT_OUT_DIR = WORKSPACE_ROOT / "tmp" / "pgp_paper_validation"


@dataclass
class B1SeedResult:
    mode: str
    seed: int
    occupancy_trace: list[int]
    production_trace: list[int]
    insertion_attempts: int
    insertion_accepts: int
    insertion_accept_rate: float
    deletion_attempts: int
    deletion_accepts: int
    deletion_accept_rate: float


def _write_text(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content.strip() + "\n")


def _parse_int_list_csv(values: str) -> list[int]:
    parsed: list[int] = []
    for token in values.split(","):
        clean = token.strip()
        if not clean:
            continue
        parsed.append(int(clean))
    if not parsed:
        raise ValueError("At least one integer is required.")
    return parsed


def _parse_float_quad(values: str) -> tuple[float, float, float, float]:
    tokens = values.replace(",", " ").split()
    if len(tokens) != 4:
        raise ValueError(
            f"Expected exactly 4 values for move probabilities, got {len(tokens)}."
        )
    val0, val1, val2, val3 = (float(token) for token in tokens)
    return (val0, val1, val2, val3)


def _find_gcmc_cpu(explicit_path: str | None) -> Path:
    if explicit_path:
        candidate = Path(explicit_path).resolve()
        if candidate.exists():
            return candidate
        raise FileNotFoundError(f"gcmc_cpu not found: {candidate}")

    candidate = PYGCMC_DEV_ROOT / "build" / "bin" / "gcmc_cpu"
    if candidate.exists():
        return candidate

    which_result = subprocess.run(
        ["which", "gcmc_cpu"], capture_output=True, text=True, check=False
    )
    if which_result.returncode == 0:
        return Path(which_result.stdout.strip())

    raise FileNotFoundError(
        "gcmc_cpu executable not found. Build target gcmc_cpu first."
    )


def _sanitize_json_line(raw_line: str) -> str:
    sanitized = raw_line
    sanitized = sanitized.replace("-nan", "0.0")
    sanitized = sanitized.replace("nan", "0.0")
    sanitized = sanitized.replace("-inf", "-1.0e308")
    sanitized = sanitized.replace("inf", "1.0e308")
    return sanitized


def _load_jsonl_records(path: Path) -> list[dict]:
    records: list[dict] = []
    for raw_line in path.read_text().splitlines():
        if not raw_line.strip():
            continue
        try:
            records.append(json.loads(raw_line))
        except json.JSONDecodeError:
            records.append(json.loads(_sanitize_json_line(raw_line)))
    return records


def _find_active_trace_file(run_output_dir: Path, fragment_name: str) -> Path:
    want = fragment_name.strip().lower()
    candidates = sorted(run_output_dir.glob("active_*.dat"))
    for candidate in candidates:
        stem = candidate.stem  # active_<name>
        suffix = stem[len("active_") :] if stem.startswith("active_") else stem
        if suffix.lower() == want:
            return candidate
    if candidates:
        return candidates[0]
    raise FileNotFoundError(
        f"No active_*.dat file found in {run_output_dir}. "
        "Set nprint/print-freq and ensure stats writing is enabled."
    )


def _read_integer_trace(path: Path) -> list[int]:
    values: list[int] = []
    for line in path.read_text().splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        values.append(int(float(stripped)))
    return values


def _acceptance_counts(records: Iterable[dict], move: str, species: str) -> tuple[int, int]:
    move_key = move.strip().lower()
    species_key = species.strip().upper()
    attempts = 0
    accepts = 0
    for record in records:
        if str(record.get("move", "")).strip().lower() != move_key:
            continue
        if str(record.get("species", "")).strip().upper() != species_key:
            continue
        attempts += 1
        if bool(record.get("accepted", False)):
            accepts += 1
    return attempts, accepts


def _rate_from_records(records: Iterable[dict], move: str, species: str) -> float:
    attempts, accepts = _acceptance_counts(records, move, species)
    if attempts == 0:
        return 0.0
    return accepts / float(attempts)


def _build_b1_lysozyme_inp(
    *,
    mode: str,
    seed: int,
    total_steps: int,
    print_freq: int,
    sol_fragconc: float,
    sol_muex: float,
    cosolvent_mctime_weight: float,
    sol_mctime_weight: float,
    use_cavity_bias: bool,
    use_conf_bias: bool,
    use_insdel_frac: bool,
    insdel_frac: float,
    move_probabilities: tuple[float, float, float, float],
    use_number_water_nbar: bool,
    energy_cutoff_frag: float | None,
    energy_cutoff_prot: float | None,
    num_conf_bias_trial: int,
    run_dir: Path,
) -> Path:
    ff_dir = TEST_DATA_ROOT / "charmm36.ff"
    lysozyme_dir = TEST_DATA_ROOT / "gcmc_opencl_examples" / "lysozyme"
    inp_path = run_dir / "run.inp"
    if use_insdel_frac:
        move_control_line = f"insdel_frac:{insdel_frac:.6f}"
    else:
        move_control_line = (
            "mc_move_prob:"
            f"{move_probabilities[0]:.6f} {move_probabilities[1]:.6f} "
            f"{move_probabilities[2]:.6f} {move_probabilities[3]:.6f}"
        )

    optional_lines: list[str] = []
    if energy_cutoff_frag is not None:
        optional_lines.append(f"energy_cutoff_frag:{energy_cutoff_frag:.6f}")
    if energy_cutoff_prot is not None:
        optional_lines.append(f"energy_cutoff_prot:{energy_cutoff_prot:.6f}")
    optional_lines.append(
        f"use_number_water_nbar:{'yes' if use_number_water_nbar else 'no'}"
    )
    if use_conf_bias:
        optional_lines.append(f"num_conf_bias_trial:{max(1, int(num_conf_bias_trial))}")
    optional_block = "\n".join(optional_lines)

    _write_text(
        inp_path,
        f"""
random_seed:{seed}
energy_method:{mode}
par:{ff_dir / "ffnonbonded.itp"}
par:{ff_dir / "silcs.itp"}
par:{ff_dir / "nbfix.itp"}
fragitp:{ff_dir / "mol" / "benx.itp"}
fragitp:{ff_dir / "mol" / "prpx.itp"}
fragitp:{ff_dir / "mol" / "dmee.itp"}
fragitp:{ff_dir / "mol" / "meoh.itp"}
fragitp:{ff_dir / "mol" / "form.itp"}
fragitp:{ff_dir / "mol" / "imia.itp"}
fragitp:{ff_dir / "mol" / "acey.itp"}
fragitp:{ff_dir / "mol" / "mamy.itp"}
fragitp:{ff_dir / "mol" / "sol.itp"}
atomtypes:{ff_dir / "atomtypes.atp"}
monomerdir:{ff_dir / "mol"}

top:{lysozyme_dir / "181L_apo_silcs.1.top"}
pdb:{lysozyme_dir / "181L_apo_silcs.1.equil.rec.pdb"}
protitp:{lysozyme_dir / "181L_apo_silcs.1.top"}

box_size:36.736 40.850 49.379
gc_center:33.368 35.425 39.690
sys_center:33.368 35.425 39.690
cutoff:12.0
temperature:300.0
{optional_block}

fragname:benx prpx dmee meoh form imia acey mamy sol
fragconc:0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 {sol_fragconc:.6f}
fragmuex:-0.79 1.96 -1.79 -5.36 -10.92 -14.18 -97.31 -68.49 {sol_muex:.6f}
mctime:{cosolvent_mctime_weight:.6f} {cosolvent_mctime_weight:.6f} {cosolvent_mctime_weight:.6f} {cosolvent_mctime_weight:.6f} {cosolvent_mctime_weight:.6f} {cosolvent_mctime_weight:.6f} {cosolvent_mctime_weight:.6f} {cosolvent_mctime_weight:.6f} {sol_mctime_weight:.6f}

moves_per_step:1
mcsteps:{total_steps}
nprint:{print_freq}
{move_control_line}
max_translation:3.0
max_rotation:20.0

use_cavity_bias:{"yes" if use_cavity_bias else "no"}
use_conf_bias:{"yes" if use_conf_bias else "no"}
map_generation:no
""",
    )
    return inp_path


def _run_subprocess(
    args: list[str], *, cwd: Path, timeout_sec: int
) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        args,
        cwd=str(cwd),
        capture_output=True,
        text=True,
        timeout=timeout_sec,
        check=False,
    )


def _run_b1_seed(
    *,
    gcmc_cpu: Path,
    out_dir: Path,
    mode: str,
    seed: int,
    warmup_steps: int,
    production_steps: int,
    print_freq: int,
    sol_fragconc: float,
    sol_muex: float,
    cosolvent_mctime_weight: float,
    sol_mctime_weight: float,
    use_cavity_bias: bool,
    use_conf_bias: bool,
    use_insdel_frac: bool,
    insdel_frac: float,
    move_probabilities: tuple[float, float, float, float],
    use_number_water_nbar: bool,
    energy_cutoff_frag: float | None,
    energy_cutoff_prot: float | None,
    num_conf_bias_trial: int,
    max_molecules_per_type: int,
    timeout_sec: int,
) -> B1SeedResult:
    run_dir = out_dir / "b1" / mode / f"seed_{seed}"
    run_output_dir = run_dir / "out"
    run_output_dir.mkdir(parents=True, exist_ok=True)

    total_steps = warmup_steps + production_steps
    inp_path = _build_b1_lysozyme_inp(
        mode=mode,
        seed=seed,
        total_steps=total_steps,
        print_freq=print_freq,
        sol_fragconc=sol_fragconc,
        sol_muex=sol_muex,
        cosolvent_mctime_weight=cosolvent_mctime_weight,
        sol_mctime_weight=sol_mctime_weight,
        use_cavity_bias=use_cavity_bias,
        use_conf_bias=use_conf_bias,
        use_insdel_frac=use_insdel_frac,
        insdel_frac=insdel_frac,
        move_probabilities=move_probabilities,
        use_number_water_nbar=use_number_water_nbar,
        energy_cutoff_frag=energy_cutoff_frag,
        energy_cutoff_prot=energy_cutoff_prot,
        num_conf_bias_trial=num_conf_bias_trial,
        run_dir=run_dir,
    )
    out_prefix = run_output_dir / "gcmc"
    accept_log = run_output_dir / "acceptance.jsonl"
    params_json = run_output_dir / "params.json"

    result = _run_subprocess(
        [
            str(gcmc_cpu),
            "--inp",
            str(inp_path),
            "--prefix",
            str(out_prefix),
            "--dump-accept",
            str(accept_log),
            "--dump-params",
            str(params_json),
            "--print-freq",
            str(print_freq),
            "--traj-freq",
            str(max(total_steps + 1, 1000000)),
            "--checkpoint-freq",
            "0",
            "--max-molecules-per-type",
            str(max_molecules_per_type),
        ],
        cwd=run_dir,
        timeout_sec=timeout_sec,
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"B1 run failed for mode={mode}, seed={seed}\n"
            f"STDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"
        )

    records = _load_jsonl_records(accept_log)
    active_path = _find_active_trace_file(run_output_dir, "sol")
    occupancy_trace = _read_integer_trace(active_path)
    if not occupancy_trace:
        raise RuntimeError(f"Empty occupancy trace in {active_path}")

    warmup_points = min(len(occupancy_trace), warmup_steps // max(print_freq, 1))
    production_trace = occupancy_trace[warmup_points:]
    if not production_trace:
        production_trace = occupancy_trace[-1:]

    insertion_attempts, insertion_accepts = _acceptance_counts(records, "insertion", "SOL")
    deletion_attempts, deletion_accepts = _acceptance_counts(records, "deletion", "SOL")
    insertion_rate = 0.0 if insertion_attempts == 0 else insertion_accepts / float(insertion_attempts)
    deletion_rate = 0.0 if deletion_attempts == 0 else deletion_accepts / float(deletion_attempts)
    return B1SeedResult(
        mode=mode,
        seed=seed,
        occupancy_trace=occupancy_trace,
        production_trace=production_trace,
        insertion_attempts=insertion_attempts,
        insertion_accepts=insertion_accepts,
        insertion_accept_rate=insertion_rate,
        deletion_attempts=deletion_attempts,
        deletion_accepts=deletion_accepts,
        deletion_accept_rate=deletion_rate,
    )


def _mean_confidence_interval95(values: list[float]) -> dict:
    if not values:
        return {"mean": 0.0, "std": 0.0, "ci95_half_width": 0.0}
    mean_value = float(statistics.mean(values))
    if len(values) == 1:
        return {"mean": mean_value, "std": 0.0, "ci95_half_width": 0.0}
    std_value = float(statistics.stdev(values))
    ci_half_width = 1.96 * std_value / math.sqrt(len(values))
    return {
        "mean": mean_value,
        "std": std_value,
        "ci95_half_width": float(ci_half_width),
    }


def _ks_2samp(values_a: list[int], values_b: list[int]) -> dict:
    if not values_a or not values_b:
        return {"statistic": 0.0, "pvalue": 1.0}

    sorted_a = np.sort(np.asarray(values_a, dtype=float))
    sorted_b = np.sort(np.asarray(values_b, dtype=float))
    all_values = np.sort(np.concatenate([sorted_a, sorted_b]))

    cdf_a = np.searchsorted(sorted_a, all_values, side="right") / float(sorted_a.size)
    cdf_b = np.searchsorted(sorted_b, all_values, side="right") / float(sorted_b.size)
    statistic = float(np.max(np.abs(cdf_a - cdf_b)))

    n_a = float(sorted_a.size)
    n_b = float(sorted_b.size)
    effective_n = math.sqrt(n_a * n_b / (n_a + n_b))
    lam = (effective_n + 0.12 + 0.11 / effective_n) * statistic

    if lam <= 0.0:
        return {"statistic": statistic, "pvalue": 1.0}

    series_sum = 0.0
    for k_index in range(1, 200):
        term = 2.0 * ((-1.0) ** (k_index - 1)) * math.exp(
            -2.0 * (k_index**2) * (lam**2)
        )
        series_sum += term
        if abs(term) < 1e-12:
            break
    pvalue = max(0.0, min(1.0, series_sum))
    return {"statistic": statistic, "pvalue": pvalue}


def _running_mean(values: list[float]) -> list[float]:
    running: list[float] = []
    accumulator = 0.0
    for index, value in enumerate(values, start=1):
        accumulator += value
        running.append(accumulator / float(index))
    return running


def _seedwise_trace_mean(seed_traces: list[list[int]]) -> list[float]:
    if not seed_traces:
        return []
    min_len = min(len(trace) for trace in seed_traces)
    if min_len <= 0:
        return []
    matrix = np.asarray([trace[:min_len] for trace in seed_traces], dtype=float)
    return matrix.mean(axis=0).tolist()


def _estimate_iat_ess(values: list[float], max_lag: int | None = None) -> dict:
    sample_count = len(values)
    if sample_count < 4:
        return {
            "n_samples": sample_count,
            "iat": 1.0,
            "ess": float(sample_count),
            "max_lag_used": 0,
        }

    signal = np.asarray(values, dtype=float)
    centered = signal - signal.mean()
    variance = float(np.var(centered))
    if variance <= 0.0:
        return {
            "n_samples": sample_count,
            "iat": 1.0,
            "ess": float(sample_count),
            "max_lag_used": 0,
        }

    max_lag_final = min(sample_count // 2, 5000) if max_lag is None else min(max_lag, sample_count - 1)
    tau_int = 1.0
    lag_used = 0
    for lag in range(1, max_lag_final + 1):
        numerator = float(np.dot(centered[:-lag], centered[lag:]))
        denominator = float((sample_count - lag) * variance)
        if denominator <= 0.0:
            break
        autocorr = numerator / denominator
        if autocorr <= 0.0:
            break
        tau_int += 2.0 * autocorr
        lag_used = lag

    tau_int = max(1.0, tau_int)
    ess = max(1.0, min(float(sample_count), sample_count / tau_int))
    return {
        "n_samples": sample_count,
        "iat": float(tau_int),
        "ess": float(ess),
        "max_lag_used": int(lag_used),
    }


def _trace_block_means(values: list[float], block_size: int) -> list[dict]:
    if not values:
        return []
    block_size_final = max(1, int(block_size))
    rows: list[dict] = []
    for block_index, start in enumerate(range(0, len(values), block_size_final), start=1):
        stop = min(len(values), start + block_size_final)
        window = values[start:stop]
        rows.append(
            {
                "block_index": block_index,
                "start_snapshot": start + 1,
                "end_snapshot": stop,
                "count": len(window),
                "block_mean": float(statistics.mean(window)),
            }
        )
    return rows


def _maybe_plot_convergence(
    csv_rows: list[dict], figure_path: Path, disable_plot: bool
) -> bool:
    if disable_plot:
        return False
    try:
        import matplotlib.pyplot as plt  # pylint: disable=import-outside-toplevel
    except Exception:
        return False

    modes = sorted({row["mode"] for row in csv_rows})
    plt.figure(figsize=(7.2, 4.8))
    for mode in modes:
        mode_rows = [row for row in csv_rows if row["mode"] == mode]
        plt.plot(
            [row["step"] for row in mode_rows],
            [row["running_mean"] for row in mode_rows],
            label=mode,
            linewidth=1.8,
        )
    plt.xlabel("MC step")
    plt.ylabel("Running mean occupancy")
    plt.title("B4 convergence from B1 occupancy traces")
    plt.grid(True, alpha=0.3)
    plt.legend()
    figure_path.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(figure_path, dpi=180)
    plt.close()
    return True


def _configure_pme(
    mesh_size: int,
    cutoff_nm: float,
    box_nm: list[float],
    *,
    alpha_value: float,
    spline_order: int,
    tolerance: float,
) -> None:
    mesh = [mesh_size, mesh_size, mesh_size]
    if hasattr(pygcmc, "resetPGPState"):
        pygcmc.resetPGPState()
    pygcmc.setPMEParameters(alpha_value, mesh, spline_order, tolerance)
    pygcmc.initializePMEParameters(
        cutoff_nm,
        box_nm,
        alpha_value,
        mesh,
        spline_order,
        tolerance,
    )


def _build_single_charge_state(
    box_nm: float = 5.0,
    cutoff_nm: float = 1.2,
) -> "pygcmc.MCState":
    state = pygcmc.MCState()
    state.info.box = [box_nm, box_nm, box_nm]
    state.info.cutoff = cutoff_nm

    forcefield = pygcmc.MCForceField()
    forcefield.numTotalTypes = 1
    forcefield.numMovementTypes = 1
    forcefield.maxTypes = 1
    forcefield.ljEps = [0.0]
    forcefield.ljSigma = [0.0]
    state.forcefield = forcefield

    atom = pygcmc.MCAtom()
    atom.x = box_nm / 2.0
    atom.y = box_nm / 2.0
    atom.z = box_nm / 2.0
    atom.charge = 1.0
    atom.type = 0

    residue = pygcmc.MCResidue()
    residue.active = True
    residue.fixed = False
    residue.atomStart = 0
    residue.atomCount = 1
    residue.type = 0

    movement = pygcmc.MCMovementResidueInfo()
    movement.startIndex = 0
    movement.activeCount = 1

    state.atoms = [atom]
    state.residues = [residue]
    state.movementResidues = [movement]
    state.activeAtomCount = 1
    state.activeResidueCount = 1
    return state


def _build_tip3p_probe_state(
    box_nm: float = 5.0,
    cutoff_nm: float = 1.2,
) -> "pygcmc.MCState":
    state = pygcmc.MCState()
    state.info.box = [box_nm, box_nm, box_nm]
    state.info.cutoff = cutoff_nm

    forcefield = pygcmc.MCForceField()
    forcefield.numTotalTypes = 2
    forcefield.numMovementTypes = 2
    forcefield.maxTypes = 2
    forcefield.ljEps = [0.0] * 4
    forcefield.ljSigma = [0.0] * 4
    state.forcefield = forcefield

    coordinates = [
        (box_nm / 2.0, box_nm / 2.0, box_nm / 2.0, -0.834, 0),
        (box_nm / 2.0 + 0.09572, box_nm / 2.0, box_nm / 2.0, 0.417, 1),
        (box_nm / 2.0 - 0.03199, box_nm / 2.0 + 0.09268, box_nm / 2.0, 0.417, 1),
    ]

    atoms = []
    for pos_x, pos_y, pos_z, charge, atom_type in coordinates:
        atom = pygcmc.MCAtom()
        atom.x = pos_x
        atom.y = pos_y
        atom.z = pos_z
        atom.charge = charge
        atom.type = atom_type
        atoms.append(atom)

    residue = pygcmc.MCResidue()
    residue.active = True
    residue.fixed = False
    residue.atomStart = 0
    residue.atomCount = 3
    residue.type = 0

    movement = pygcmc.MCMovementResidueInfo()
    movement.startIndex = 0
    movement.activeCount = 1

    state.atoms = atoms
    state.residues = [residue]
    state.movementResidues = [movement]
    state.activeAtomCount = 3
    state.activeResidueCount = 1
    return state


def _apply_random_translation(
    atoms: list["pygcmc.MCAtom"],
    box_nm: list[float],
    rng: random.Random,
    max_displacement_nm: float = 0.2,
) -> None:
    shift_x = (rng.random() - 0.5) * max_displacement_nm
    shift_y = (rng.random() - 0.5) * max_displacement_nm
    shift_z = (rng.random() - 0.5) * max_displacement_nm

    for atom in atoms:
        atom.x = (atom.x + shift_x) % box_nm[0]
        atom.y = (atom.y + shift_y) % box_nm[1]
        atom.z = (atom.z + shift_z) % box_nm[2]


def _center_of_mass(atoms: list["pygcmc.MCAtom"]) -> tuple[float, float, float]:
    count = float(len(atoms))
    center_x = sum(atom.x for atom in atoms) / count
    center_y = sum(atom.y for atom in atoms) / count
    center_z = sum(atom.z for atom in atoms) / count
    return center_x, center_y, center_z


def _random_unit_quaternion(rng: random.Random) -> tuple[float, float, float, float]:
    value_u1 = rng.random()
    value_u2 = rng.random()
    value_u3 = rng.random()
    return (
        math.sqrt(1.0 - value_u1) * math.sin(2.0 * math.pi * value_u2),
        math.sqrt(1.0 - value_u1) * math.cos(2.0 * math.pi * value_u2),
        math.sqrt(value_u1) * math.sin(2.0 * math.pi * value_u3),
        math.sqrt(value_u1) * math.cos(2.0 * math.pi * value_u3),
    )


def _apply_random_rotation(
    atoms: list["pygcmc.MCAtom"], box_nm: list[float], rng: random.Random
) -> None:
    center_x, center_y, center_z = _center_of_mass(atoms)
    quat_x, quat_y, quat_z, quat_w = _random_unit_quaternion(rng)

    m00 = 1.0 - 2.0 * (quat_y * quat_y + quat_z * quat_z)
    m01 = 2.0 * (quat_x * quat_y - quat_z * quat_w)
    m02 = 2.0 * (quat_x * quat_z + quat_y * quat_w)
    m10 = 2.0 * (quat_x * quat_y + quat_z * quat_w)
    m11 = 1.0 - 2.0 * (quat_x * quat_x + quat_z * quat_z)
    m12 = 2.0 * (quat_y * quat_z - quat_x * quat_w)
    m20 = 2.0 * (quat_x * quat_z - quat_y * quat_w)
    m21 = 2.0 * (quat_y * quat_z + quat_x * quat_w)
    m22 = 1.0 - 2.0 * (quat_x * quat_x + quat_y * quat_y)

    for atom in atoms:
        offset_x = atom.x - center_x
        offset_y = atom.y - center_y
        offset_z = atom.z - center_z

        rotated_x = m00 * offset_x + m01 * offset_y + m02 * offset_z
        rotated_y = m10 * offset_x + m11 * offset_y + m12 * offset_z
        rotated_z = m20 * offset_x + m21 * offset_y + m22 * offset_z

        atom.x = (center_x + rotated_x) % box_nm[0]
        atom.y = (center_y + rotated_y) % box_nm[1]
        atom.z = (center_z + rotated_z) % box_nm[2]


def _collect_mesh_self_gap_series(
    *,
    state: "pygcmc.MCState",
    mesh_size: int,
    move_kind: str,
    steps: int,
    seed: int,
    alpha_value: float,
    spline_order: int,
    tolerance: float,
) -> list[float]:
    box_nm = [float(value) for value in state.info.box]
    cutoff_nm = float(state.info.cutoff)
    _configure_pme(
        mesh_size,
        cutoff_nm,
        box_nm,
        alpha_value=alpha_value,
        spline_order=spline_order,
        tolerance=tolerance,
    )

    pme_delta_values: list[float] = []

    random_stream = random.Random(seed)
    for _ in range(steps):
        pme_before = float(pygcmc.computeMovementEnergyPME(state)[0])

        if move_kind == "translation":
            _apply_random_translation(state.atoms, box_nm, random_stream)
        elif move_kind == "rotation":
            _apply_random_rotation(state.atoms, box_nm, random_stream)
        else:
            raise ValueError(f"Unsupported move kind: {move_kind}")

        pme_after = float(pygcmc.computeMovementEnergyPME(state)[0])
        pme_delta = pme_after - pme_before
        pme_delta_values.append(pme_delta)
    return pme_delta_values


def _series_summary(values: list[float]) -> dict:
    if not values:
        return {
            "count": 0,
            "max_abs": 0.0,
            "mean_abs": 0.0,
            "std": 0.0,
            "mean": 0.0,
        }
    abs_values = [abs(value) for value in values]
    return {
        "count": len(values),
        "max_abs": float(max(abs_values)),
        "mean_abs": float(statistics.mean(abs_values)),
        "std": float(statistics.pstdev(values)),
        "mean": float(statistics.mean(values)),
    }


def run_b1_and_b4(
    *,
    gcmc_cpu: Path,
    out_dir: Path,
    seeds: list[int],
    warmup_steps: int,
    production_steps: int,
    print_freq: int,
    sol_fragconc: float,
    sol_muex: float,
    cosolvent_mctime_weight: float,
    sol_mctime_weight: float,
    use_cavity_bias: bool,
    use_conf_bias: bool,
    use_insdel_frac: bool,
    insdel_frac: float,
    move_probabilities: tuple[float, float, float, float],
    use_number_water_nbar: bool,
    energy_cutoff_frag: float | None,
    energy_cutoff_prot: float | None,
    num_conf_bias_trial: int,
    max_molecules_per_type: int,
    accept_target: int,
    report_iat_ess: bool,
    b4_block_size: int,
    timeout_sec: int,
    disable_plot: bool,
) -> dict:
    mode_to_seed_results: dict[str, list[B1SeedResult]] = {"pgp_full": [], "pgp_full_pme": []}
    for mode in ["pgp_full", "pgp_full_pme"]:
        for seed in seeds:
            mode_to_seed_results[mode].append(
                _run_b1_seed(
                    gcmc_cpu=gcmc_cpu,
                    out_dir=out_dir,
                    mode=mode,
                    seed=seed,
                    warmup_steps=warmup_steps,
                    production_steps=production_steps,
                    print_freq=print_freq,
                    sol_fragconc=sol_fragconc,
                    sol_muex=sol_muex,
                    cosolvent_mctime_weight=cosolvent_mctime_weight,
                    sol_mctime_weight=sol_mctime_weight,
                    use_cavity_bias=use_cavity_bias,
                    use_conf_bias=use_conf_bias,
                    use_insdel_frac=use_insdel_frac,
                    insdel_frac=insdel_frac,
                    move_probabilities=move_probabilities,
                    use_number_water_nbar=use_number_water_nbar,
                    energy_cutoff_frag=energy_cutoff_frag,
                    energy_cutoff_prot=energy_cutoff_prot,
                    num_conf_bias_trial=num_conf_bias_trial,
                    max_molecules_per_type=max_molecules_per_type,
                    timeout_sec=timeout_sec,
                )
            )

    mode_summary: dict[str, dict] = {}
    pooled_samples: dict[str, list[int]] = {}
    b1_acceptance_targets: dict[str, dict] = {}
    for mode, seed_results in mode_to_seed_results.items():
        pooled = [value for result in seed_results for value in result.production_trace]
        pooled_samples[mode] = pooled
        per_seed_mean = [float(statistics.mean(result.production_trace)) for result in seed_results]
        per_seed_ins = [result.insertion_accept_rate for result in seed_results]
        per_seed_del = [result.deletion_accept_rate for result in seed_results]
        per_seed_ins_accepts = [result.insertion_accepts for result in seed_results]
        per_seed_del_accepts = [result.deletion_accepts for result in seed_results]

        mode_summary[mode] = {
            "seeds": [result.seed for result in seed_results],
            "n_total_samples": len(pooled),
            "mean_ci95": _mean_confidence_interval95(per_seed_mean),
            "insertion_acceptance_ci95": _mean_confidence_interval95(per_seed_ins),
            "deletion_acceptance_ci95": _mean_confidence_interval95(per_seed_del),
            "insertion_accepts_ci95": _mean_confidence_interval95(
                [float(value) for value in per_seed_ins_accepts]
            ),
            "deletion_accepts_ci95": _mean_confidence_interval95(
                [float(value) for value in per_seed_del_accepts]
            ),
            "per_seed_mean_occupancy": per_seed_mean,
            "per_seed_acceptance_counts": [
                {
                    "seed": result.seed,
                    "insertion": {
                        "attempts": result.insertion_attempts,
                        "accepts": result.insertion_accepts,
                        "rate": result.insertion_accept_rate,
                    },
                    "deletion": {
                        "attempts": result.deletion_attempts,
                        "accepts": result.deletion_accepts,
                        "rate": result.deletion_accept_rate,
                    },
                }
                for result in seed_results
            ],
        }
        b1_acceptance_targets[mode] = {
            "accept_target_per_seed": int(accept_target),
            "per_seed_pass": [
                {
                    "seed": result.seed,
                    "passes": (
                        result.insertion_accepts >= accept_target
                        and result.deletion_accepts >= accept_target
                    ),
                    "insertion_accepts": result.insertion_accepts,
                    "deletion_accepts": result.deletion_accepts,
                }
                for result in seed_results
            ],
        }
        b1_acceptance_targets[mode]["all_seeds_pass"] = all(
            row["passes"] for row in b1_acceptance_targets[mode]["per_seed_pass"]
        )

    ks_result = _ks_2samp(pooled_samples["pgp_full"], pooled_samples["pgp_full_pme"])

    convergence_rows: list[dict] = []
    block_rows: list[dict] = []
    diagnostics: dict[str, dict] = {}
    for mode, seed_results in mode_to_seed_results.items():
        mean_trace = _seedwise_trace_mean([result.occupancy_trace for result in seed_results])
        running_trace = _running_mean(mean_trace)
        for index, (occupancy_value, running_value) in enumerate(
            zip(mean_trace, running_trace), start=1
        ):
            step = index * print_freq
            convergence_rows.append(
                {
                    "mode": mode,
                    "snapshot": index,
                    "step": step,
                    "occupancy_mean": float(occupancy_value),
                    "running_mean": float(running_value),
                }
            )
        block_means = _trace_block_means(mean_trace, b4_block_size)
        for block_row in block_means:
            block_rows.append(
                {
                    "mode": mode,
                    "block_index": block_row["block_index"],
                    "start_snapshot": block_row["start_snapshot"],
                    "end_snapshot": block_row["end_snapshot"],
                    "count": block_row["count"],
                    "block_mean": block_row["block_mean"],
                }
            )

        mode_production = [value for result in seed_results for value in result.production_trace]
        mode_diag: dict = {
            "n_production_samples": len(mode_production),
            "block_size_snapshots": int(b4_block_size),
            "n_blocks": len(block_means),
        }
        if report_iat_ess:
            mode_diag["iat_ess"] = _estimate_iat_ess([float(value) for value in mode_production])
        diagnostics[mode] = mode_diag

    convergence_csv = out_dir / "b4_convergence.csv"
    convergence_csv.parent.mkdir(parents=True, exist_ok=True)
    with convergence_csv.open("w", encoding="utf-8") as handle:
        handle.write("mode,snapshot,step,occupancy_mean,running_mean\n")
        for row in convergence_rows:
            handle.write(
                f"{row['mode']},{row['snapshot']},{row['step']},"
                f"{row['occupancy_mean']:.8f},{row['running_mean']:.8f}\n"
            )

    block_csv = out_dir / "b4_block_means.csv"
    with block_csv.open("w", encoding="utf-8") as handle:
        handle.write("mode,block_index,start_snapshot,end_snapshot,count,block_mean\n")
        for row in block_rows:
            handle.write(
                f"{row['mode']},{row['block_index']},{row['start_snapshot']},"
                f"{row['end_snapshot']},{row['count']},{row['block_mean']:.8f}\n"
            )

    plot_generated = _maybe_plot_convergence(
        convergence_rows,
        out_dir / "b4_convergence.png",
        disable_plot=disable_plot,
    )

    b1_summary = {
        "warmup_steps": warmup_steps,
        "production_steps": production_steps,
        "print_freq": print_freq,
        "sol_fragconc": float(sol_fragconc),
        "sol_muex": float(sol_muex),
        "cosolvent_mctime_weight": float(cosolvent_mctime_weight),
        "sol_mctime_weight": float(sol_mctime_weight),
        "use_cavity_bias": bool(use_cavity_bias),
        "use_conf_bias": bool(use_conf_bias),
        "use_insdel_frac": bool(use_insdel_frac),
        "insdel_frac": float(insdel_frac),
        "move_probabilities": [float(value) for value in move_probabilities],
        "use_number_water_nbar": bool(use_number_water_nbar),
        "energy_cutoff_frag": (
            None if energy_cutoff_frag is None else float(energy_cutoff_frag)
        ),
        "energy_cutoff_prot": (
            None if energy_cutoff_prot is None else float(energy_cutoff_prot)
        ),
        "num_conf_bias_trial": int(num_conf_bias_trial),
        "max_molecules_per_type": int(max_molecules_per_type),
        "mode_summary": mode_summary,
        "acceptance_target_check": b1_acceptance_targets,
        "ks_test": ks_result,
        "b4_diagnostics": diagnostics,
        "histogram_mode_c": dict(
            sorted(Counter(pooled_samples["pgp_full"]).items(), key=lambda item: item[0])
        ),
        "histogram_mode_e": dict(
            sorted(Counter(pooled_samples["pgp_full_pme"]).items(), key=lambda item: item[0])
        ),
    }

    (out_dir / "b1_occupancy_summary.json").write_text(
        json.dumps(b1_summary, indent=2, sort_keys=True) + "\n"
    )
    return {
        "b1_summary": b1_summary,
        "b4_csv": str(convergence_csv),
        "b4_block_csv": str(block_csv),
        "b4_plot_generated": plot_generated,
    }


def run_b2(
    *,
    out_dir: Path,
    mesh_sizes: list[int],
    box_nm: float,
    cutoff_nm: float,
    alpha_value: float,
    spline_order: int,
    tolerance: float,
    steps: int,
    seed: int,
) -> dict:
    rows: list[dict] = []
    for mesh_size in mesh_sizes:
        translation_state = _build_tip3p_probe_state(box_nm=box_nm, cutoff_nm=cutoff_nm)
        rotation_state = _build_tip3p_probe_state(box_nm=box_nm, cutoff_nm=cutoff_nm)
        translation_gap = _collect_mesh_self_gap_series(
            state=translation_state,
            mesh_size=mesh_size,
            move_kind="translation",
            steps=steps,
            seed=seed + mesh_size * 31,
            alpha_value=alpha_value,
            spline_order=spline_order,
            tolerance=tolerance,
        )
        rotation_gap = _collect_mesh_self_gap_series(
            state=rotation_state,
            mesh_size=mesh_size,
            move_kind="rotation",
            steps=steps,
            seed=seed + mesh_size * 31 + 7919,
            alpha_value=alpha_value,
            spline_order=spline_order,
            tolerance=tolerance,
        )
        rows.append(
            {
                "mesh_size": mesh_size,
                "translation": _series_summary(translation_gap),
                "rotation": _series_summary(rotation_gap),
            }
        )

    primary_row = rows[0]
    b2_summary = {
        "mesh_size": primary_row["mesh_size"],
        "probe": "rigid_tip3p_like_water",
        "estimator": "isolated_probe_pme_delta",
        "box_nm": float(box_nm),
        "cutoff_nm": float(cutoff_nm),
        "alpha": float(alpha_value),
        "spline_order": int(spline_order),
        "tolerance": float(tolerance),
        "translation": primary_row["translation"],
        "rotation": primary_row["rotation"],
        "mesh_rows": rows,
    }
    (out_dir / "b2_multicharge_mesh_self.json").write_text(
        json.dumps(b2_summary, indent=2, sort_keys=True) + "\n"
    )
    return b2_summary


def run_b3(
    *,
    out_dir: Path,
    mesh_sizes: list[int],
    box_nm: float,
    cutoff_nm: float,
    alpha_value: float,
    spline_order: int,
    tolerance: float,
    steps: int,
    seed: int,
) -> dict:
    rows: list[dict] = []
    for mesh_size in mesh_sizes:
        single_charge_state = _build_single_charge_state(box_nm=box_nm, cutoff_nm=cutoff_nm)
        gap_series = _collect_mesh_self_gap_series(
            state=single_charge_state,
            mesh_size=mesh_size,
            move_kind="translation",
            steps=steps,
            seed=seed + mesh_size,
            alpha_value=alpha_value,
            spline_order=spline_order,
            tolerance=tolerance,
        )
        row = {
            "requested_mesh": mesh_size,
            "summary": _series_summary(gap_series),
        }
        rows.append(row)

    monotonic_nonincreasing = True
    for idx in range(1, len(rows)):
        if rows[idx]["summary"]["max_abs"] > rows[idx - 1]["summary"]["max_abs"]:
            monotonic_nonincreasing = False
            break

    b3_summary = {
        "probe": "single_charge_translation",
        "box_nm": float(box_nm),
        "cutoff_nm": float(cutoff_nm),
        "alpha": float(alpha_value),
        "spline_order": int(spline_order),
        "tolerance": float(tolerance),
        "mesh_rows": rows,
        "max_abs_monotonic_nonincreasing": monotonic_nonincreasing,
    }
    (out_dir / "b3_mesh_sensitivity.json").write_text(
        json.dumps(b3_summary, indent=2, sort_keys=True) + "\n"
    )
    return b3_summary


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run PGP paper validation package (B1/B2/B3/B4)."
    )
    parser.add_argument(
        "--gcmc-cpu",
        default=None,
        help="Path to gcmc_cpu executable (default: pygcmc_dev/build/bin/gcmc_cpu).",
    )
    parser.add_argument(
        "--out-dir",
        default=str(DEFAULT_OUT_DIR),
        help=f"Output directory (default: {DEFAULT_OUT_DIR}).",
    )
    parser.add_argument(
        "--b1-profile",
        choices=["default", "paper_calibrated", "old_lysozyme", "old_protein"],
        default="default",
        help=(
            "B1 parameter profile. default keeps current behavior; "
            "paper_calibrated applies the validated manuscript settings; "
            "old_lysozyme/old_protein align with legacy opencl input styles."
        ),
    )
    parser.add_argument(
        "--b1-seeds",
        default="101,202,303",
        help="Comma-separated RNG seeds for B1.",
    )
    parser.add_argument(
        "--b1-warmup-steps",
        type=int,
        default=2000,
        help="Warmup steps per B1 seed.",
    )
    parser.add_argument(
        "--b1-production-steps",
        type=int,
        default=8000,
        help="Production steps per B1 seed.",
    )
    parser.add_argument(
        "--b1-print-freq",
        type=int,
        default=20,
        help="Print frequency for occupancy trace extraction in B1.",
    )
    parser.add_argument(
        "--b1-sol-conc",
        type=float,
        default=55.0,
        help="Override SOL concentration in B1 input deck.",
    )
    parser.add_argument(
        "--b1-sol-muex",
        type=float,
        default=-5.60,
        help="Override SOL excess chemical potential in B1 input deck.",
    )
    parser.add_argument(
        "--b1-cosolvent-mctime-weight",
        type=float,
        default=1.0,
        help="MC-time selection weight for each non-SOL fragment in B1.",
    )
    parser.add_argument(
        "--b1-sol-mctime-weight",
        type=float,
        default=1.0,
        help="MC-time selection weight for SOL in B1.",
    )
    parser.add_argument(
        "--b1-use-cavity-bias",
        choices=["yes", "no"],
        default="no",
        help="Enable cavity-biased insertions in B1 input deck.",
    )
    parser.add_argument(
        "--b1-use-conf-bias",
        choices=["yes", "no"],
        default="no",
        help="Enable configuration-bias insertions in B1 input deck.",
    )
    parser.add_argument(
        "--b1-num-conf-bias-trial",
        type=int,
        default=10,
        help="num_conf_bias_trial written when configuration bias is enabled.",
    )
    parser.add_argument(
        "--b1-use-insdel-frac",
        choices=["yes", "no"],
        default="no",
        help="Use insdel_frac instead of explicit mc_move_prob in B1 input deck.",
    )
    parser.add_argument(
        "--b1-insdel-frac",
        type=float,
        default=0.8,
        help="Insertion/deletion fraction used when --b1-use-insdel-frac=yes.",
    )
    parser.add_argument(
        "--b1-mc-move-prob",
        default="0.5,0.5,0.0,0.0",
        help="Move probabilities (ins,del,trn,rot) when --b1-use-insdel-frac=no.",
    )
    parser.add_argument(
        "--b1-use-number-water-nbar",
        choices=["yes", "no"],
        default="no",
        help="Write use_number_water_nbar in B1 input deck.",
    )
    parser.add_argument(
        "--b1-energy-cutoff-frag",
        type=float,
        default=None,
        help="Optional energy_cutoff_frag (legacy-style; in input units).",
    )
    parser.add_argument(
        "--b1-energy-cutoff-prot",
        type=float,
        default=None,
        help="Optional energy_cutoff_prot (legacy-style; in input units).",
    )
    parser.add_argument(
        "--b1-max-molecules-per-type",
        type=int,
        default=-1,
        help=(
            "Pass-through value for gcmc_cpu --max-molecules-per-type in B1 runs. "
            "Use -1 to disable hard cap."
        ),
    )
    parser.add_argument(
        "--b1-accept-target",
        type=int,
        default=30,
        help="Target accepted insertions/deletions per seed for B1 completeness checks.",
    )
    parser.add_argument(
        "--b2-steps",
        type=int,
        default=400,
        help="Number of translation/rotation moves for B2 diagnostics.",
    )
    parser.add_argument(
        "--b3-steps",
        type=int,
        default=400,
        help="Number of translation moves per mesh for B3 diagnostics.",
    )
    parser.add_argument(
        "--b3-mesh-sizes",
        default="32,64,128",
        help=(
            "Comma-separated mesh sizes for B3. "
            "Power-of-two values are recommended by current PME implementation."
        ),
    )
    parser.add_argument(
        "--b2-mesh-sizes",
        default="64",
        help="Comma-separated mesh sizes for B2 (default: 64).",
    )
    parser.add_argument(
        "--probe-box-nm",
        type=float,
        default=5.0,
        help="Cubic probe box length in nm for B2/B3 diagnostics.",
    )
    parser.add_argument(
        "--probe-cutoff-nm",
        type=float,
        default=1.2,
        help="Electrostatic cutoff in nm for B2/B3 diagnostics.",
    )
    parser.add_argument(
        "--pme-alpha",
        type=float,
        default=5.0,
        help="PME alpha (nm^-1) for B2/B3 diagnostics.",
    )
    parser.add_argument(
        "--pme-spline-order",
        type=int,
        default=4,
        help="PME B-spline order for B2/B3 diagnostics.",
    )
    parser.add_argument(
        "--pme-tolerance",
        type=float,
        default=1e-6,
        help="PME tolerance for B2/B3 diagnostics.",
    )
    parser.add_argument(
        "--b4-block-size",
        type=int,
        default=20,
        help="Block size in snapshots for B4 block-mean export.",
    )
    parser.add_argument(
        "--report-iat-ess",
        action="store_true",
        help="Compute IAT/ESS diagnostics for B4.",
    )
    parser.add_argument(
        "--timeout-sec",
        type=int,
        default=1200,
        help="Timeout per B1 gcmc_cpu run.",
    )
    parser.add_argument(
        "--disable-plot",
        action="store_true",
        help="Disable matplotlib convergence figure generation.",
    )
    parser.add_argument(
        "--skip-b1",
        action="store_true",
        help="Skip B1/B4 CLI runs (still runs B2/B3 Python diagnostics).",
    )
    parser.add_argument(
        "--quick",
        action="store_true",
        help="Quick smoke mode with reduced steps/seeds.",
    )
    args = parser.parse_args()

    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    b1_seeds = _parse_int_list_csv(args.b1_seeds)
    b3_mesh_sizes = _parse_int_list_csv(args.b3_mesh_sizes)
    b2_mesh_sizes = _parse_int_list_csv(args.b2_mesh_sizes)

    b1_warmup = args.b1_warmup_steps
    b1_production = args.b1_production_steps
    b2_steps = args.b2_steps
    b3_steps = args.b3_steps
    b1_sol_muex = args.b1_sol_muex
    b1_cosolvent_mctime_weight = args.b1_cosolvent_mctime_weight
    b1_sol_mctime_weight = args.b1_sol_mctime_weight
    b1_use_cavity_bias = args.b1_use_cavity_bias
    b1_use_conf_bias = args.b1_use_conf_bias
    b1_num_conf_bias_trial = args.b1_num_conf_bias_trial
    b1_use_insdel_frac = args.b1_use_insdel_frac
    b1_insdel_frac = args.b1_insdel_frac
    b1_use_number_water_nbar = args.b1_use_number_water_nbar
    b1_energy_cutoff_frag = args.b1_energy_cutoff_frag
    b1_energy_cutoff_prot = args.b1_energy_cutoff_prot
    b1_move_probabilities = _parse_float_quad(args.b1_mc_move_prob)

    if args.b1_profile == "old_lysozyme":
        b1_sol_muex = -5.60
        b1_cosolvent_mctime_weight = 0.125
        b1_sol_mctime_weight = 0.125
        b1_use_cavity_bias = "no"
        b1_use_conf_bias = "no"
        b1_use_insdel_frac = "no"
        b1_move_probabilities = (0.5, 0.5, 0.0, 0.0)
        b1_use_number_water_nbar = "no"
        b1_energy_cutoff_frag = None
        b1_energy_cutoff_prot = None
    elif args.b1_profile == "paper_calibrated":
        b1_sol_muex = -2.50
        b1_cosolvent_mctime_weight = 1.0
        b1_sol_mctime_weight = 100.0
        b1_use_cavity_bias = "yes"
        b1_use_conf_bias = "no"
        b1_use_insdel_frac = "no"
        b1_move_probabilities = (0.5, 0.5, 0.0, 0.0)
        b1_use_number_water_nbar = "no"
        b1_energy_cutoff_frag = None
        b1_energy_cutoff_prot = None
    elif args.b1_profile == "old_protein":
        b1_sol_muex = -5.60
        b1_cosolvent_mctime_weight = 0.70
        b1_sol_mctime_weight = 0.30
        b1_use_cavity_bias = "yes"
        b1_use_conf_bias = "yes"
        b1_num_conf_bias_trial = max(1, b1_num_conf_bias_trial)
        b1_use_insdel_frac = "yes"
        b1_insdel_frac = 0.8
        b1_use_number_water_nbar = "yes"
        b1_energy_cutoff_frag = 8.0
        b1_energy_cutoff_prot = 8.0

    if args.quick:
        b1_seeds = b1_seeds[:2]
        b1_warmup = min(b1_warmup, 200)
        b1_production = min(b1_production, 600)
        b2_steps = min(b2_steps, 80)
        b3_steps = min(b3_steps, 120)

    summary: dict = {
        "config": {
            "quick": bool(args.quick),
            "b1_profile": args.b1_profile,
            "b1_seeds": b1_seeds,
            "b1_warmup_steps": b1_warmup,
            "b1_production_steps": b1_production,
            "b1_sol_conc": args.b1_sol_conc,
            "b1_sol_muex": b1_sol_muex,
            "b1_cosolvent_mctime_weight": b1_cosolvent_mctime_weight,
            "b1_sol_mctime_weight": b1_sol_mctime_weight,
            "b1_use_cavity_bias": b1_use_cavity_bias,
            "b1_use_conf_bias": b1_use_conf_bias,
            "b1_num_conf_bias_trial": b1_num_conf_bias_trial,
            "b1_use_insdel_frac": b1_use_insdel_frac,
            "b1_insdel_frac": b1_insdel_frac,
            "b1_mc_move_prob": [float(value) for value in b1_move_probabilities],
            "b1_use_number_water_nbar": b1_use_number_water_nbar,
            "b1_energy_cutoff_frag": b1_energy_cutoff_frag,
            "b1_energy_cutoff_prot": b1_energy_cutoff_prot,
            "b1_max_molecules_per_type": args.b1_max_molecules_per_type,
            "b1_accept_target": args.b1_accept_target,
            "b2_steps": b2_steps,
            "b2_mesh_sizes": b2_mesh_sizes,
            "b3_steps": b3_steps,
            "b3_mesh_sizes": b3_mesh_sizes,
            "probe_box_nm": args.probe_box_nm,
            "probe_cutoff_nm": args.probe_cutoff_nm,
            "pme_alpha": args.pme_alpha,
            "pme_spline_order": args.pme_spline_order,
            "pme_tolerance": args.pme_tolerance,
            "b4_block_size": args.b4_block_size,
            "report_iat_ess": bool(args.report_iat_ess),
        }
    }

    if args.skip_b1:
        summary["b1"] = {"skipped": True}
    else:
        gcmc_cpu = _find_gcmc_cpu(args.gcmc_cpu)
        summary["b1"] = run_b1_and_b4(
            gcmc_cpu=gcmc_cpu,
            out_dir=out_dir,
            seeds=b1_seeds,
            warmup_steps=b1_warmup,
            production_steps=b1_production,
            print_freq=args.b1_print_freq,
            sol_fragconc=args.b1_sol_conc,
            sol_muex=b1_sol_muex,
            cosolvent_mctime_weight=b1_cosolvent_mctime_weight,
            sol_mctime_weight=b1_sol_mctime_weight,
            use_cavity_bias=(b1_use_cavity_bias == "yes"),
            use_conf_bias=(b1_use_conf_bias == "yes"),
            use_insdel_frac=(b1_use_insdel_frac == "yes"),
            insdel_frac=b1_insdel_frac,
            move_probabilities=b1_move_probabilities,
            use_number_water_nbar=(b1_use_number_water_nbar == "yes"),
            energy_cutoff_frag=b1_energy_cutoff_frag,
            energy_cutoff_prot=b1_energy_cutoff_prot,
            num_conf_bias_trial=b1_num_conf_bias_trial,
            max_molecules_per_type=args.b1_max_molecules_per_type,
            accept_target=args.b1_accept_target,
            report_iat_ess=bool(args.report_iat_ess),
            b4_block_size=args.b4_block_size,
            timeout_sec=args.timeout_sec,
            disable_plot=args.disable_plot,
        )

    summary["b2"] = run_b2(
        out_dir=out_dir,
        mesh_sizes=b2_mesh_sizes,
        box_nm=args.probe_box_nm,
        cutoff_nm=args.probe_cutoff_nm,
        alpha_value=args.pme_alpha,
        spline_order=args.pme_spline_order,
        tolerance=args.pme_tolerance,
        steps=b2_steps,
        seed=12345,
    )
    summary["b3"] = run_b3(
        out_dir=out_dir,
        mesh_sizes=b3_mesh_sizes,
        box_nm=args.probe_box_nm,
        cutoff_nm=args.probe_cutoff_nm,
        alpha_value=args.pme_alpha,
        spline_order=args.pme_spline_order,
        tolerance=args.pme_tolerance,
        steps=b3_steps,
        seed=24680,
    )

    summary_path = out_dir / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(f"Validation summary written to: {summary_path}")


if __name__ == "__main__":
    main()

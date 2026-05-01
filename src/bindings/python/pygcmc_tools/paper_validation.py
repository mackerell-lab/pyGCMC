from __future__ import annotations

import argparse
import json
import math
import shutil
import subprocess
from pathlib import Path
from typing import Any

import pygcmc


PGP_PAPER_METHODS = ("pgp_full", "pgp_full_pme")
PME_PAPER_METHOD = "pme"
PAPER_MODES = ("pgp_full", "pme", "pgp_full_pme")
MODE_LABELS = {"pgp_full": "C", "pme": "D", "pgp_full_pme": "E"}
PAIR_ORDER = (("pgp_full_pme", "pgp_full"), ("pme", "pgp_full"), ("pgp_full_pme", "pme"))


def load_jsonl(path: str | Path) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for line in Path(path).read_text().splitlines():
        if line.strip():
            records.append(json.loads(line))
    return records


def _stats(values: list[float]) -> dict[str, float | int]:
    if not values:
        return {"n": 0}
    n = len(values)
    return {
        "n": n,
        "mean": sum(values) / n,
        "mean_abs": sum(abs(value) for value in values) / n,
        "rmse": math.sqrt(sum(value * value for value in values) / n),
        "min": min(values),
        "max": max(values),
        "max_abs": max(abs(value) for value in values),
    }


def _resolve_value_path(base_dir: Path, value: str) -> str:
    candidate = Path(value)
    if candidate.is_absolute() or not value:
        return value
    return str((base_dir / candidate).resolve())


def rewrite_inp_for_mode(
    base_inp: str | Path,
    output_inp: str | Path,
    *,
    mode: str,
    seed: int | None = None,
    steps: int | None = None,
    move_prob: str | None = None,
    use_conf_bias: str | None = None,
    use_cavity_bias: str | None = None,
    num_conf_bias_trial: int | None = None,
    temperature: float | None = None,
) -> Path:
    base_inp = Path(base_inp).resolve()
    output_inp = Path(output_inp).resolve()
    base_dir = base_inp.parent
    path_like_keys = {"par", "fragitp", "atomtypes", "monomerdir", "top", "pdb", "protitp"}

    replacements: dict[str, str] = {"energy_method": f"energy_method:{mode}"}
    if seed is not None:
        replacements["random_seed"] = f"random_seed:{seed}"
    if steps is not None:
        replacements["mcsteps"] = f"mcsteps:{steps}"
        replacements["nprint"] = f"nprint:{steps}"
    if move_prob is not None:
        replacements["mc_move_prob"] = f"mc_move_prob:{move_prob}"
    if use_conf_bias is not None:
        replacements["use_conf_bias"] = f"use_conf_bias:{use_conf_bias}"
    if use_cavity_bias is not None:
        replacements["use_cavity_bias"] = f"use_cavity_bias:{use_cavity_bias}"
    if num_conf_bias_trial is not None:
        replacements["num_conf_bias_trial"] = f"num_conf_bias_trial:{num_conf_bias_trial}"
    if temperature is not None:
        replacements["temperature"] = f"temperature:{temperature:.8g}"

    seen = {key: False for key in replacements}
    output_lines: list[str] = []
    for line in base_inp.read_text().splitlines():
        if ":" not in line:
            output_lines.append(line)
            continue
        key, raw_value = line.split(":", 1)
        key_norm = key.strip().lower()
        if key_norm in replacements:
            output_lines.append(replacements[key_norm])
            seen[key_norm] = True
            continue
        value = raw_value.strip()
        if key_norm in path_like_keys:
            output_lines.append(f"{key.strip()}:{_resolve_value_path(base_dir, value)}")
        else:
            output_lines.append(line)

    for key, replacement in replacements.items():
        if not seen[key]:
            output_lines.append(replacement)

    output_inp.parent.mkdir(parents=True, exist_ok=True)
    output_inp.write_text("\n".join(output_lines) + "\n")
    return output_inp


def _run_mode_pybind(
    inp_path: Path,
    run_dir: Path,
    *,
    seed: int | None,
    diagnostics_buffer: int,
) -> dict[str, Any]:
    config = pygcmc.GCMCCPUConfig()
    config.inputFile = str(inp_path)
    config.outputPrefix = str(run_dir / "gcmc")
    config.randomSeed = -1 if seed is None else seed
    config.printFrequency = 1000000000
    config.trajectoryFrequency = 1000000000
    config.checkpointFrequency = 0
    config.enableStatistics = False

    return dict(
        pygcmc.run_gcmc_cpu(
            config,
            "",
            str(run_dir / "acceptance.jsonl"),
            str(run_dir / "params.json"),
            diagnostics_buffer,
        )
    )


def _run_mode_cli(
    inp_path: Path,
    run_dir: Path,
    *,
    executable: str | None,
    timeout: int,
) -> dict[str, Any]:
    gcmc_cpu = executable or shutil.which("gcmc_cpu")
    if not gcmc_cpu:
        raise FileNotFoundError("gcmc_cpu executable not found in PATH")
    command = [
        gcmc_cpu,
        "--inp",
        str(inp_path),
        "--prefix",
        str(run_dir / "gcmc"),
        "--dump-accept",
        str(run_dir / "acceptance.jsonl"),
        "--dump-params",
        str(run_dir / "params.json"),
        "--print-freq",
        "1000000000",
        "--traj-freq",
        "1000000000",
        "--checkpoint-freq",
        "0",
        "--no-stats",
    ]
    completed = subprocess.run(
        command,
        cwd=str(run_dir),
        capture_output=True,
        text=True,
        timeout=timeout,
        check=False,
    )
    (run_dir / "stdout.log").write_text(completed.stdout)
    (run_dir / "stderr.log").write_text(completed.stderr)
    return {"returncode": completed.returncode, "ran": completed.returncode == 0}


def run_mode(
    base_inp: str | Path,
    out_dir: str | Path,
    *,
    mode: str,
    seed: int | None = None,
    steps: int | None = None,
    move_prob: str | None = None,
    use_conf_bias: str | None = None,
    use_cavity_bias: str | None = None,
    num_conf_bias_trial: int | None = None,
    temperature: float | None = None,
    backend: str = "pybind",
    executable: str | None = None,
    diagnostics_buffer: int = 65536,
    timeout: int = 120,
) -> dict[str, Any]:
    if mode not in PAPER_MODES:
        raise ValueError(f"Unsupported paper mode: {mode}")
    run_dir = Path(out_dir).resolve() / mode
    run_dir.mkdir(parents=True, exist_ok=True)
    inp_path = rewrite_inp_for_mode(
        base_inp,
        run_dir / "run.inp",
        mode=mode,
        seed=seed,
        steps=steps,
        move_prob=move_prob,
        use_conf_bias=use_conf_bias,
        use_cavity_bias=use_cavity_bias,
        num_conf_bias_trial=num_conf_bias_trial,
        temperature=temperature,
    )

    if backend == "pybind":
        result = _run_mode_pybind(inp_path, run_dir, seed=seed, diagnostics_buffer=diagnostics_buffer)
    elif backend == "cli":
        result = _run_mode_cli(inp_path, run_dir, executable=executable, timeout=timeout)
    else:
        raise ValueError("backend must be 'pybind' or 'cli'")

    result["mode"] = mode
    result["label"] = MODE_LABELS[mode]
    result["run_dir"] = str(run_dir)
    result["acceptance_log"] = str(run_dir / "acceptance.jsonl")
    result["params_json"] = str(run_dir / "params.json")
    if (run_dir / "acceptance.jsonl").exists():
        result["record_count"] = len(load_jsonl(run_dir / "acceptance.jsonl"))
    return result


def _first_divergence(rows_a: list[dict[str, Any]], rows_b: list[dict[str, Any]]) -> tuple[int | None, str | None, int]:
    pair_count = min(len(rows_a), len(rows_b))
    for idx in range(pair_count):
        row_a = rows_a[idx]
        row_b = rows_b[idx]
        for key in ("move", "species"):
            if str(row_a.get(key, "")).strip().lower() != str(row_b.get(key, "")).strip().lower():
                return idx + 1, f"{key}_mismatch", idx
        if bool(row_a.get("accepted", False)) != bool(row_b.get("accepted", False)):
            return idx + 1, "accepted_mismatch", idx
        if int(row_a.get("nBefore", -1)) != int(row_b.get("nBefore", -1)):
            return idx + 1, "n_before_mismatch", idx
    return None, None, pair_count


def compare_mode_records(rows_a: list[dict[str, Any]], rows_b: list[dict[str, Any]]) -> dict[str, Any]:
    first_step, first_reason, prefix_limit = _first_divergence(rows_a, rows_b)
    deltas: list[float] = []
    by_move: dict[str, list[float]] = {}
    for idx in range(prefix_limit):
        row_a = rows_a[idx]
        row_b = rows_b[idx]
        delta = float(row_a.get("deltaU", 0.0)) - float(row_b.get("deltaU", 0.0))
        deltas.append(delta)
        move_key = f"{str(row_a.get('move', '')).strip().lower()}:{str(row_a.get('species', '')).strip().lower()}"
        by_move.setdefault(move_key, []).append(delta)
    return {
        "pair_count": min(len(rows_a), len(rows_b)),
        "first_divergence_step": first_step,
        "first_divergence_reason": first_reason,
        "prefix_steps": prefix_limit,
        "all": _stats(deltas),
        "by_move_species": {key: _stats(values) for key, values in sorted(by_move.items())},
    }


def summarize_paper_modes(out_dir: str | Path, modes: tuple[str, ...] = PAPER_MODES) -> dict[str, Any]:
    out_dir = Path(out_dir).resolve()
    records = {mode: load_jsonl(out_dir / mode / "acceptance.jsonl") for mode in modes}
    pairs: dict[str, Any] = {}
    for mode_a, mode_b in PAIR_ORDER:
        if mode_a not in records or mode_b not in records:
            continue
        label = f"{MODE_LABELS[mode_a]}-{MODE_LABELS[mode_b]}"
        pairs[label] = compare_mode_records(records[mode_a], records[mode_b])
    return {
        "modes": list(modes),
        "labels": {mode: MODE_LABELS[mode] for mode in modes},
        "record_counts": {mode: len(rows) for mode, rows in records.items()},
        "pairs": pairs,
    }


def run_paper_modes(
    base_inp: str | Path,
    out_dir: str | Path,
    *,
    modes: tuple[str, ...] = PAPER_MODES,
    seed: int | None = None,
    steps: int | None = None,
    move_prob: str | None = None,
    use_conf_bias: str | None = None,
    use_cavity_bias: str | None = None,
    num_conf_bias_trial: int | None = None,
    temperature: float | None = None,
    backend: str = "pybind",
    executable: str | None = None,
    diagnostics_buffer: int = 65536,
    timeout: int = 120,
) -> dict[str, Any]:
    out_dir = Path(out_dir).resolve()
    runs = {
        mode: run_mode(
            base_inp,
            out_dir,
            mode=mode,
            seed=seed,
            steps=steps,
            move_prob=move_prob,
            use_conf_bias=use_conf_bias,
            use_cavity_bias=use_cavity_bias,
            num_conf_bias_trial=num_conf_bias_trial,
            temperature=temperature,
            backend=backend,
            executable=executable,
            diagnostics_buffer=diagnostics_buffer,
            timeout=timeout,
        )
        for mode in modes
    }
    summary = summarize_paper_modes(out_dir, modes=modes)
    summary["runs"] = runs
    summary_path = out_dir / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    summary["summary_path"] = str(summary_path)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description="Run installed PGP/PME paper validation modes.")
    parser.add_argument("--base-inp", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--backend", choices=("pybind", "cli"), default="pybind")
    parser.add_argument("--executable", default=None)
    parser.add_argument("--seed", type=int, default=None)
    parser.add_argument("--steps", type=int, default=None)
    parser.add_argument("--move-prob", default=None)
    parser.add_argument("--use-conf-bias", default=None)
    parser.add_argument("--use-cavity-bias", default=None)
    parser.add_argument("--num-conf-bias-trial", type=int, default=None)
    parser.add_argument("--temperature", type=float, default=None)
    parser.add_argument("--modes", default=",".join(PAPER_MODES))
    args = parser.parse_args()

    modes = tuple(mode.strip() for mode in args.modes.split(",") if mode.strip())
    summary = run_paper_modes(
        args.base_inp,
        args.out_dir,
        modes=modes,
        seed=args.seed,
        steps=args.steps,
        move_prob=args.move_prob,
        use_conf_bias=args.use_conf_bias,
        use_cavity_bias=args.use_cavity_bias,
        num_conf_bias_trial=args.num_conf_bias_trial,
        temperature=args.temperature,
        backend=args.backend,
        executable=args.executable,
    )
    print(summary["summary_path"])


if __name__ == "__main__":
    main()

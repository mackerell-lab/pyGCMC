"""
Paired C/E deltaU audit runner for fixed-input gcmc_cpu comparisons.

Goal:
- run Mode C (`pgp_full`) and Mode E (`pgp_full_pme`) with identical seeds/proposals
- pair accept-log records by index (and validate move/species identity)
- quantify deltaU(E-C) statistics by move/species buckets

This is a reproducible diagnostic script (not a pytest test).
"""

from __future__ import annotations

import argparse
import json
import subprocess
from collections import defaultdict
from pathlib import Path


def _parse_int_list_csv(values: str) -> list[int]:
    parsed: list[int] = []
    for token in values.split(","):
        stripped = token.strip()
        if stripped:
            parsed.append(int(stripped))
    if not parsed:
        raise ValueError("At least one seed is required.")
    return parsed


def _find_gcmc_cpu(explicit_path: str | None) -> Path:
    if explicit_path:
        candidate = Path(explicit_path).resolve()
        if candidate.exists():
            return candidate
        raise FileNotFoundError(f"gcmc_cpu not found: {candidate}")
    fallback = Path("pygcmc_dev/build/bin/gcmc_cpu").resolve()
    if fallback.exists():
        return fallback
    result = subprocess.run(["which", "gcmc_cpu"], capture_output=True, text=True, check=False)
    if result.returncode == 0 and result.stdout.strip():
        return Path(result.stdout.strip())
    raise FileNotFoundError("gcmc_cpu executable not found. Build target gcmc_cpu first.")


def _load_records(path: Path) -> list[dict]:
    records: list[dict] = []
    for line in path.read_text().splitlines():
        if line.strip():
            records.append(json.loads(line))
    return records


def _stats(values: list[float]) -> dict:
    if not values:
        return {"n": 0}
    abs_values = sorted(abs(v) for v in values)

    def _percentile_abs(p: float) -> float:
        if len(abs_values) == 1:
            return abs_values[0]
        index = int(round((len(abs_values) - 1) * p))
        index = max(0, min(len(abs_values) - 1, index))
        return abs_values[index]

    return {
        "n": len(values),
        "mean": sum(values) / len(values),
        "mean_abs": sum(abs(v) for v in values) / len(values),
        "p95_abs": _percentile_abs(0.95),
        "p99_abs": _percentile_abs(0.99),
        "max_abs": abs_values[-1],
        "min": min(values),
        "max": max(values),
    }


def _rewrite_inp(
    *,
    base_inp: Path,
    output_inp: Path,
    mode: str,
    seed: int,
    steps: int,
    use_cavity_bias: str,
    temperature: float,
) -> None:
    lines_out: list[str] = []
    for line in base_inp.read_text().splitlines():
        if line.startswith("energy_method:"):
            lines_out.append(f"energy_method:{mode}")
        elif line.startswith("random_seed:"):
            lines_out.append(f"random_seed:{seed}")
        elif line.startswith("mcsteps:"):
            lines_out.append(f"mcsteps:{steps}")
        elif line.startswith("nprint:"):
            lines_out.append(f"nprint:{steps}")
        elif line.startswith("mc_move_prob:"):
            lines_out.append("mc_move_prob:0.500000 0.500000 0.000000 0.000000")
        elif line.startswith("use_cavity_bias:"):
            lines_out.append(f"use_cavity_bias:{use_cavity_bias}")
        elif line.startswith("temperature:"):
            lines_out.append(f"temperature:{temperature:.8g}")
        else:
            lines_out.append(line)
    output_inp.parent.mkdir(parents=True, exist_ok=True)
    output_inp.write_text("\n".join(lines_out) + "\n")


def _run_mode(
    *,
    gcmc_cpu: Path,
    inp_path: Path,
    run_dir: Path,
) -> None:
    run_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        str(gcmc_cpu),
        "--inp",
        str(inp_path),
        "--prefix",
        str(run_dir / "gcmc"),
        "--dump-accept",
        str(run_dir / "acceptance.jsonl"),
        "--dump-params",
        str(run_dir / "params.json"),
        "--print-freq",
        "1000000",
        "--traj-freq",
        "1000000",
        "--checkpoint-freq",
        "0",
        "--max-molecules-per-type",
        "-1",
    ]
    result = subprocess.run(
        cmd,
        cwd=str(Path.cwd()),
        capture_output=True,
        text=True,
        check=False,
    )
    (run_dir / "stdout.log").write_text(result.stdout)
    (run_dir / "stderr.log").write_text(result.stderr)
    if result.returncode != 0:
        raise RuntimeError(
            f"gcmc_cpu failed for {inp_path}\n"
            f"stdout:\n{result.stdout}\n"
            f"stderr:\n{result.stderr}"
        )


def _summarize_seed(seed_dir: Path) -> dict:
    rows_c = _load_records(seed_dir / "pgp_full" / "acceptance.jsonl")
    rows_e = _load_records(seed_dir / "pgp_full_pme" / "acceptance.jsonl")
    pair_count = min(len(rows_c), len(rows_e))
    pairs: list[dict] = []
    for index in range(pair_count):
        rec_c = rows_c[index]
        rec_e = rows_e[index]
        move_c = str(rec_c.get("move", "")).strip().lower()
        move_e = str(rec_e.get("move", "")).strip().lower()
        species_c = str(rec_c.get("species", "")).strip().lower()
        species_e = str(rec_e.get("species", "")).strip().lower()
        pairs.append(
            {
                "move_c": move_c,
                "move_e": move_e,
                "species_c": species_c,
                "species_e": species_e,
                "delta_e_minus_c": float(rec_e.get("deltaU", 0.0)) - float(rec_c.get("deltaU", 0.0)),
            }
        )

    same_move_count = sum(1 for row in pairs if row["move_c"] == row["move_e"])
    same_species_count = sum(1 for row in pairs if row["species_c"] == row["species_e"])
    matched = [row for row in pairs if row["move_c"] == row["move_e"] and row["species_c"] == row["species_e"]]

    by_move_species: dict[str, list[float]] = defaultdict(list)
    all_values: list[float] = []
    for row in matched:
        key = f"{row['move_c']}:{row['species_c']}"
        value = row["delta_e_minus_c"]
        by_move_species[key].append(value)
        all_values.append(value)

    move_species_stats = {key: _stats(values) for key, values in sorted(by_move_species.items())}
    sol_insertion = by_move_species.get("insertion:sol", [])
    sol_deletion = by_move_species.get("deletion:sol", [])

    first_state_divergence_step: int | None = None
    first_state_divergence_reason: str | None = None
    for index in range(pair_count):
        rec_c = rows_c[index]
        rec_e = rows_e[index]
        if str(rec_c.get("move", "")).strip().lower() != str(rec_e.get("move", "")).strip().lower():
            first_state_divergence_step = index + 1
            first_state_divergence_reason = "move_mismatch"
            break
        if str(rec_c.get("species", "")).strip().lower() != str(rec_e.get("species", "")).strip().lower():
            first_state_divergence_step = index + 1
            first_state_divergence_reason = "species_mismatch"
            break
        if bool(rec_c.get("accepted", False)) != bool(rec_e.get("accepted", False)):
            first_state_divergence_step = index + 1
            first_state_divergence_reason = "accepted_mismatch"
            break
        if int(rec_c.get("nBefore", -1)) != int(rec_e.get("nBefore", -1)):
            first_state_divergence_step = index + 1
            first_state_divergence_reason = "n_before_mismatch"
            break

    prefix_limit = pair_count if first_state_divergence_step is None else max(0, first_state_divergence_step - 1)
    prefix_values: list[float] = []
    prefix_by_move_species: dict[str, list[float]] = defaultdict(list)
    for index in range(prefix_limit):
        rec_c = rows_c[index]
        rec_e = rows_e[index]
        move_c = str(rec_c.get("move", "")).strip().lower()
        move_e = str(rec_e.get("move", "")).strip().lower()
        species_c = str(rec_c.get("species", "")).strip().lower()
        species_e = str(rec_e.get("species", "")).strip().lower()
        if move_c != move_e or species_c != species_e:
            continue
        delta = float(rec_e.get("deltaU", 0.0)) - float(rec_c.get("deltaU", 0.0))
        prefix_values.append(delta)
        prefix_by_move_species[f"{move_c}:{species_c}"].append(delta)

    return {
        "counts": {
            "pgp_full": len(rows_c),
            "pgp_full_pme": len(rows_e),
            "paired_by_index": pair_count,
            "same_move_count": same_move_count,
            "same_species_count": same_species_count,
            "same_move_ratio": same_move_count / pair_count if pair_count else 0.0,
            "same_species_ratio": same_species_count / pair_count if pair_count else 0.0,
            "matched_count": len(matched),
            "matched_ratio": len(matched) / pair_count if pair_count else 0.0,
        },
        "pairing_window": {
            "first_state_divergence_step": first_state_divergence_step,
            "first_state_divergence_reason": first_state_divergence_reason,
            "prefix_paired_steps": prefix_limit,
            "prefix_ratio": prefix_limit / pair_count if pair_count else 0.0,
        },
        "stats": {
            "all_matched": _stats(all_values),
            "sol_insertion": _stats(sol_insertion),
            "sol_deletion": _stats(sol_deletion),
            "by_move_species": move_species_stats,
            "prefix_all_paired": _stats(prefix_values),
            "prefix_sol_insertion": _stats(prefix_by_move_species.get("insertion:sol", [])),
            "prefix_sol_deletion": _stats(prefix_by_move_species.get("deletion:sol", [])),
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Run paired C/E deltaU audit for gcmc_cpu.")
    parser.add_argument(
        "--base-inp",
        default="tmp/pgp_mu_recenter_bracket_m14p0_seed101_ref/run.inp",
        help="Base INP path to clone and rewrite.",
    )
    parser.add_argument(
        "--out-dir",
        default="tmp/pgp_ce_deltau_audit_runner",
        help="Output directory.",
    )
    parser.add_argument(
        "--seeds",
        default="246810",
        help="Comma-separated seed list (e.g., 101,202,303).",
    )
    parser.add_argument("--steps", type=int, default=800, help="MC steps per run.")
    parser.add_argument(
        "--use-cavity-bias",
        choices=["yes", "no"],
        default="no",
        help="Override use_cavity_bias in generated INP.",
    )
    parser.add_argument("--temperature", type=float, default=300.0, help="Override temperature.")
    parser.add_argument("--gcmc-cpu", default=None, help="Path to gcmc_cpu executable.")
    args = parser.parse_args()

    base_inp = Path(args.base_inp).resolve()
    if not base_inp.exists():
        raise FileNotFoundError(f"Base INP not found: {base_inp}")
    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    gcmc_cpu = _find_gcmc_cpu(args.gcmc_cpu)
    seeds = _parse_int_list_csv(args.seeds)

    summary = {
        "config": {
            "base_inp": str(base_inp),
            "steps": args.steps,
            "use_cavity_bias": args.use_cavity_bias,
            "temperature": args.temperature,
            "seeds": seeds,
            "gcmc_cpu": str(gcmc_cpu),
            "mc_move_prob_override": [0.5, 0.5, 0.0, 0.0],
        },
        "seeds": {},
    }

    for seed in seeds:
        seed_dir = out_dir / f"seed_{seed}"
        for mode in ["pgp_full", "pgp_full_pme"]:
            mode_dir = seed_dir / mode
            inp_path = mode_dir / "run.inp"
            _rewrite_inp(
                base_inp=base_inp,
                output_inp=inp_path,
                mode=mode,
                seed=seed,
                steps=args.steps,
                use_cavity_bias=args.use_cavity_bias,
                temperature=args.temperature,
            )
            _run_mode(gcmc_cpu=gcmc_cpu, inp_path=inp_path, run_dir=mode_dir)
        summary["seeds"][str(seed)] = _summarize_seed(seed_dir)

    out_path = out_dir / "summary.json"
    out_path.write_text(json.dumps(summary, indent=2) + "\n")
    print(f"Paired deltaU audit summary written to: {out_path}")


if __name__ == "__main__":
    main()

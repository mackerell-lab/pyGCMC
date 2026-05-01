#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any


MODE_LABELS = {"pgp_full": "C", "pme": "D", "pgp_full_pme": "E"}
MODE_DESCRIPTIONS = {
    "pgp_full": "PGP cross-only; PME mesh-self excluded",
    "pme": "standard PME reference",
    "pgp_full_pme": "PGP with mesh-self restored; PME-equivalent",
}
PAIRS = (("pgp_full_pme", "pme"), ("pgp_full", "pme"), ("pgp_full_pme", "pgp_full"))


def load_jsonl(path: Path) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    if not path.exists():
        return records
    for line in path.read_text().splitlines():
        if line.strip():
            records.append(json.loads(line))
    return records


def percentile_abs(values: list[float], percentile: float) -> float:
    if not values:
        return float("nan")
    ordered = sorted(abs(value) for value in values)
    index = min(len(ordered) - 1, max(0, math.ceil(percentile / 100.0 * len(ordered)) - 1))
    return ordered[index]


def stats(values: list[float]) -> dict[str, float | int]:
    if not values:
        return {"n": 0}
    return {
        "n": len(values),
        "mean": sum(values) / len(values),
        "mean_abs": sum(abs(value) for value in values) / len(values),
        "rmse": math.sqrt(sum(value * value for value in values) / len(values)),
        "p95_abs": percentile_abs(values, 95.0),
        "max_abs": max(abs(value) for value in values),
        "min": min(values),
        "max": max(values),
    }


def common_prefix(rows_a: list[dict[str, Any]], rows_b: list[dict[str, Any]]) -> tuple[int, str | None]:
    limit = min(len(rows_a), len(rows_b))
    for index in range(limit):
        row_a = rows_a[index]
        row_b = rows_b[index]
        for key in ("move", "species", "nBefore", "accepted"):
            if row_a.get(key) != row_b.get(key):
                return index, f"{key} mismatch at paired row {index + 1}"
    return limit, None


def compare(rows_a: list[dict[str, Any]], rows_b: list[dict[str, Any]]) -> dict[str, Any]:
    limit, divergence = common_prefix(rows_a, rows_b)
    deltas = [float(rows_a[i].get("deltaU", 0.0)) - float(rows_b[i].get("deltaU", 0.0)) for i in range(limit)]
    by_move: dict[str, list[float]] = {}
    for i, delta in enumerate(deltas):
        key = f"{rows_a[i].get('move', 'unknown')}:{rows_a[i].get('species', 'unknown')}"
        by_move.setdefault(key, []).append(delta)
    return {
        "paired_rows": limit,
        "first_divergence": divergence,
        "all": stats(deltas),
        "by_move_species": {key: stats(values) for key, values in sorted(by_move.items())},
    }


def load_wall_times(runs_dir: Path) -> dict[str, float]:
    path = runs_dir / "wall_time.tsv"
    times: dict[str, float] = {}
    if not path.exists():
        return times
    for line in path.read_text().splitlines():
        parts = line.split()
        if len(parts) == 2:
            times[parts[0]] = float(parts[1])
    return times


def load_reported_performance(runs_dir: Path) -> dict[str, float]:
    values: dict[str, float] = {}
    for mode in MODE_LABELS:
        path = runs_dir / mode / "stdout.log"
        if not path.exists():
            continue
        last_value: float | None = None
        for line in path.read_text(errors="replace").splitlines():
            if "Performance:" not in line or "steps/s" not in line:
                continue
            try:
                last_value = float(line.split("Performance:", 1)[1].split("steps/s", 1)[0].strip())
            except ValueError:
                continue
        if last_value is not None:
            values[mode] = last_value
    return values


def summarize(runs_dir: Path, system_name: str) -> dict[str, Any]:
    records = {
        mode: load_jsonl(runs_dir / mode / "acceptance.jsonl")
        for mode in MODE_LABELS
        if (runs_dir / mode).exists()
    }
    pairs: dict[str, Any] = {}
    for mode_a, mode_b in PAIRS:
        if mode_a in records and mode_b in records:
            pairs[f"{MODE_LABELS[mode_a]}-{MODE_LABELS[mode_b]}"] = compare(records[mode_a], records[mode_b])
    return {
        "system": system_name,
        "mode_labels": MODE_LABELS,
        "mode_descriptions": MODE_DESCRIPTIONS,
        "record_counts": {mode: len(rows) for mode, rows in records.items()},
        "wall_times_s": load_wall_times(runs_dir),
        "reported_steps_per_s": load_reported_performance(runs_dir),
        "pairs": pairs,
    }


def format_float(value: Any) -> str:
    if isinstance(value, int):
        return str(value)
    if not isinstance(value, float) or math.isnan(value):
        return "NA"
    if value == 0:
        return "0"
    if abs(value) < 1.0e-4 or abs(value) >= 1.0e5:
        return f"{value:.6e}"
    return f"{value:.6f}"


def render_log(summary: dict[str, Any]) -> str:
    lines: list[str] = []
    lines.append(f"PGP/PME paper example summary: {summary['system']}")
    lines.append("")
    lines.append("Mode labels:")
    for mode, label in summary["mode_labels"].items():
        if mode in summary["record_counts"]:
            lines.append(f"  {label} = {mode}: {summary['mode_descriptions'][mode]}")
    lines.append("")
    lines.append("Record counts:")
    for mode, count in summary["record_counts"].items():
        lines.append(f"  {summary['mode_labels'][mode]} ({mode}): {count}")
    if summary["wall_times_s"]:
        lines.append("")
        lines.append("Observed wrapper wall timing:")
        for mode, seconds in summary["wall_times_s"].items():
            count = summary["record_counts"].get(mode, 0)
            steps_per_s = count / seconds if seconds > 0 else float("nan")
            lines.append(
                f"  {summary['mode_labels'].get(mode, '?')} ({mode}): "
                f"{seconds:.3f} s, {steps_per_s:.6f} accepted-log records/s"
            )
    if summary["reported_steps_per_s"]:
        lines.append("")
        lines.append("gcmc_cpu reported performance:")
        for mode, value in summary["reported_steps_per_s"].items():
            lines.append(f"  {summary['mode_labels'].get(mode, '?')} ({mode}): {value:.6f} steps/s")
        for numerator in ("pgp_full", "pgp_full_pme"):
            denominator = "pme"
            if numerator in summary["reported_steps_per_s"] and denominator in summary["reported_steps_per_s"]:
                pme_value = summary["reported_steps_per_s"][denominator]
                if pme_value > 0:
                    ratio = summary["reported_steps_per_s"][numerator] / pme_value
                    lines.append(
                        f"  speedup {summary['mode_labels'][numerator]}/{summary['mode_labels'][denominator]}: "
                        f"{ratio:.1f}x"
                    )
    lines.append("")
    lines.append("Paired trial-energy differences, ΔΔU = ΔU(first mode) - ΔU(second mode):")
    for pair, data in summary["pairs"].items():
        all_stats = data["all"]
        lines.append(
            f"  {pair}: n={all_stats['n']}, mean={format_float(all_stats.get('mean'))}, "
            f"RMSE={format_float(all_stats.get('rmse'))}, "
            f"P95|.|={format_float(all_stats.get('p95_abs'))}, "
            f"Max|.|={format_float(all_stats.get('max_abs'))} kJ/mol"
        )
        if data["first_divergence"]:
            lines.append(f"    paired-prefix warning: {data['first_divergence']}")
        for move_key, move_stats in data["by_move_species"].items():
            lines.append(
                f"    {move_key}: n={move_stats['n']}, "
                f"mean={format_float(move_stats.get('mean'))}, "
                f"Max|.|={format_float(move_stats.get('max_abs'))} kJ/mol"
            )
    lines.append("")
    lines.append("Paper interpretation:")
    lines.append("  C is the cross-only PGP surface used to remove position-dependent PME mesh-self artifacts.")
    lines.append("  E-C reports the mesh-self-restoration term inside the PGP implementation for the same trial.")
    lines.append("  C-D and E-D compare the PGP paths with the slow raw PME reference path for the identical proposal.")
    lines.append("  PGP timing should be much faster than PME on these manuscript-scale systems.")
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser(description="Summarize PGP/PME example acceptance logs.")
    parser.add_argument("runs_dir", type=Path)
    parser.add_argument("--system", default="unknown")
    parser.add_argument("--json", type=Path, default=None)
    parser.add_argument("--log", type=Path, default=None)
    args = parser.parse_args()

    summary = summarize(args.runs_dir.resolve(), args.system)
    text = render_log(summary)
    if args.json:
        args.json.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    if args.log:
        args.log.write_text(text)
    print(text, end="")


if __name__ == "__main__":
    main()

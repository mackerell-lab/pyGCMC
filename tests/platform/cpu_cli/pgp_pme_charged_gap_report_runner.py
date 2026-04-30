"""
Charged D/E gap report runner for PGP paper tracking.

This utility compares Mode D (`pme`) and Mode E (`pgp_full_pme`) under
charged insertion/deletion scenarios, then writes per-seed and aggregate
difference statistics for manuscript/SI reporting.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
import subprocess
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path


MODES = ("pme", "pgp_full_pme")


@dataclass(frozen=True)
class ScenarioConfig:
    name: str
    move: str
    mcsteps: int
    initial_frg_count: int
    fragconc: float
    temperature: float
    max_translation: float = 0.0
    max_rotation: float = 0.0


SCENARIOS = (
    ScenarioConfig(
        name="charged_insertion",
        move="insertion",
        mcsteps=20,
        initial_frg_count=0,
        fragconc=1.0,
        temperature=1.0e30,
    ),
    ScenarioConfig(
        name="charged_deletion",
        move="deletion",
        mcsteps=20,
        initial_frg_count=20,
        fragconc=1.0e-4,
        temperature=1.0e30,
    ),
)


def _find_gcmc_cpu(explicit_path: str | None) -> Path:
    if explicit_path:
        candidate = Path(explicit_path).resolve()
        if candidate.exists():
            return candidate
        raise FileNotFoundError(f"gcmc_cpu not found: {candidate}")

    default = Path("pygcmc_dev/build/bin/gcmc_cpu").resolve()
    if default.exists():
        return default
    raise FileNotFoundError("gcmc_cpu executable not found (expected at pygcmc_dev/build/bin/gcmc_cpu).")


def _write_text(path: Path, content: str) -> None:
    path.write_text(content.strip() + "\n")


def _write_pdb(path: Path, scenario: ScenarioConfig) -> None:
    lines = [
        "CRYST1   40.000   40.000   40.000  90.00  90.00  90.00 P 1           1",
        "ATOM      1  QN  HST A   1      10.000  20.000  20.000  1.00  0.00           C",
    ]
    for index in range(scenario.initial_frg_count):
        atom_id = index + 2
        resid = index + 1
        x = 25.0 + 0.4 * index
        y = 20.0 + 0.2 * (index % 5)
        z = 20.0 + 0.2 * (index % 7)
        lines.append(
            f"ATOM  {atom_id:5d}  QP  FRG B{resid:4d}      "
            f"{x:6.3f}  {y:6.3f}  {z:6.3f}  1.00  0.00           C"
        )
    lines.append("END")
    _write_text(path, "\n".join(lines))


def _write_top(path: Path, scenario: ScenarioConfig) -> None:
    _write_text(
        path,
        f"""
[ defaults ]
1 2 yes 1.0 1.0

[ moleculetype ]
HST  1

[ atoms ]
; nr  type  resnr  residue  atom  cgnr  charge  mass
1   QN    1      HST      QN    1    -1.000   12.011

[ moleculetype ]
FRG  1

[ atoms ]
; nr  type  resnr  residue  atom  cgnr  charge  mass
1   QP    1      FRG      QP    1     1.000   12.011

[ system ]
ChargedGapAudit

[ molecules ]
HST  1
FRG  {scenario.initial_frg_count}
""",
    )


def _write_par(path: Path) -> None:
    _write_text(
        path,
        """
[ defaults ]
1 2 yes 1.0 1.0

[ atomtypes ]
; name  at.num  mass   charge  ptype  sigma   epsilon
QP      0       12.011 0.000   A      0.300   0.000
QN      0       12.011 0.000   A      0.300   0.000
""",
    )


def _write_frag(path: Path) -> None:
    _write_text(
        path,
        """
[ moleculetype ]
FRG  1

[ atoms ]
; nr  type  resnr  residue  atom  cgnr  charge  mass
1   QP    1      FRG      QP    1     1.000  12.011
""",
    )


def _move_prob_vector(move: str) -> str:
    if move == "insertion":
        return "1 0 0 0"
    if move == "deletion":
        return "0 1 0 0"
    raise ValueError(f"Unsupported move type: {move}")


def _write_inp(
    path: Path,
    *,
    scenario: ScenarioConfig,
    method: str,
    seed: int,
    pdb: Path,
    top: Path,
    par: Path,
    frag: Path,
) -> None:
    _write_text(
        path,
        f"""
random_seed:{seed}
energy_method:{method}
par:{par}
fragitp:{frag}
fragname:FRG
fragconc:{scenario.fragconc}
fragmuex:0.0

pdb:{pdb}
top:{top}
box_size:40.0 40.0 40.0
cutoff:6.0
temperature:{scenario.temperature:.8g}
moves_per_step:1
mcsteps:{scenario.mcsteps}
nprint:{scenario.mcsteps}
mc_move_prob:{_move_prob_vector(scenario.move)}
max_translation:{scenario.max_translation}
max_rotation:{scenario.max_rotation}

use_cavity_bias:no
use_conf_bias:no
""",
    )


def _run_gcmc_cpu(
    *,
    gcmc_cpu: Path,
    inp: Path,
    out_prefix: Path,
    accept_log: Path,
    timeout_s: int,
) -> subprocess.CompletedProcess[str]:
    cmd = [
        str(gcmc_cpu),
        "--inp",
        str(inp),
        "--prefix",
        str(out_prefix),
        "--dump-accept",
        str(accept_log),
        "--no-stats",
        "--print-freq",
        "1000000000",
        "--traj-freq",
        "1000000000",
        "--checkpoint-freq",
        "0",
    ]
    return subprocess.run(cmd, capture_output=True, text=True, timeout=timeout_s, check=False)


def _load_records(path: Path, *, move: str) -> tuple[list[dict], int]:
    rows: list[dict] = []
    invalid_lines = 0
    for raw in path.read_text().splitlines():
        if not raw.strip():
            continue
        try:
            rec = json.loads(raw)
        except json.JSONDecodeError:
            invalid_lines += 1
            continue
        if str(rec.get("move", "")).strip().lower() != move:
            continue
        if str(rec.get("species", "")).strip().upper() != "FRG":
            continue
        rows.append(rec)
    return rows, invalid_lines


def _quantile_abs(values: list[float], q: float) -> float:
    if not values:
        return float("nan")
    ordered = sorted(abs(v) for v in values)
    if len(ordered) == 1:
        return ordered[0]
    pos = q * (len(ordered) - 1)
    low = int(math.floor(pos))
    high = int(math.ceil(pos))
    if low == high:
        return ordered[low]
    frac = pos - low
    return ordered[low] * (1.0 - frac) + ordered[high] * frac


def _stats(values: list[float]) -> dict:
    if not values:
        return {
            "n": 0,
            "mean": None,
            "std": None,
            "min": None,
            "max": None,
            "mean_abs": None,
            "p50_abs": None,
            "p95_abs": None,
            "max_abs": None,
        }
    return {
        "n": len(values),
        "mean": statistics.fmean(values),
        "std": statistics.pstdev(values) if len(values) > 1 else 0.0,
        "min": min(values),
        "max": max(values),
        "mean_abs": statistics.fmean(abs(v) for v in values),
        "p50_abs": _quantile_abs(values, 0.50),
        "p95_abs": _quantile_abs(values, 0.95),
        "max_abs": max(abs(v) for v in values),
    }


def _run_scenario_seed(
    *,
    gcmc_cpu: Path,
    out_dir: Path,
    scenario: ScenarioConfig,
    seed: int,
    timeout_s: int,
) -> dict:
    seed_dir = out_dir / scenario.name / f"seed_{seed}"
    seed_dir.mkdir(parents=True, exist_ok=True)

    pdb = seed_dir / "sys.pdb"
    top = seed_dir / "sys.top"
    par = seed_dir / "par.itp"
    frag = seed_dir / "frag.itp"
    _write_pdb(pdb, scenario)
    _write_top(top, scenario)
    _write_par(par)
    _write_frag(frag)

    run_outputs: dict[str, dict] = {}
    for mode in MODES:
        mode_dir = seed_dir / mode
        mode_dir.mkdir(parents=True, exist_ok=True)
        inp = mode_dir / "run.inp"
        out_prefix = mode_dir / "out" / "gcmc"
        out_prefix.parent.mkdir(parents=True, exist_ok=True)
        accept_log = mode_dir / "out" / "accept.jsonl"
        _write_inp(
            inp,
            scenario=scenario,
            method=mode,
            seed=seed,
            pdb=pdb,
            top=top,
            par=par,
            frag=frag,
        )
        result = _run_gcmc_cpu(
            gcmc_cpu=gcmc_cpu,
            inp=inp,
            out_prefix=out_prefix,
            accept_log=accept_log,
            timeout_s=timeout_s,
        )
        (mode_dir / "stdout.log").write_text(result.stdout)
        (mode_dir / "stderr.log").write_text(result.stderr)
        if result.returncode != 0:
            raise RuntimeError(
                f"gcmc_cpu failed for {scenario.name}/seed_{seed}/{mode}\n"
                f"stdout:\n{result.stdout}\n\nstderr:\n{result.stderr}"
            )

        rows, invalid = _load_records(accept_log, move=scenario.move)
        run_outputs[mode] = {"rows": rows, "invalid_lines": invalid}

    rows_d = run_outputs["pme"]["rows"]
    rows_e = run_outputs["pgp_full_pme"]["rows"]
    pair_count = min(len(rows_d), len(rows_e))

    pair_rows = []
    diffs = []
    accepted_mismatch = 0
    for i in range(pair_count):
        row_d = rows_d[i]
        row_e = rows_e[i]
        delta_d = float(row_d["deltaU"])
        delta_e = float(row_e["deltaU"])
        diff_de = delta_d - delta_e
        diffs.append(diff_de)
        acc_d = bool(row_d.get("accepted", False))
        acc_e = bool(row_e.get("accepted", False))
        if acc_d != acc_e:
            accepted_mismatch += 1
        pair_rows.append(
            {
                "scenario": scenario.name,
                "seed": seed,
                "pair_index": i + 1,
                "delta_pme": delta_d,
                "delta_mode_e": delta_e,
                "delta_diff_d_minus_e": diff_de,
                "accepted_pme": acc_d,
                "accepted_mode_e": acc_e,
            }
        )

    return {
        "scenario": scenario.name,
        "seed": seed,
        "move": scenario.move,
        "records_pme": len(rows_d),
        "records_mode_e": len(rows_e),
        "invalid_json_lines_pme": run_outputs["pme"]["invalid_lines"],
        "invalid_json_lines_mode_e": run_outputs["pgp_full_pme"]["invalid_lines"],
        "pair_count": pair_count,
        "accepted_mismatch_count": accepted_mismatch,
        "delta_diff_stats": _stats(diffs),
        "pair_rows": pair_rows,
    }


def _write_seed_csv(path: Path, rows: list[dict]) -> None:
    fields = [
        "scenario",
        "seed",
        "move",
        "records_pme",
        "records_mode_e",
        "invalid_json_lines_pme",
        "invalid_json_lines_mode_e",
        "pair_count",
        "accepted_mismatch_count",
        "n",
        "mean",
        "std",
        "min",
        "max",
        "mean_abs",
        "p50_abs",
        "p95_abs",
        "max_abs",
    ]
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            stats = row["delta_diff_stats"]
            writer.writerow(
                {
                    "scenario": row["scenario"],
                    "seed": row["seed"],
                    "move": row["move"],
                    "records_pme": row["records_pme"],
                    "records_mode_e": row["records_mode_e"],
                    "invalid_json_lines_pme": row["invalid_json_lines_pme"],
                    "invalid_json_lines_mode_e": row["invalid_json_lines_mode_e"],
                    "pair_count": row["pair_count"],
                    "accepted_mismatch_count": row["accepted_mismatch_count"],
                    "n": stats["n"],
                    "mean": stats["mean"],
                    "std": stats["std"],
                    "min": stats["min"],
                    "max": stats["max"],
                    "mean_abs": stats["mean_abs"],
                    "p50_abs": stats["p50_abs"],
                    "p95_abs": stats["p95_abs"],
                    "max_abs": stats["max_abs"],
                }
            )


def _write_pair_csv(path: Path, rows: list[dict]) -> None:
    fields = [
        "scenario",
        "seed",
        "pair_index",
        "delta_pme",
        "delta_mode_e",
        "delta_diff_d_minus_e",
        "accepted_pme",
        "accepted_mode_e",
    ]
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate charged D/E gap report by seed.")
    parser.add_argument("--gcmc-cpu", default=None, help="Path to gcmc_cpu binary.")
    parser.add_argument(
        "--out-dir",
        default=None,
        help="Output directory (default: tmp/pgp_pme_charged_gap_report_<timestamp>).",
    )
    parser.add_argument("--seeds", default="12345,246810,97531,135791,864209")
    parser.add_argument("--timeout", type=int, default=120)
    args = parser.parse_args()

    gcmc_cpu = _find_gcmc_cpu(args.gcmc_cpu)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_dir = Path(args.out_dir) if args.out_dir else Path("tmp") / f"pgp_pme_charged_gap_report_{timestamp}"
    out_dir = out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    seeds = [int(item.strip()) for item in args.seeds.split(",") if item.strip()]

    seed_rows: list[dict] = []
    pair_rows: list[dict] = []
    scenario_aggregate: dict[str, list[float]] = {s.name: [] for s in SCENARIOS}

    for scenario in SCENARIOS:
        print(f"[scenario] {scenario.name}")
        for seed in seeds:
            print(f"  - seed={seed}")
            result = _run_scenario_seed(
                gcmc_cpu=gcmc_cpu,
                out_dir=out_dir,
                scenario=scenario,
                seed=seed,
                timeout_s=args.timeout,
            )
            seed_rows.append(result)
            pair_rows.extend(result["pair_rows"])
            scenario_aggregate[scenario.name].extend(
                r["delta_diff_d_minus_e"] for r in result["pair_rows"]
            )

    aggregate = {
        name: _stats(values)
        for name, values in scenario_aggregate.items()
    }

    summary = {
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "gcmc_cpu": str(gcmc_cpu),
        "output_dir": str(out_dir),
        "seeds": seeds,
        "scenarios": [s.name for s in SCENARIOS],
        "seed_results": [
            {
                k: v for k, v in row.items() if k != "pair_rows"
            }
            for row in seed_rows
        ],
        "aggregate_by_scenario": aggregate,
    }

    summary_path = out_dir / "summary.json"
    seed_csv_path = out_dir / "seed_stats.csv"
    pair_csv_path = out_dir / "pair_diffs.csv"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")
    _write_seed_csv(seed_csv_path, seed_rows)
    _write_pair_csv(pair_csv_path, pair_rows)

    print("\nReport written:")
    print(f"- {summary_path}")
    print(f"- {seed_csv_path}")
    print(f"- {pair_csv_path}")
    print("\nAggregate |D-E| stats:")
    for scenario_name, stats_dict in aggregate.items():
        print(
            f"- {scenario_name}: n={stats_dict['n']}, "
            f"mean={stats_dict['mean']}, max_abs={stats_dict['max_abs']}"
        )


if __name__ == "__main__":
    main()

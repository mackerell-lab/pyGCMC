"""
Three-way C/D/E paired deltaU audit runner.

Runs Mode C (pgp_full), Mode D (pme), and Mode E (pgp_full_pme) with
identical seeds/proposals and compares all pairwise deltaU differences.

Expected:
- E-C ~ +1.99 (insertion), -1.99 (deletion) — mesh-self offset
- D-C ~ +1.99 (insertion), -1.99 (deletion) — same mesh-self offset
- E-D ~ 0 (machine precision) — PGP+mesh-self matches PME exactly
"""

from __future__ import annotations

import argparse
import json
import subprocess
from collections import defaultdict
from pathlib import Path

MODES = ["pgp_full", "pme", "pgp_full_pme"]
MODE_LABELS = {"pgp_full": "C", "pme": "D", "pgp_full_pme": "E"}
PAIRS = [("pgp_full_pme", "pgp_full"), ("pme", "pgp_full"), ("pgp_full_pme", "pme")]  # E-C, D-C, E-D


def _find_gcmc_cpu(explicit_path: str | None) -> Path:
    if explicit_path:
        candidate = Path(explicit_path).resolve()
        if candidate.exists():
            return candidate
        raise FileNotFoundError(f"gcmc_cpu not found: {candidate}")
    fallback = Path("pygcmc_dev/build/bin/gcmc_cpu").resolve()
    if fallback.exists():
        return fallback
    raise FileNotFoundError("gcmc_cpu executable not found.")


def _load_records(path: Path) -> list[dict]:
    records: list[dict] = []
    for line in path.read_text().splitlines():
        if line.strip():
            records.append(json.loads(line))
    return records


def _stats(values: list[float]) -> dict:
    if not values:
        return {"n": 0}
    return {
        "n": len(values),
        "mean": sum(values) / len(values),
        "min": min(values),
        "max": max(values),
        "max_abs": max(abs(v) for v in values),
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


def _run_mode(*, gcmc_cpu: Path, inp_path: Path, run_dir: Path) -> None:
    run_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        str(gcmc_cpu),
        "--inp", str(inp_path),
        "--prefix", str(run_dir / "gcmc"),
        "--dump-accept", str(run_dir / "acceptance.jsonl"),
        "--dump-params", str(run_dir / "params.json"),
        "--print-freq", "1000000",
        "--traj-freq", "1000000",
        "--checkpoint-freq", "0",
        "--max-molecules-per-type", "-1",
    ]
    result = subprocess.run(cmd, cwd=str(Path.cwd()), capture_output=True, text=True, check=False)
    (run_dir / "stdout.log").write_text(result.stdout)
    (run_dir / "stderr.log").write_text(result.stderr)
    if result.returncode != 0:
        raise RuntimeError(
            f"gcmc_cpu failed for {inp_path}\nstdout:\n{result.stdout}\nstderr:\n{result.stderr}"
        )


def _pair_analysis(rows_a: list[dict], rows_b: list[dict], label: str) -> dict:
    """Compare rows_a - rows_b (e.g. E-C means rows_a=E, rows_b=C)."""
    pair_count = min(len(rows_a), len(rows_b))

    # Find first state divergence
    first_div_step: int | None = None
    first_div_reason: str | None = None
    for i in range(pair_count):
        ra, rb = rows_a[i], rows_b[i]
        if str(ra.get("move", "")).strip().lower() != str(rb.get("move", "")).strip().lower():
            first_div_step = i + 1
            first_div_reason = "move_mismatch"
            break
        if str(ra.get("species", "")).strip().lower() != str(rb.get("species", "")).strip().lower():
            first_div_step = i + 1
            first_div_reason = "species_mismatch"
            break
        if bool(ra.get("accepted", False)) != bool(rb.get("accepted", False)):
            first_div_step = i + 1
            first_div_reason = "accepted_mismatch"
            break
        if int(ra.get("nBefore", -1)) != int(rb.get("nBefore", -1)):
            first_div_step = i + 1
            first_div_reason = "n_before_mismatch"
            break

    prefix_limit = pair_count if first_div_step is None else max(0, first_div_step - 1)

    # Prefix-window stats (primary evidence)
    by_ms: dict[str, list[float]] = defaultdict(list)
    for i in range(prefix_limit):
        ra, rb = rows_a[i], rows_b[i]
        move = str(ra.get("move", "")).strip().lower()
        species = str(ra.get("species", "")).strip().lower()
        delta = float(ra.get("deltaU", 0.0)) - float(rb.get("deltaU", 0.0))
        by_ms[f"{move}:{species}"].append(delta)

    all_prefix = []
    for vals in by_ms.values():
        all_prefix.extend(vals)

    return {
        "label": label,
        "pair_count": pair_count,
        "first_divergence_step": first_div_step,
        "first_divergence_reason": first_div_reason,
        "prefix_steps": prefix_limit,
        "prefix_all": _stats(all_prefix),
        "prefix_sol_insertion": _stats(by_ms.get("insertion:sol", [])),
        "prefix_sol_deletion": _stats(by_ms.get("deletion:sol", [])),
        "prefix_by_move_species": {k: _stats(v) for k, v in sorted(by_ms.items())},
    }


def _summarize_seed(seed_dir: Path) -> dict:
    all_rows = {}
    for mode in MODES:
        path = seed_dir / mode / "acceptance.jsonl"
        all_rows[mode] = _load_records(path)

    pair_results = {}
    for mode_a, mode_b in PAIRS:
        label_a = MODE_LABELS[mode_a]
        label_b = MODE_LABELS[mode_b]
        label = f"{label_a}-{label_b}"
        pair_results[label] = _pair_analysis(all_rows[mode_a], all_rows[mode_b], label)

    return {
        "record_counts": {MODE_LABELS[m]: len(all_rows[m]) for m in MODES},
        "pairs": pair_results,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Three-way C/D/E paired deltaU audit.")
    parser.add_argument("--base-inp", default="tmp/pgp_mu_recenter_bracket_m14p0_seed101_ref/run.inp")
    parser.add_argument("--out-dir", default="tmp/pgp_cde_deltau_audit")
    parser.add_argument("--seeds", default="246810")
    parser.add_argument("--steps", type=int, default=500)
    parser.add_argument("--use-cavity-bias", choices=["yes", "no"], default="no")
    parser.add_argument("--temperature", type=float, default=300.0)
    parser.add_argument("--gcmc-cpu", default=None)
    args = parser.parse_args()

    base_inp = Path(args.base_inp).resolve()
    if not base_inp.exists():
        raise FileNotFoundError(f"Base INP not found: {base_inp}")
    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    gcmc_cpu = _find_gcmc_cpu(args.gcmc_cpu)
    seeds = [int(s.strip()) for s in args.seeds.split(",") if s.strip()]

    summary = {
        "config": {
            "base_inp": str(base_inp),
            "steps": args.steps,
            "use_cavity_bias": args.use_cavity_bias,
            "temperature": args.temperature,
            "seeds": seeds,
            "modes": {m: MODE_LABELS[m] for m in MODES},
        },
        "seeds": {},
    }

    for seed in seeds:
        print(f"=== Seed {seed} ===")
        seed_dir = out_dir / f"seed_{seed}"
        for mode in MODES:
            label = MODE_LABELS[mode]
            print(f"  Running Mode {label} ({mode})...")
            mode_dir = seed_dir / mode
            inp_path = mode_dir / "run.inp"
            _rewrite_inp(
                base_inp=base_inp, output_inp=inp_path, mode=mode,
                seed=seed, steps=args.steps,
                use_cavity_bias=args.use_cavity_bias, temperature=args.temperature,
            )
            _run_mode(gcmc_cpu=gcmc_cpu, inp_path=inp_path, run_dir=mode_dir)
        summary["seeds"][str(seed)] = _summarize_seed(seed_dir)

    out_path = out_dir / "summary.json"
    out_path.write_text(json.dumps(summary, indent=2) + "\n")
    print(f"\nSummary written to: {out_path}")

    # Print quick table
    for seed_key, seed_data in summary["seeds"].items():
        print(f"\n--- Seed {seed_key} ---")
        for pair_key, pair_data in seed_data["pairs"].items():
            sol_ins = pair_data.get("prefix_sol_insertion", {})
            sol_del = pair_data.get("prefix_sol_deletion", {})
            ins_mean = sol_ins.get("mean", "N/A")
            del_mean = sol_del.get("mean", "N/A")
            ins_n = sol_ins.get("n", 0)
            del_n = sol_del.get("n", 0)
            prefix = pair_data.get("prefix_steps", "?")
            print(f"  {pair_key}: prefix={prefix} steps")
            if isinstance(ins_mean, float):
                print(f"    SOL insertion: mean={ins_mean:+.6f}  (n={ins_n})")
            if isinstance(del_mean, float):
                print(f"    SOL deletion:  mean={del_mean:+.6f}  (n={del_n})")


if __name__ == "__main__":
    main()

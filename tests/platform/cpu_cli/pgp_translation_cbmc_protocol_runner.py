"""
Translation + CBMC protocol runner for Mode D/E paper validation.

Scope:
- Translation deltaU consistency: Mode D (pme) vs Mode E (pgp_full_pme)
- CBMC trial-energy consistency: compare relative trial energetics and
  normalized Boltzmann weights, not raw absolute energies
"""

from __future__ import annotations

import argparse
import json
import math
import subprocess
from pathlib import Path

MODES = ["pme", "pgp_full_pme"]
MODE_LABEL = {"pme": "D", "pgp_full_pme": "E"}


def _find_gcmc_cpu(explicit_path: str | None) -> Path:
    if explicit_path:
        candidate = Path(explicit_path).resolve()
        if candidate.exists():
            return candidate
        raise FileNotFoundError(f"gcmc_cpu not found: {candidate}")
    fallback = Path("pygcmc_dev/build/bin/gcmc_cpu").resolve()
    if fallback.exists():
        return fallback
    raise FileNotFoundError("gcmc_cpu executable not found. Build target gcmc_cpu first.")


def _parse_seed_csv(values: str) -> list[int]:
    parsed = [int(token.strip()) for token in values.split(",") if token.strip()]
    if not parsed:
        raise ValueError("At least one seed is required.")
    return parsed


def _load_jsonl(path: Path) -> list[dict]:
    rows: list[dict] = []
    for line in path.read_text().splitlines():
        if line.strip():
            rows.append(json.loads(line))
    return rows


def _stats(values: list[float]) -> dict:
    if not values:
        return {"n": 0}
    n = len(values)
    mean = sum(values) / n
    max_abs = max(abs(v) for v in values)
    rmse = math.sqrt(sum(v * v for v in values) / n)
    abs_sorted = sorted(abs(v) for v in values)
    p95 = abs_sorted[min(n - 1, int(round((n - 1) * 0.95)))]
    return {
        "n": n,
        "mean": mean,
        "mean_abs": sum(abs(v) for v in values) / n,
        "rmse": rmse,
        "max_abs": max_abs,
        "p95_abs": p95,
        "min": min(values),
        "max": max(values),
    }


def _rewrite_inp(
    *,
    base_inp: Path,
    out_inp: Path,
    mode: str,
    seed: int,
    steps: int,
    move_prob: str,
    use_conf_bias: str,
    use_cavity_bias: str,
    num_conf_bias_trial: int | None,
    temperature: float,
) -> None:
    base_dir = base_inp.parent
    path_like_keys = {"par", "fragitp", "atomtypes", "monomerdir", "top", "pdb", "protitp"}

    replacements = {
        "energy_method:": f"energy_method:{mode}",
        "random_seed:": f"random_seed:{seed}",
        "mcsteps:": f"mcsteps:{steps}",
        "nprint:": f"nprint:{steps}",
        "mc_move_prob:": f"mc_move_prob:{move_prob}",
        "use_conf_bias:": f"use_conf_bias:{use_conf_bias}",
        "use_cavity_bias:": f"use_cavity_bias:{use_cavity_bias}",
        "temperature:": f"temperature:{temperature:.8g}",
    }
    if num_conf_bias_trial is not None:
        replacements["num_conf_bias_trial:"] = f"num_conf_bias_trial:{num_conf_bias_trial}"

    seen: dict[str, bool] = {key: False for key in replacements}
    out_lines: list[str] = []
    for line in base_inp.read_text().splitlines():
        replaced = False
        for prefix, replacement in replacements.items():
            if line.startswith(prefix):
                out_lines.append(replacement)
                seen[prefix] = True
                replaced = True
                break
        if not replaced:
            if ":" in line:
                key, raw_value = line.split(":", 1)
                key_norm = key.strip().lower()
                value = raw_value.strip()
                if key_norm in path_like_keys and value:
                    candidate = Path(value)
                    if not candidate.is_absolute():
                        abs_candidate = (base_dir / candidate).resolve()
                        out_lines.append(f"{key.strip()}:{abs_candidate}")
                        continue
            out_lines.append(line)

    for prefix, replacement in replacements.items():
        if not seen[prefix]:
            out_lines.append(replacement)

    out_inp.parent.mkdir(parents=True, exist_ok=True)
    out_inp.write_text("\n".join(out_lines) + "\n")


def _run_mode(gcmc_cpu: Path, inp_path: Path, run_dir: Path, timeout_sec: int) -> None:
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
        timeout=timeout_sec,
    )
    (run_dir / "stdout.log").write_text(result.stdout)
    (run_dir / "stderr.log").write_text(result.stderr)
    if result.returncode != 0:
        raise RuntimeError(
            f"gcmc_cpu failed for {inp_path}\nstdout:\n{result.stdout}\nstderr:\n{result.stderr}"
        )


def _ranks(values: list[float]) -> list[float]:
    indexed = sorted((value, idx) for idx, value in enumerate(values))
    ranks = [0.0] * len(values)
    i = 0
    while i < len(indexed):
        j = i
        while j + 1 < len(indexed) and indexed[j + 1][0] == indexed[i][0]:
            j += 1
        rank = 0.5 * (i + j) + 1.0
        for k in range(i, j + 1):
            ranks[indexed[k][1]] = rank
        i = j + 1
    return ranks


def _spearman(a: list[float], b: list[float]) -> float:
    if len(a) != len(b) or len(a) < 2:
        return 1.0 if len(a) == len(b) else 0.0
    ra = _ranks(a)
    rb = _ranks(b)
    ma = sum(ra) / len(ra)
    mb = sum(rb) / len(rb)
    va = sum((x - ma) * (x - ma) for x in ra)
    vb = sum((x - mb) * (x - mb) for x in rb)
    if va <= 0.0 or vb <= 0.0:
        return 1.0
    cov = sum((ra[i] - ma) * (rb[i] - mb) for i in range(len(ra)))
    return cov / math.sqrt(va * vb)


def _softmax_weights(relative_energies: list[float], beta: float) -> list[float]:
    if not relative_energies:
        return []
    exps = [math.exp(-beta * value) for value in relative_energies]
    total = sum(exps)
    if total <= 0.0:
        return [0.0 for _ in exps]
    return [value / total for value in exps]


def _first_divergence(rows_d: list[dict], rows_e: list[dict]) -> tuple[int | None, str | None, int]:
    limit = min(len(rows_d), len(rows_e))
    for i in range(limit):
        d = rows_d[i]
        e = rows_e[i]
        if str(d.get("move", "")).strip().lower() != str(e.get("move", "")).strip().lower():
            return i + 1, "move_mismatch", max(0, i)
        if str(d.get("species", "")).strip().lower() != str(e.get("species", "")).strip().lower():
            return i + 1, "species_mismatch", max(0, i)
        if bool(d.get("accepted", False)) != bool(e.get("accepted", False)):
            return i + 1, "accepted_mismatch", max(0, i)
        if int(d.get("nBefore", -1)) != int(e.get("nBefore", -1)):
            return i + 1, "n_before_mismatch", max(0, i)
    return None, None, limit


def _summarize_translation(rows_d: list[dict], rows_e: list[dict]) -> dict:
    limit = min(len(rows_d), len(rows_e))
    deltas_diff: list[float] = []
    accepted_mismatch = 0
    used = 0
    for i in range(limit):
        d = rows_d[i]
        e = rows_e[i]
        if str(d.get("move", "")).strip().lower() != "translation":
            continue
        if str(e.get("move", "")).strip().lower() != "translation":
            continue
        used += 1
        deltas_diff.append(float(e.get("deltaU", 0.0)) - float(d.get("deltaU", 0.0)))
        if bool(d.get("accepted", False)) != bool(e.get("accepted", False)):
            accepted_mismatch += 1
    return {
        "paired_records": limit,
        "used_translation_pairs": used,
        "accepted_mismatch_count": accepted_mismatch,
        "delta_e_minus_d": _stats(deltas_diff),
    }


def _summarize_cbmc(rows_d: list[dict], rows_e: list[dict], temperature: float) -> dict:
    beta = 1.0 / (8.314e-3 * temperature)
    div_step, div_reason, prefix_limit = _first_divergence(rows_d, rows_e)

    top1_matches = 0
    compared = 0
    spearman_values: list[float] = []
    tvd_values: list[float] = []
    kl_values: list[float] = []
    logw_diff: list[float] = []
    selected_diff: list[float] = []
    trial_count_diff: list[float] = []

    for i in range(prefix_limit):
        d = rows_d[i]
        e = rows_e[i]
        if str(d.get("move", "")).strip().lower() != "insertion":
            continue
        if str(e.get("move", "")).strip().lower() != "insertion":
            continue
        trials_d = [float(x) for x in d.get("cbmcTrialEnergies", [])]
        trials_e = [float(x) for x in e.get("cbmcTrialEnergies", [])]
        if not trials_d or not trials_e:
            continue
        if len(trials_d) != len(trials_e):
            continue

        compared += 1
        rel_d = [value - min(trials_d) for value in trials_d]
        rel_e = [value - min(trials_e) for value in trials_e]

        idx_d = min(range(len(trials_d)), key=trials_d.__getitem__)
        idx_e = min(range(len(trials_e)), key=trials_e.__getitem__)
        if idx_d == idx_e:
            top1_matches += 1

        spearman_values.append(_spearman(rel_d, rel_e))

        wd = _softmax_weights(rel_d, beta)
        we = _softmax_weights(rel_e, beta)
        tvd_values.append(0.5 * sum(abs(wd[j] - we[j]) for j in range(len(wd))))

        eps = 1e-300
        kl_values.append(sum(wd[j] * math.log((wd[j] + eps) / (we[j] + eps)) for j in range(len(wd))))

        logw_diff.append(float(e.get("cbmcLogWOverK", 0.0)) - float(d.get("cbmcLogWOverK", 0.0)))
        selected_diff.append(float(e.get("cbmcSelectedEnergy", 0.0)) - float(d.get("cbmcSelectedEnergy", 0.0)))
        trial_count_diff.append(float(len(trials_e) - len(trials_d)))

    return {
        "pairing_window": {
            "first_divergence_step": div_step,
            "first_divergence_reason": div_reason,
            "prefix_steps": prefix_limit,
        },
        "cbmc_compared_attempts": compared,
        "top1_match_rate": (top1_matches / compared) if compared else 0.0,
        "spearman": _stats(spearman_values),
        "tvd": _stats(tvd_values),
        "kl_d_to_e": _stats(kl_values),
        "delta_logWOverK_e_minus_d": _stats(logw_diff),
        "delta_selected_energy_e_minus_d": _stats(selected_diff),
        "trial_count_diff_e_minus_d": _stats(trial_count_diff),
    }


def _run_phase(
    *,
    gcmc_cpu: Path,
    base_inp: Path,
    out_root: Path,
    seeds: list[int],
    mode_move_prob: str,
    use_conf_bias: str,
    num_conf_bias_trial: int | None,
    use_cavity_bias: str,
    steps: int,
    temperature: float,
    phase_name: str,
    run_timeout_sec: int,
) -> dict:
    phase_root = out_root / phase_name
    phase_root.mkdir(parents=True, exist_ok=True)
    results: dict[str, dict] = {}
    for seed in seeds:
        seed_dir = phase_root / f"seed_{seed}"
        for mode in MODES:
            print(f"[{phase_name}] seed={seed} mode={mode} steps={steps}")
            mode_dir = seed_dir / mode
            inp_path = mode_dir / "run.inp"
            _rewrite_inp(
                base_inp=base_inp,
                out_inp=inp_path,
                mode=mode,
                seed=seed,
                steps=steps,
                move_prob=mode_move_prob,
                use_conf_bias=use_conf_bias,
                use_cavity_bias=use_cavity_bias,
                num_conf_bias_trial=num_conf_bias_trial,
                temperature=temperature,
            )
            _run_mode(
                gcmc_cpu=gcmc_cpu,
                inp_path=inp_path,
                run_dir=mode_dir,
                timeout_sec=run_timeout_sec,
            )

        rows_d = _load_jsonl(seed_dir / "pme" / "acceptance.jsonl")
        rows_e = _load_jsonl(seed_dir / "pgp_full_pme" / "acceptance.jsonl")
        results[str(seed)] = {"rows_d": len(rows_d), "rows_e": len(rows_e)}
        if phase_name == "translation":
            results[str(seed)]["summary"] = _summarize_translation(rows_d, rows_e)
        else:
            results[str(seed)]["summary"] = _summarize_cbmc(rows_d, rows_e, temperature)
    return results


def _aggregate_translation(seed_results: dict[str, dict]) -> dict:
    merged: list[float] = []
    total_pairs = 0
    mismatch = 0
    for result in seed_results.values():
        summary = result["summary"]
        total_pairs += int(summary["used_translation_pairs"])
        mismatch += int(summary["accepted_mismatch_count"])
        # Rebuild from summary is not possible; keep aggregate from per-seed mean not robust.
        # For aggregation, only report pair counts and mismatch rate.
    return {
        "total_translation_pairs": total_pairs,
        "total_accepted_mismatch": mismatch,
        "accepted_mismatch_rate": (mismatch / total_pairs) if total_pairs else 0.0,
        "note": "Use per-seed delta_e_minus_d statistics as primary quantitative evidence.",
    }


def _aggregate_cbmc(seed_results: dict[str, dict]) -> dict:
    compared = 0
    weighted_top1 = 0.0
    prefix_steps = 0
    for result in seed_results.values():
        summary = result["summary"]
        n = int(summary["cbmc_compared_attempts"])
        compared += n
        weighted_top1 += float(summary["top1_match_rate"]) * n
        prefix_steps += int(summary["pairing_window"]["prefix_steps"])
    return {
        "total_cbmc_compared_attempts": compared,
        "pooled_top1_match_rate": (weighted_top1 / compared) if compared else 0.0,
        "sum_prefix_steps": prefix_steps,
        "note": "Use per-seed spearman/tvd/kl statistics as primary quantitative evidence.",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Run protocol tests for translation + CBMC (Mode D vs E).")
    parser.add_argument("--base-inp", default="tmp/pgp_mu_recenter_bracket_m14p0_seed101_ref/run.inp")
    parser.add_argument("--out-dir", default="tmp/pgp_translation_cbmc_protocol")
    parser.add_argument("--seeds", default="246810,97531")
    parser.add_argument("--translation-steps", type=int, default=600)
    parser.add_argument("--cbmc-steps", type=int, default=240)
    parser.add_argument("--cbmc-trials", type=int, default=10)
    parser.add_argument("--run-timeout-sec", type=int, default=360)
    parser.add_argument("--temperature", type=float, default=300.0)
    parser.add_argument("--use-cavity-bias", choices=["yes", "no"], default="no")
    parser.add_argument("--gcmc-cpu", default=None)
    args = parser.parse_args()

    base_inp = Path(args.base_inp).resolve()
    if not base_inp.exists():
        raise FileNotFoundError(f"Base INP not found: {base_inp}")

    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    gcmc_cpu = _find_gcmc_cpu(args.gcmc_cpu)
    seeds = _parse_seed_csv(args.seeds)

    translation_results = _run_phase(
        gcmc_cpu=gcmc_cpu,
        base_inp=base_inp,
        out_root=out_dir,
        seeds=seeds,
        mode_move_prob="0.000000 0.000000 1.000000 0.000000",
        use_conf_bias="no",
        num_conf_bias_trial=None,
        use_cavity_bias=args.use_cavity_bias,
        steps=args.translation_steps,
        temperature=args.temperature,
        phase_name="translation",
        run_timeout_sec=args.run_timeout_sec,
    )

    cbmc_results = _run_phase(
        gcmc_cpu=gcmc_cpu,
        base_inp=base_inp,
        out_root=out_dir,
        seeds=seeds,
        mode_move_prob="1.000000 0.000000 0.000000 0.000000",
        use_conf_bias="yes",
        num_conf_bias_trial=args.cbmc_trials,
        use_cavity_bias=args.use_cavity_bias,
        steps=args.cbmc_steps,
        temperature=args.temperature,
        phase_name="cbmc",
        run_timeout_sec=args.run_timeout_sec,
    )

    summary = {
        "config": {
            "base_inp": str(base_inp),
            "gcmc_cpu": str(gcmc_cpu),
            "seeds": seeds,
            "temperature": args.temperature,
            "use_cavity_bias": args.use_cavity_bias,
            "translation_steps": args.translation_steps,
            "cbmc_steps": args.cbmc_steps,
            "cbmc_trials": args.cbmc_trials,
            "run_timeout_sec": args.run_timeout_sec,
            "modes": MODE_LABEL,
        },
        "translation": {
            "per_seed": translation_results,
            "aggregate": _aggregate_translation(translation_results),
        },
        "cbmc": {
            "per_seed": cbmc_results,
            "aggregate": _aggregate_cbmc(cbmc_results),
        },
    }

    out_path = out_dir / "summary.json"
    out_path.write_text(json.dumps(summary, indent=2) + "\n")
    print(f"Protocol summary written to: {out_path}")


if __name__ == "__main__":
    main()

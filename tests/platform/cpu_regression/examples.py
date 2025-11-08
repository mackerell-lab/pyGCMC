"""
High-level regression tests that reuse the sample GCMC INP decks under
tests/data/gcmc_examples. Each test launches the gcmc_cpu CLI with those inputs,
captures the acceptance log, and verifies physical consistency.
"""

from __future__ import annotations

import json
import math
import subprocess
from pathlib import Path
from typing import Dict, List

import pytest


GCMC_CPU = Path(__file__).resolve().parents[3] / "build" / "bin" / "gcmc_cpu"
DATA_DIR = Path(__file__).resolve().parents[2] / "data"
TEMPLATE_DIR = DATA_DIR / "gcmc_examples"


def _ensure_binary() -> None:
    if not GCMC_CPU.exists():
        pytest.skip(f"gcmc_cpu not built: {GCMC_CPU}")


def _render_template(template_name: str, tmp_path: Path) -> Path:
    template_path = TEMPLATE_DIR / template_name
    content = template_path.read_text().replace("{DATA_DIR}", str(DATA_DIR))
    rendered = tmp_path / template_name.replace("_template", "")
    rendered.write_text(content)
    return rendered


def _run_gcmc(tmp_path: Path, template_name: str, prefix: str) -> Dict[str, Path]:
    _ensure_binary()
    inp_path = _render_template(template_name, tmp_path)
    accept_log = tmp_path / f"{prefix}_accept.jsonl"
    prefix_path = tmp_path / prefix
    cmd = [
        str(GCMC_CPU),
        "--inp",
        str(inp_path),
        "--prefix",
        str(prefix_path),
        "--dump-accept",
        str(accept_log),
        "--store-probabilities",
    ]
    result = subprocess.run(
        cmd,
        cwd=tmp_path,
        text=True,
        capture_output=True,
        timeout=180,
    )
    if result.returncode != 0:
        raise AssertionError(
            "gcmc_cpu failed\n"
            f"STDOUT:\n{result.stdout}\n"
            f"STDERR:\n{result.stderr}"
        )
    stats_path = tmp_path / f"{prefix}_statistics.dat"
    assert accept_log.exists(), "acceptance log missing"
    assert stats_path.exists(), "statistics file missing"
    return {"accept_log": accept_log, "stats": stats_path}


def _load_acceptance_records(path: Path) -> List[Dict]:
    with path.open() as fh:
        return [json.loads(line) for line in fh if line.strip()]


def _probability_diffs(records: List[Dict], max_samples: int = 400) -> List[float]:
    diffs: List[float] = []
    for rec in records:
        move = rec.get("move")
        if move not in {"insertion", "deletion"}:
            continue
        actual = rec.get("pAcc", -1.0)
        if not (0.0 <= actual <= 1.0):
            continue

        z = rec.get("z", 1.0)
        v_eff = rec.get("vEff", rec.get("vBox", 1.0))
        n_before = rec.get("nBefore", 0)
        beta_delta = rec.get("betaDeltaU", 0.0)
        proposal_ratio = rec.get("proposalRatio", 1.0) or 1.0
        if move == "insertion":
            if n_before + 1 <= 0:
                continue
            q_forward = rec.get("qForward", 1.0)
            expected = (
                z * v_eff / (n_before + 1.0) * math.exp(-beta_delta) * q_forward * proposal_ratio
            )
        else:
            if n_before <= 0 or z * v_eff <= 0:
                continue
            q_reverse = rec.get("qReverse", 1.0) or 1.0
            expected = (
                n_before / (z * v_eff) * math.exp(beta_delta) / q_reverse * proposal_ratio
            )
        diffs.append(abs(min(1.0, expected) - actual))
        if len(diffs) >= max_samples:
            break
    return diffs


def _read_last_stats(stats_path: Path) -> Dict[str, float]:
    last_line = None
    with stats_path.open() as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            last_line = line
    assert last_line, "statistics file was empty"
    parts = last_line.split()
    return {
        "step": int(parts[0]),
        "n_total": int(parts[2]),
        "ins_attempts": int(parts[4]),
        "ins_accepted": int(parts[5]),
        "del_attempts": int(parts[6]),
        "del_accepted": int(parts[7]),
    }


def _species_frequencies(records: List[Dict]) -> Dict[str, float]:
    counts: Dict[str, int] = {}
    total = 0
    for rec in records:
        if rec.get("move") not in {"insertion", "deletion"}:
            continue
        species = rec.get("species", "unknown")
        counts[species] = counts.get(species, 0) + 1
        total += 1
    if total == 0:
        return {}
    return {name: cnt / total for name, cnt in counts.items()}


class TestGCMCCPUExamples:
    """Integration-style checks that reuse curated gcmc_cpu example inputs."""

    def test_water_tip3p_probabilities(self, tmp_path):
        paths = _run_gcmc(tmp_path, "water_tip3p_template.inp", "water_case")
        records = _load_acceptance_records(paths["accept_log"])
        diffs = _probability_diffs(records, max_samples=200)
        assert diffs, "No insertion/deletion records captured"
        assert max(diffs) < 5e-4

        stats = _read_last_stats(paths["stats"])
        assert 5 <= stats["n_total"] <= 60, f"N_total out of range: {stats['n_total']}"

    def test_multi_salt_distribution_and_probabilities(self, tmp_path):
        paths = _run_gcmc(tmp_path, "multi_salt_template.inp", "multi_case")
        records = _load_acceptance_records(paths["accept_log"])

        diffs = _probability_diffs(records, max_samples=400)
        assert diffs, "Missing GCMC attempts in multi-component log"
        assert max(diffs) < 5e-4

        freqs = _species_frequencies(records)
        water_freq = freqs.get("WAT") or freqs.get("water")
        assert water_freq is not None, "Water species missing from log"
        assert abs(water_freq - 0.8) < 0.1
        for ion in ("NA", "CL"):
            assert ion in freqs, f"{ion} species missing from log"
            assert abs(freqs[ion] - 0.1) < 0.1

    def test_benzene_balances_insert_delete(self, tmp_path):
        paths = _run_gcmc(tmp_path, "benzene_template.inp", "benz_case")
        records = _load_acceptance_records(paths["accept_log"])
        diffs = _probability_diffs(records, max_samples=300)
        assert diffs, "No benzene moves recorded"
        assert max(diffs) < 5e-4

        stats = _read_last_stats(paths["stats"])
        ins = stats["ins_accepted"]
        dels = stats["del_accepted"]
        assert ins > 0 and dels > 0, "Benzene run did not attempt both move types"
        imbalance = abs(ins - dels) / max(1, ins + dels)
        assert imbalance < 0.15, f"Insertion/deletion imbalance too high: {imbalance:.2f}"

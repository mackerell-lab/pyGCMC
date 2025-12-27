"""
Regression: asymmetric insertion/deletion proposal probabilities must enter MH acceptance.

This test targets a common P0 bug pattern:
- Tests that only use symmetric p(ins)==p(del) cannot detect missing MH correction.
- Here we make p(ins)!=p(del) via attempt_prob_{ins,del} and assert that:
  1) JSONL `proposalRatio` matches the expected p_rev/p_fwd for each move type
  2) The logged `pAcc` can be reconstructed from the structured fields (no log parsing)
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

from .inp_units_compat_helpers import _run_gcmc_cpu, _write_inp


def _read_accept_records(path: Path) -> list[dict]:
    return [json.loads(line) for line in path.read_text().splitlines() if line.strip()]


def test_asymmetric_attempt_prob_enters_mh_proposal_ratio_and_pacc(
    gcmc_cpu, test_data_dir, temp_dir
):
    itp = test_data_dir / "charmm36.ff" / "mol" / "na.itp"
    assert itp.exists()

    work = Path(temp_dir) / "asymmetric_proposal_bias"
    work.mkdir(parents=True, exist_ok=True)

    out_prefix = work / "out" / "gcmc"
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    accept_log = work / "out" / "acceptance.jsonl"

    # Choose a small box so the run stays fast, but large enough to allow multiple particles.
    # INP is gcmc_gpu/opencl style: Å for lengths and kcal/mol for fragmuex.
    _write_inp(
        work / "run.inp",
        f"""
random_seed:123
fragitp:{itp}
fragname:NA
fragconc:1.0
fragmuex:0.0

box_size:20.0 20.0 20.0
cutoff:10.0
temperature:300.0
moves_per_step:1
mcsteps:2500
nprint:2500

use_cavity_bias:no
use_conf_bias:no

attempt_prob_ins:0.9
attempt_prob_del:0.1
attempt_prob_trn:0.0
attempt_prob_rot:0.0
""",
    )

    result = _run_gcmc_cpu(
        gcmc_cpu,
        workdir=work,
        inp=work / "run.inp",
        out_prefix=out_prefix,
        extra_args=["--dump-accept", str(accept_log)],
        timeout=60,
    )
    assert result.returncode == 0, result.stdout + result.stderr

    records = _read_accept_records(accept_log)
    assert records, f"No acceptance records found in {accept_log}"

    # Filter to a state where both insertion and deletion are allowed.
    ins = [
        r
        for r in records
        if r.get("move") == "insertion" and r.get("species") == "NA" and int(r.get("nBefore", 0)) > 0
    ]
    dels = [r for r in records if r.get("move") == "deletion" and r.get("species") == "NA"]

    assert ins, "Need insertion records with nBefore>0 (so deletion is also allowed)"
    assert dels, "Need deletion records"

    expected_ins_ratio = 0.1 / 0.9  # p(delete)/p(insert)
    expected_del_ratio = 0.9 / 0.1  # p(insert)/p(delete)

    # Use a modest tolerance: proposalRatio is computed from double CDF values.
    for rec in ins[:50]:
        assert float(rec["proposalRatio"]) == pytest.approx(expected_ins_ratio, rel=1e-6, abs=1e-12)
    for rec in dels[:50]:
        assert float(rec["proposalRatio"]) == pytest.approx(expected_del_ratio, rel=1e-6, abs=1e-12)

    # Reconstruct pAcc from the JSONL fields (behavior-driven, log-independent).
    def pacc_from_record(rec: dict) -> float:
        proposal = float(rec.get("proposalRatio", 1.0))
        beta_delta_u = float(rec.get("betaDeltaU", 0.0))
        z = float(rec.get("z", 1.0))
        v_eff = float(rec.get("vEff", 1.0))
        n_before = int(rec.get("nBefore", 0))
        cavity = float(rec.get("cavityFraction", 1.0))
        rosen = float(rec.get("rosenbluthWeight", 1.0))

        if rec.get("move") == "insertion":
            log_ratio = (
                math.log(max(proposal, 1e-300))
                - beta_delta_u
                + math.log(max(z, 1e-300))
                + math.log(max(v_eff, 1e-300))
                + math.log(max(cavity, 1e-300))
                - math.log(n_before + 1.0)
                + math.log(max(rosen, 1e-300))
            )
        else:
            log_ratio = (
                math.log(max(proposal, 1e-300))
                - beta_delta_u
                - math.log(max(z, 1e-300))
                + math.log(max(n_before, 1e-300))
                - math.log(max(v_eff, 1e-300))
                - math.log(max(cavity, 1e-300))
                - math.log(max(rosen, 1e-300))
            )

        if log_ratio >= 0.0:
            return 1.0
        return math.exp(max(log_ratio, -700.0))

    for rec in (ins[:50] + dels[:50]):
        expected = pacc_from_record(rec)
        actual = float(rec.get("pAcc", -1.0))
        assert actual >= 0.0
        assert actual == pytest.approx(expected, rel=1e-10, abs=1e-12)

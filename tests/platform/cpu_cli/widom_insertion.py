"""
Widom test-particle insertion cross-check for interacting systems (μVT contract).

For a single species in the grand canonical ensemble, the exact identity holds:

    ρ = z * < exp(-β ΔU_insert) >

where:
- ρ = <N>/V is the number density,
- z is the activity used by the acceptance formula,
- ΔU_insert is the energy change for inserting a test particle at a uniformly random position.

We estimate the Widom factor using the --dump-accept JSONL stream from a normal GCMC run:
each insertion attempt provides ΔU for a uniformly proposed insertion (no cavity/CBMC bias here).

Assertions are file-driven only (no stdout/stderr log matching).
"""

from __future__ import annotations

import json
import math
import statistics
from pathlib import Path

import pytest

from .inp_units_compat_helpers import _run_gcmc_cpu, _write_inp


def _write_text(path: Path, content: str) -> None:
    path.write_text(content.strip() + "\n")


def _load_accept_records(path: Path) -> list[dict]:
    return [json.loads(line) for line in path.read_text().splitlines() if line.strip()]


def test_widom_insertion_mu_rho_consistency_interacting_lj(gcmc_cpu, temp_dir):
    work = Path(temp_dir) / "widom_insertion"
    work.mkdir(parents=True, exist_ok=True)

    par = work / "par.itp"
    frag = work / "frag.itp"
    _write_text(
        par,
        """
[ defaults ]
1 2 yes 0.5 0.8333

[ atomtypes ]
; name  at.num  mass   charge  ptype  sigma   epsilon
X       0       1.000  0.000   A      0.250   0.200
""",
    )
    _write_text(
        frag,
        """
[ moleculetype ]
X   1

[ atoms ]
; nr  type  resnr  residue  atom  cgnr  charge  mass
1   X     1      X        X     1     0.000   1.000
""",
    )

    out_prefix = work / "out" / "gcmc"
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    accept_log = work / "out" / "acceptance.jsonl"

    # Use default gcmc_gpu INP semantics (Å, kcal) for box/cutoff; forcefield remains GROMACS (nm, kJ).
    # With fragmuex=0, activity z is set by the ideal-gas concentration: z = conc[M] * 0.6022 (nm^-3).
    inp = work / "run.inp"
    _write_inp(
        inp,
        f"""
random_seed:20251229
par:{par}
fragitp:{frag}
fragname:X
fragconc:2.0
fragmuex:0.0

box_size:20.0 20.0 20.0
gcmc_region:box 0 0 0 20 20 20
cutoff:10.0
temperature:300.0
moves_per_step:1
mcsteps:4000
nprint:1

attempt_prob_ins:0.5
attempt_prob_del:0.5
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
        timeout=60,
    )
    assert result.returncode == 0, result.stdout + result.stderr

    records = _load_accept_records(accept_log)
    insertions = [
        r
        for r in records
        if str(r.get("move", "")).strip().lower() == "insertion"
        and str(r.get("species", "")).strip().upper() == "X"
    ]
    assert len(insertions) >= 800, f"Insufficient insertion attempts for Widom estimate: {len(insertions)}"

    # Guardrails: this test assumes unbiased uniform insertion proposals (no CBMC/cavity bias).
    assert max(int(r.get("cbmcTrials", 0)) for r in insertions) <= 1
    for r in insertions[:10]:
        assert float(r.get("qForward", 1.0)) == pytest.approx(1.0, abs=1e-12)
        assert float(r.get("qReverse", 1.0)) == pytest.approx(1.0, abs=1e-12)
        assert float(r.get("proposalRatio", 1.0)) == pytest.approx(1.0, abs=1e-12)
        assert float(r.get("vEff", r.get("vBox", 0.0))) == pytest.approx(float(r["vBox"]), abs=1e-12)

    burnin = len(insertions) // 5
    insertions = insertions[burnin:]
    assert len(insertions) >= 500

    beta = float(insertions[0]["beta"])
    z = float(insertions[0]["z"])
    v_box = float(insertions[0]["vBox"])
    assert beta > 0.0
    assert z > 0.0
    assert v_box > 0.0

    widom_weights = [math.exp(-beta * float(r["deltaU"])) for r in insertions]
    widom_factor = statistics.mean(widom_weights)
    rho_widom = z * widom_factor

    n_values = [int(r["nBefore"]) for r in insertions]
    mean_n = statistics.mean(n_values)
    rho_measured = mean_n / v_box

    # In the μVT ensemble this identity is exact; we allow finite-sample noise with a moderate tolerance.
    assert rho_measured == pytest.approx(rho_widom, rel=0.25, abs=0.05)


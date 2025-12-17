"""
INP compatibility tests for gcmc_gpu-style units (fragmuex acceptance scaling).
"""

from __future__ import annotations

import math
from pathlib import Path

import pytest

from .inp_units_compat_helpers import _first_accept_record, _run_gcmc_cpu, _write_inp


def test_fragmuex_scales_activity_and_acceptance_in_ideal_gas_limit(
    gcmc_cpu, test_data_dir, temp_dir
):
    """
    Behavior-driven check that fragmuex influences acceptance probability via activity.

    Use an empty 1 nm^3 box with a single-atom fragment so deltaU ~= 0 and:
        pAcc = min(1, z * V / (nBefore+1))
    For small z (z<1), this reduces to pAcc == z, and changing μ should scale z by exp(beta*μ).
    """
    itp = test_data_dir / "charmm36.ff" / "mol" / "na.itp"
    assert itp.exists()

    conc_m = 0.1
    temperature_k = 300.0

    def run_once(muex_kcal_mol: float) -> dict:
        tag = "mu0" if muex_kcal_mol == 0.0 else "mu_nonzero"
        work = Path(temp_dir) / "muex_pacc_effect" / tag
        work.mkdir(parents=True, exist_ok=True)

        out_prefix = work / "out" / "gcmc"
        out_prefix.parent.mkdir(parents=True, exist_ok=True)
        accept_log = work / "out" / "acceptance.jsonl"

        inp = work / "test.inp"
        _write_inp(
            inp,
            f"""
random_seed:123
fragitp:{itp}
fragname:NA
fragconc:{conc_m}
fragmuex:{muex_kcal_mol}

box_size:10.0 10.0 10.0
cutoff:4.0
temperature:{temperature_k}
moves_per_step:1
mcsteps:1
nprint:1
mc_move_prob:1 0 0 0
""",
        )

        result = _run_gcmc_cpu(
            gcmc_cpu,
            workdir=work,
            inp=inp,
            out_prefix=out_prefix,
            extra_args=["--dump-accept", str(accept_log)],
            timeout=30,
        )
        assert result.returncode == 0, result.stdout + result.stderr

        rec = _first_accept_record(accept_log, move="insertion", species="NA")
        return rec

    rec_mu0 = run_once(0.0)
    rec_mu1 = run_once(-1.0)

    # Ideal-gas-ish: in an empty box, insertion should have no interaction energy.
    assert float(rec_mu0["deltaU"]) == pytest.approx(0.0, abs=1e-12)
    assert float(rec_mu1["deltaU"]) == pytest.approx(0.0, abs=1e-12)
    assert int(rec_mu0["nBefore"]) == 0
    assert int(rec_mu1["nBefore"]) == 0
    assert float(rec_mu0["vBox"]) == pytest.approx(1.0, abs=1e-12)
    assert float(rec_mu1["vBox"]) == pytest.approx(1.0, abs=1e-12)
    assert float(rec_mu0["cavityFraction"]) == pytest.approx(1.0, abs=1e-12)
    assert float(rec_mu1["cavityFraction"]) == pytest.approx(1.0, abs=1e-12)
    assert float(rec_mu0["rosenbluthWeight"]) == pytest.approx(1.0, abs=1e-12)
    assert float(rec_mu1["rosenbluthWeight"]) == pytest.approx(1.0, abs=1e-12)
    assert float(rec_mu0["proposalRatio"]) == pytest.approx(1.0, abs=1e-12)
    assert float(rec_mu1["proposalRatio"]) == pytest.approx(1.0, abs=1e-12)

    # Activity should follow z = conc * (M->nm^-3) * exp(beta*mu).
    # (The implementation uses a rounded conversion constant, so keep tolerances loose for absolute z.)
    beta = 1.0 / (8.314e-3 * temperature_k)  # mol/kJ
    expected_z0 = conc_m * 0.6022
    assert float(rec_mu0["mu"]) == pytest.approx(0.0, abs=1e-12)
    assert float(rec_mu0["z"]) == pytest.approx(expected_z0, rel=1e-4, abs=1e-12)

    expected_mu1_kj = -4.184  # -1.0 kcal/mol -> kJ/mol
    assert float(rec_mu1["mu"]) == pytest.approx(expected_mu1_kj, rel=1e-6, abs=1e-6)
    expected_z1 = float(rec_mu0["z"]) * math.exp(beta * expected_mu1_kj)
    assert float(rec_mu1["z"]) == pytest.approx(expected_z1, rel=1e-6, abs=1e-12)

    # In this limit, pAcc should reduce to z*V/(n+1). With V=1 and n=0, pAcc==z.
    assert float(rec_mu0["pAcc"]) == pytest.approx(float(rec_mu0["z"]), rel=1e-12, abs=1e-12)
    assert float(rec_mu1["pAcc"]) == pytest.approx(float(rec_mu1["z"]), rel=1e-12, abs=1e-12)

    # And changing μ should scale pAcc by exp(beta*Δμ) (and here Δμ == μ1 since μ0==0).
    pacc_ratio = float(rec_mu1["pAcc"]) / float(rec_mu0["pAcc"])
    assert float(rec_mu0["pAcc"]) < 0.2
    assert float(rec_mu1["pAcc"]) < float(rec_mu0["pAcc"])
    assert pacc_ratio == pytest.approx(math.exp(beta * float(rec_mu1["mu"])), rel=1e-6, abs=1e-12)

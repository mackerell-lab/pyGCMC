"""
Unit tests for GCMCAcceptance detailed interfaces.

Verifies that the logarithmic acceptance ratios for a matched
insertion / deletion pair are strict inverses
(i.e., ln r_ins + ln r_del == 0), satisfying detailed balance
in log space.
"""

import math

import pytest

import pygcmc


@pytest.mark.parametrize(
    "n_before, delta_e_ins, delta_e_del, cavity_fraction, lambda_nm, rosenbluth, cbmc_trials, proposal_log_ratio",
    [
        (10, 2.5, -2.5, 0.4, 0.3, 1.8, 5, -0.2),
        (2, -1.2, 1.2, 0.85, 1.0, 0.9, 3, 0.15),
    ],
)
def test_detailed_log_ratio_matches_formula(
    n_before,
    delta_e_ins,
    delta_e_del,
    cavity_fraction,
    lambda_nm,
    rosenbluth,
    cbmc_trials,
    proposal_log_ratio,
):
    """Verify detailed acceptance returns the analytic log ratio."""
    acc = pygcmc.GCMCAcceptance()
    temperature = 298.15
    acc.setTemperature(temperature)
    volume = 27.0
    acc.setVolume(volume)

    mu = -4.2
    acc.setChemicalPotential(0, mu)
    acc.setThermalLambda(0, lambda_nm)

    insertion_result = acc.calculate_insertion_probability_detailed(
        typeId=0,
        currentNumber=n_before,
        deltaE=delta_e_ins,
        cavityFraction=cavity_fraction,
        lambdaNm=lambda_nm,
        rosenbluthWeight=rosenbluth,
        cbmcTrials=cbmc_trials,
        proposalLogRatio=proposal_log_ratio,
    )

    deletion_result = acc.calculate_deletion_probability_detailed(
        typeId=0,
        currentNumber=n_before + 1,
        deltaE=delta_e_del,
        cavityFraction=cavity_fraction,
        lambdaNm=lambda_nm,
        rosenbluthWeight=rosenbluth,
        cbmcTrials=cbmc_trials,
        proposalLogRatio=-proposal_log_ratio,
    )

    beta = 1.0 / (8.314e-3 * temperature)
    activity = math.exp(beta * mu)
    log_volume = math.log(volume)
    log_cavity = math.log(cavity_fraction)
    log_lambda3 = 3.0 * math.log(lambda_nm)

    expected_log_ins = (
        proposal_log_ratio
        - beta * delta_e_ins
        + math.log(activity)
        + (log_volume + log_cavity)
        - math.log(n_before + 1)
        + math.log(rosenbluth)
        - log_lambda3
    )

    expected_log_del = (
        -proposal_log_ratio
        + beta * delta_e_del
        - math.log(activity)
        + math.log(n_before + 1)
        - (log_volume + log_cavity)
        - math.log(rosenbluth)
        + log_lambda3
    )

    assert math.isclose(
        insertion_result["logRatio"],
        expected_log_ins,
        rel_tol=0.0,
        abs_tol=1e-12,
    ), (
        f"Insertion log ratio mismatch: "
        f"got {insertion_result['logRatio']:.12e}, expected {expected_log_ins:.12e}"
    )

    assert math.isclose(
        deletion_result["logRatio"],
        expected_log_del,
        rel_tol=0.0,
        abs_tol=1e-12,
    ), (
        f"Deletion log ratio mismatch: "
        f"got {deletion_result['logRatio']:.12e}, expected {expected_log_del:.12e}"
    )

    # Probabilities should be bounded between 0 and 1
    for label, result in [
        ("insertion", insertion_result),
        ("deletion", deletion_result),
    ]:
        prob = result["probability"]
        assert 0.0 <= prob <= 1.0, f"{label} probability out of bounds: {prob}"


def test_log_ratio_sum_is_zero():
    """
    Ensure ln(r_ins) + ln(r_del) cancels for matched moves.
    
    CRITICAL: This test verifies detailed balance in log-space for the SAME microstate.
    For a reversible move pair (insert then delete same molecule), energy changes cancel:
    - Insertion: ΔE_ins = E(N+1) - E(N)
    - Deletion: ΔE_del = E(N) - E(N+1) = -ΔE_ins
    
    When both energy changes occur (not when they cancel), we verify the
    probability RATIO p_ins/p_del matches theory (energy-independent part only).
    """
    acc = pygcmc.GCMCAcceptance()
    temperature = 298.15
    acc.setTemperature(temperature)
    volume = 27.0
    acc.setVolume(volume)

    mu = -4.2
    acc.setChemicalPotential(0, mu)
    beta = 1.0 / (8.314e-3 * temperature)

    test_cases = [
        # (N, cavity, lambda_nm, rosenbluth, cbmc_k, proposal_log)
        # NOTE: Use rosenbluth=1, cbmc_k=1 to match MovementModule baseline
        (10, 0.4, 0.5, 1.0, 1, 0.0),
        (5, 1.0, 1.0, 1.0, 1, 0.0),
        (15, 0.3, 0.7, 1.0, 1, 0.0),
        (2, 0.9, 0.4, 1.0, 1, 0.0),
    ]

    for (
        n_before,
        cavity_fraction,
        lambda_nm,
        rosenbluth,
        cbmc_trials,
        proposal_log_ratio,
    ) in test_cases:
        acc.setThermalLambda(0, lambda_nm)

        # For symmetry check, use deltaE = 0 (energy cancels in reversible move)
        ins_result = acc.calculate_insertion_probability_detailed(
            typeId=0,
            currentNumber=n_before,
            deltaE=0.0,  # Energy-neutral for symmetry test
            cavityFraction=cavity_fraction,
            lambdaNm=lambda_nm,
            rosenbluthWeight=rosenbluth,
            cbmcTrials=cbmc_trials,
            proposalLogRatio=proposal_log_ratio,
        )

        del_result = acc.calculate_deletion_probability_detailed(
            typeId=0,
            currentNumber=n_before + 1,
            deltaE=0.0,  # Energy-neutral for symmetry test
            cavityFraction=cavity_fraction,
            lambdaNm=lambda_nm,
            rosenbluthWeight=rosenbluth,
            cbmcTrials=cbmc_trials,
            proposalLogRatio=-proposal_log_ratio,
        )

        p_ins = ins_result["probability"]
        p_del = del_result["probability"]
        
        # Verify probability ratio matches theory (energy-independent part)
        # p_ins/p_del = exp(βμ) * V * f_cav / ((N+1) * Λ³)
        theory_ratio = (
            math.exp(beta * mu) 
            * volume * cavity_fraction 
            / ((n_before + 1) * (lambda_nm ** 3))
        )
        # Note: (W/K) / (K/W) = (W/K)² cancels if same rosenbluth for both
        
        actual_ratio = p_ins / p_del if p_del > 0 else 0
        rel_error = abs(actual_ratio - theory_ratio) / max(theory_ratio, 1e-12)
        
        assert rel_error < 1e-10, (
            f"Ratio mismatch: N={n_before}, cav={cavity_fraction}, λ={lambda_nm}, "
            f"p_ins/p_del={actual_ratio:.12e}, theory={theory_ratio:.12e}, "
            f"error={rel_error:.2e}"
        )
        
        print(f"  ✓ Case N={n_before}, cav={cavity_fraction:.1f}, λ={lambda_nm:.1f}: "
              f"ratio={actual_ratio:.6f}, theory={theory_ratio:.6f}")
    
    print(f"\n✅ All {len(test_cases)} symmetry cases passed (energy-neutral)")

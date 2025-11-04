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

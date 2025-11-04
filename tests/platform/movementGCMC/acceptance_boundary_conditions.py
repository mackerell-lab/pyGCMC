"""
Boundary condition tests for grand canonical acceptance calculations.
"""

import pytest

import pygcmc


def _make_acceptance():
    acc = pygcmc.GCMCAcceptance()
    acc.setChemicalPotential(0, -2.0)
    acc.setThermalLambda(0, 1.0)
    return acc


def test_zero_cavity_rejects_insertion():
    """Insertion probability should vanish when cavity fraction is ~0."""
    acc = _make_acceptance()
    acc.setTemperature(300.0)
    acc.setVolume(50.0)

    result = acc.calculate_insertion_probability_detailed(
        typeId=0,
        currentNumber=5,
        deltaE=-10.0,
        cavityFraction=1e-12,
        lambdaNm=1.0,
        rosenbluthWeight=1.0,
        cbmcTrials=1,
        proposalLogRatio=0.0,
    )

    assert result["probability"] < 1e-8, (
        "Cavity fraction approaching zero should suppress insertions, "
        f"got probability={result['probability']}"
    )


def test_zero_molecules_rejects_deletion():
    """Deletion probability must be zero when the system is empty."""
    acc = _make_acceptance()
    acc.setTemperature(300.0)
    acc.setVolume(50.0)

    result = acc.calculate_deletion_probability_detailed(
        typeId=0,
        currentNumber=0,
        deltaE=-10.0,
        cavityFraction=1.0,
        lambdaNm=1.0,
        rosenbluthWeight=1.0,
        cbmcTrials=1,
        proposalLogRatio=0.0,
    )

    assert result["probability"] == 0.0, "Cannot delete from an empty system"


def test_large_lambda_suppresses_insertion():
    """Increasing thermal wavelength should reduce insertion probability."""
    acc = _make_acceptance()
    acc.setTemperature(300.0)
    acc.setVolume(1.0)

    lambda_small = 0.1
    lambda_large = 10.0
    n_before = 5
    delta_e = 5.0

    result_small = acc.calculate_insertion_probability_detailed(
        typeId=0,
        currentNumber=n_before,
        deltaE=delta_e,
        cavityFraction=1.0,
        lambdaNm=lambda_small,
        rosenbluthWeight=1.0,
        cbmcTrials=1,
        proposalLogRatio=0.0,
    )

    result_large = acc.calculate_insertion_probability_detailed(
        typeId=0,
        currentNumber=n_before,
        deltaE=delta_e,
        cavityFraction=1.0,
        lambdaNm=lambda_large,
        rosenbluthWeight=1.0,
        cbmcTrials=1,
        proposalLogRatio=0.0,
    )

    p_small = result_small["probability"]
    p_large = result_large["probability"]

    assert p_small > 0.0, "Small lambda case should yield a finite probability"
    assert p_large > 0.0, "Large lambda case should still be finite"
    assert p_large < p_small, "Larger lambda should reduce acceptance probability"

    expected_ratio = (lambda_small / lambda_large) ** 3
    actual_ratio = p_large / p_small
    assert actual_ratio == pytest.approx(expected_ratio, rel=0.1), (
        f"Expected probability ratio {expected_ratio:.2e}, "
        f"observed {actual_ratio:.2e}"
    )


def test_zero_temperature_accepts_only_favorable():
    """At very low temperature, only energetically favorable insertions survive."""
    acc = _make_acceptance()
    acc.setTemperature(1.0)
    acc.setVolume(10.0)

    favorable = acc.calculate_insertion_probability_detailed(
        typeId=0,
        currentNumber=5,
        deltaE=-10.0,
        cavityFraction=1.0,
        lambdaNm=1.0,
        rosenbluthWeight=1.0,
        cbmcTrials=1,
        proposalLogRatio=0.0,
    )

    unfavorable = acc.calculate_insertion_probability_detailed(
        typeId=0,
        currentNumber=5,
        deltaE=10.0,
        cavityFraction=1.0,
        lambdaNm=1.0,
        rosenbluthWeight=1.0,
        cbmcTrials=1,
        proposalLogRatio=0.0,
    )

    assert favorable["probability"] > 0.99, "Low temperature should accept exothermic moves"
    assert unfavorable["probability"] < 1e-250, (
        "Low temperature should reject endothermic moves, "
        f"got probability={unfavorable['probability']}"
    )

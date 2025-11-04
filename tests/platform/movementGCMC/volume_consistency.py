"""
Volume scaling tests for GCMCAcceptance probabilities.
"""

import pygcmc


def test_volume_affects_acceptance():
    """Changing the configured volume should scale insertion probability."""
    acc = pygcmc.GCMCAcceptance()
    acc.setTemperature(300.0)
    acc.setChemicalPotential(0, -8.0)
    acc.setThermalLambda(0, 1.0)

    n_before = 8
    delta_e = 4.0
    cavity_fraction = 0.5

    volumes = [10.0, 100.0]
    probabilities = []

    for volume in volumes:
        acc.setVolume(volume)
        result = acc.calculate_insertion_probability_detailed(
            typeId=0,
            currentNumber=n_before,
            deltaE=delta_e,
            cavityFraction=cavity_fraction,
            lambdaNm=1.0,
            rosenbluthWeight=1.0,
            cbmcTrials=1,
            proposalLogRatio=0.0,
        )
        probabilities.append(result["probability"])

    p_small, p_large = probabilities
    assert p_small > 0.0 and p_large > 0.0, "Probabilities should be finite"

    observed_ratio = p_large / p_small
    expected_ratio = volumes[1] / volumes[0]

    rel_error = abs(observed_ratio - expected_ratio) / expected_ratio
    assert rel_error < 0.05, (
        f"Insertion probability should scale with volume. "
        f"Observed ratio={observed_ratio:.3f}, expected={expected_ratio:.1f}"
    )

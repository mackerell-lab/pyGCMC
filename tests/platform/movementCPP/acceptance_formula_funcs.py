"""
Unit tests for the low-level log-space acceptance calculators exposed through
pygcmc.movement.  These checks guard against regressions where the cavity bias
or thermal wavelength factors are accidentally dropped from the probability.
"""

from __future__ import annotations

import math

import pygcmc


def _clamp_probability(log_ratio: float) -> float:
    """Apply the same clamping that the C++ helper uses."""
    if log_ratio >= 0.0:
        return 1.0
    return math.exp(max(log_ratio, -700.0))


def _expected_insertion_probability(
    n: int,
    delta_e: float,
    beta: float,
    chemical_potential: float,
    cavity_bias: float,
    volume_nm3: float,
    lambda_nm: float,
) -> float:
    log_ratio = (
        -beta * delta_e
        + beta * chemical_potential
        + math.log(volume_nm3)
        + math.log(cavity_bias)
        - math.log(n + 1.0)
        - 3.0 * math.log(lambda_nm)
    )
    return _clamp_probability(log_ratio)


def _expected_deletion_probability(
    n: int,
    delta_e: float,
    beta: float,
    chemical_potential: float,
    cavity_bias: float,
    volume_nm3: float,
    lambda_nm: float,
) -> float:
    log_ratio = (
        beta * delta_e
        - beta * chemical_potential
        + math.log(n)
        - math.log(volume_nm3)
        - math.log(cavity_bias)
        + 3.0 * math.log(lambda_nm)
    )
    return _clamp_probability(log_ratio)


def test_insertion_probability_includes_cavity_and_lambda() -> None:
    n = 4
    delta_e = 0.9
    beta = 0.7
    chemical_potential = -2.0
    cavity_bias = 0.3
    volume_nm3 = 40.0
    lambda_nm = 0.8

    expected = _expected_insertion_probability(
        n, delta_e, beta, chemical_potential, cavity_bias, volume_nm3, lambda_nm
    )
    prob_cpp = pygcmc.movement.calculate_insertion_probability_with_lambda(
        n, delta_e, beta, chemical_potential, cavity_bias, volume_nm3, lambda_nm
    )

    assert math.isclose(prob_cpp, expected, rel_tol=1e-12, abs_tol=1e-12)


def test_deletion_probability_includes_cavity_and_lambda() -> None:
    n = 5
    delta_e = -1.1
    beta = 0.7
    chemical_potential = -2.0
    cavity_bias = 0.3
    volume_nm3 = 40.0
    lambda_nm = 0.8

    expected = _expected_deletion_probability(
        n, delta_e, beta, chemical_potential, cavity_bias, volume_nm3, lambda_nm
    )
    prob_cpp = pygcmc.movement.calculate_deletion_probability_with_cavity_and_lambda(
        n, delta_e, beta, chemical_potential, cavity_bias, volume_nm3, lambda_nm
    )

    assert math.isclose(prob_cpp, expected, rel_tol=1e-12, abs_tol=1e-12)

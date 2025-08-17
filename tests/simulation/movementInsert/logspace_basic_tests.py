# tests/simulation/movementInsert/logspace_stability_tests.py
"""
Log-space numerical stability tests for GCMC movements

Tests the log-space implementation for avoiding numerical overflow/underflow
in GCMC acceptance calculations and Rosenbluth weight computations.
"""

import pytest
import math
import random
import numpy as np
import pygcmc
from typing import List, Optional

# Constants
kB = 0.008314463  # kJ/(mol·K)


def logsumexp(xs: List[float]) -> float:
    """
    Compute log(sum(exp(xs))) in a numerically stable way.
    
    Uses the log-sum-exp trick to avoid overflow/underflow.
    """
    if not xs:
        return -math.inf
    
    m = max(xs)
    if m == -math.inf:
        return -math.inf
    
    s = sum(math.exp(x - m) for x in xs)
    return m + math.log(s)


def prob_from_log(log_acc: float) -> float:
    """
    Convert log acceptance to probability, handling saturation naturally.
    
    Properties:
    - log_acc >= 0 → prob = 1 (saturated)
    - log_acc < -745 → prob ≈ 0 (underflow threshold for 64-bit double)
    - Otherwise → prob = exp(log_acc)
    """
    if log_acc >= 0.0:
        return 1.0
    
    # -745 is approximately the underflow threshold for exp()
    if log_acc < -745.0:
        return 0.0
    
    return math.exp(log_acc)


def stable_rosenbluth_log(orient_energies: List[float], beta: float) -> float:
    """
    Calculate log of Rosenbluth weight in log-space.
    
    log W = log(sum_i exp(-beta * E_i))
    """
    if not orient_energies:
        return 0.0  # log(1) = 0 for empty list
    
    return logsumexp([-beta * e for e in orient_energies])


def acc_prob_insertion_logspace(
    fn: float, 
    n: int, 
    B: float, 
    beta: float, 
    delta_E: float, 
    orient_energies: Optional[List[float]] = None
) -> float:
    """
    Calculate insertion acceptance probability in log-space.
    
    Formula: A(n→n+1) = min{1, fn/(n+1) * exp[B - β*ΔE] * W}
    Log form: log_acc = log(fn) - log(n+1) + B - β*ΔE + log(W)
    """
    # Handle edge cases
    if fn <= 0:
        return 0.0
    
    # Calculate log terms
    log_fn = math.log(fn)
    log_n1 = math.log(n + 1.0)
    
    # Configurational bias (Rosenbluth weight)
    logW = stable_rosenbluth_log(orient_energies, beta) if orient_energies else 0.0
    
    # Sum all log terms
    log_acc = log_fn - log_n1 + B - beta * delta_E + logW
    
    return prob_from_log(log_acc)


def acc_prob_deletion_logspace(
    fn: float,
    n: int,
    B: float,
    beta: float,
    delta_E: float,
    orient_energies: Optional[List[float]] = None
) -> float:
    """
    Calculate deletion acceptance probability in log-space.
    
    Formula: A(n→n-1) = min{1, n/fn * exp[-B - β*ΔE] / W}
    Log form: log_acc = log(n) - log(fn) - B - β*ΔE - log(W)
    """
    # Handle edge cases
    if n <= 0 or fn <= 0:
        return 0.0
    
    # Calculate log terms
    log_n = math.log(n)
    log_fn = math.log(fn)
    
    # Configurational bias
    logW = stable_rosenbluth_log(orient_energies, beta) if orient_energies else 0.0
    
    # Sum all log terms
    log_acc = log_n - log_fn - B - beta * delta_E - logW
    
    return prob_from_log(log_acc)


def test_logsumexp_extreme_values():
    """Test logsumexp with extreme values that would overflow in linear space."""
    # Test with very large values
    xs = [1000, 1001, 999, 1002]
    result = logsumexp(xs)
    
    # Should not overflow or return nan/inf
    assert not math.isnan(result)
    assert not math.isinf(result)
    
    # FIX: Result should be >= max (always) and <= max + log(len)
    assert result >= max(xs)
    assert result <= max(xs) + math.log(len(xs)) + 1e-12  # Small tolerance for numerical error
    
    # Test with very negative values
    xs = [-1000, -999, -1001, -998]
    result = logsumexp(xs)
    
    assert not math.isnan(result)
    # Result should be close to max but slightly larger due to other terms
    assert abs(result - max(xs)) < 1.0  # Allow larger tolerance for multiple similar values
    
    # Test with mixed extreme values
    xs = [-1000, 0, 1000]
    result = logsumexp(xs)
    
    assert not math.isnan(result)
    # Should be dominated by largest value
    assert abs(result - 1000) < 0.1


def test_prob_from_log_saturation():
    """Test probability conversion with saturation behavior."""
    # Positive log_acc should saturate to 1
    assert prob_from_log(10.0) == 1.0
    assert prob_from_log(0.0) == 1.0
    
    # Very negative should underflow to 0
    assert prob_from_log(-1000.0) == 0.0
    assert prob_from_log(-746.0) == 0.0
    
    # Normal range
    log_acc = -1.0
    prob = prob_from_log(log_acc)
    assert abs(prob - math.exp(-1.0)) < 1e-10
    
    # Edge of underflow
    log_acc = -744.0
    prob = prob_from_log(log_acc)
    assert 0 < prob < 1e-300  # Very small but not zero


def test_insertion_acceptance_extreme_parameters():
    """Test insertion acceptance with extreme parameters."""
    # Test cases with extreme values
    test_cases = [
        # (fn, n, B, beta, delta_E, orient_energies, description)
        (1e-10, 1, -100.0, 10.0, 1000.0, None, "Tiny cavity, huge energy"),
        (0.5, 1000000, 100.0, 0.001, -1000.0, None, "Many molecules, favorable"),
        (1.0, 0, 1000.0, 1000.0, -10000.0, None, "First insertion, very favorable"),
        (1e-6, 1, -1000.0, 1000.0, 10000.0, None, "Extreme unfavorable"),
    ]
    
    for fn, n, B, beta, delta_E, orient_energies, desc in test_cases:
        acc = acc_prob_insertion_logspace(fn, n, B, beta, delta_E, orient_energies)
        
        # Should always get valid probability
        assert 0 <= acc <= 1, f"{desc}: Invalid acceptance {acc}"
        assert not math.isnan(acc), f"{desc}: NaN acceptance"
        assert not math.isinf(acc), f"{desc}: Inf acceptance"


def test_deletion_acceptance_extreme_parameters():
    """Test deletion acceptance with extreme parameters."""
    test_cases = [
        # (fn, n, B, beta, delta_E, orient_energies, description)
        (0.001, 1000000, 100.0, 10.0, -1000.0, None, "Many molecules, favorable deletion"),
        (1.0, 1, -100.0, 0.001, 1000.0, None, "Last molecule, unfavorable"),
        (0.5, 100, 0.0, 1000.0, 0.0, None, "High temperature, neutral"),
    ]
    
    for fn, n, B, beta, delta_E, orient_energies, desc in test_cases:
        acc = acc_prob_deletion_logspace(fn, n, B, beta, delta_E, orient_energies)
        
        assert 0 <= acc <= 1, f"{desc}: Invalid acceptance {acc}"
        assert not math.isnan(acc), f"{desc}: NaN acceptance"
        assert not math.isinf(acc), f"{desc}: Inf acceptance"


def test_rosenbluth_weight_extreme_energies():
    """Test Rosenbluth weight calculation with extreme energies."""
    beta = 1.0 / (kB * 300.0)  # Room temperature
    
    # Test with very high energies
    energies = [10000, 20000, 30000, 40000, 50000]  # kJ/mol
    logW = stable_rosenbluth_log(energies, beta)
    
    # Should not overflow
    assert not math.isnan(logW)
    assert not math.isinf(logW)
    
    # Should be dominated by lowest energy
    expected = -beta * min(energies)
    assert abs(logW - expected) < 1.0  # Small correction from other terms
    
    # Test with very negative energies
    energies = [-10000, -20000, -30000]
    logW = stable_rosenbluth_log(energies, beta)
    
    assert not math.isnan(logW)
    assert not math.isinf(logW)
    
    # Test with mixed energies
    energies = [-1000, 0, 1000, 10000]
    logW = stable_rosenbluth_log(energies, beta)
    
    assert not math.isnan(logW)
    assert not math.isinf(logW)



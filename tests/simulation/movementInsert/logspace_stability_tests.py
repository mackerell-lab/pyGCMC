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


def test_configurational_bias_logspace():
    """Test configurational bias with log-space calculations."""
    # System parameters
    T = 300.0
    beta = 1.0 / (kB * T)
    fn = 0.1
    n = 10
    B = -2.0
    
    # Generate 100 trial orientations with wide energy range
    random.seed(42)
    orient_energies = [random.gauss(0, 100) for _ in range(100)]
    
    # Calculate with log-space (should be stable)
    delta_E = 10.0
    acc = acc_prob_insertion_logspace(fn, n, B, beta, delta_E, orient_energies)
    
    # Should get valid result
    assert 0 <= acc <= 1
    assert not math.isnan(acc)
    
    # Test that including orientations affects acceptance
    acc_no_orient = acc_prob_insertion_logspace(fn, n, B, beta, delta_E, None)
    
    # They should be different (unless acceptance is saturated)
    if acc < 0.99 and acc_no_orient < 0.99:
        assert acc != acc_no_orient, "Orientations should affect acceptance"


def test_detailed_balance_logspace():
    """Test detailed balance is maintained with log-space calculations."""
    # System parameters
    T = 300.0
    beta = 1.0 / (kB * T)
    fn = 0.3
    mu_ex = -5.0
    n_bar = 100
    B = beta * mu_ex + math.log(n_bar)
    
    # Test multiple state pairs
    test_cases = [
        (10, -50.0, 5.0),   # n, E_n, delta_E_insertion
        (100, -500.0, 10.0),
        (1, -5.0, 2.0),
    ]
    
    for n_A, E_A, delta_E in test_cases:
        n_B = n_A + 1
        E_B = E_A + delta_E
        
        # Calculate acceptance probabilities using log-space functions
        alpha_ins = acc_prob_insertion_logspace(fn, n_A, B, beta, delta_E)
        alpha_del = acc_prob_deletion_logspace(fn, n_B, B, beta, -delta_E)
        
        # Proposal probabilities
        # For GCMC: q_del/q_ins = f_n/(n_B)
        q_ins = 1.0
        q_del = fn / n_B
        
        # Grand canonical weights: π_GC(n,E) ∝ exp(-βE + Bn)
        pi_A = math.exp(-beta * E_A + B * n_A)
        pi_B = math.exp(-beta * E_B + B * n_B)
        
        # Full transition probabilities
        forward = pi_A * q_ins * alpha_ins
        reverse = pi_B * q_del * alpha_del
        
        # Check detailed balance with relative error
        if forward > 0 and reverse > 0:
            rel_err = abs(forward - reverse) / max(forward, reverse)
            assert rel_err < 1e-12, \
                f"Detailed balance violated: forward={forward}, reverse={reverse}, rel_err={rel_err}"


def test_extreme_temperature_logspace():
    """Test log-space calculations at extreme temperatures."""
    fn = 0.5
    n = 10
    mu_ex = -5.0
    n_bar = 100
    delta_E = 10.0
    orient_energies = [0, 10, 20, 30, 40]
    
    # Very low temperature (beta → ∞)
    T_low = 0.1  # K
    beta_low = 1.0 / (kB * T_low)
    B_low = beta_low * mu_ex + math.log(n_bar)
    
    acc_low = acc_prob_insertion_logspace(fn, n, B_low, beta_low, delta_E, orient_energies)
    assert 0 <= acc_low <= 1
    assert not math.isnan(acc_low)
    
    # Very high temperature (beta → 0)
    T_high = 100000.0  # K
    beta_high = 1.0 / (kB * T_high)
    B_high = beta_high * mu_ex + math.log(n_bar)
    
    acc_high = acc_prob_insertion_logspace(fn, n, B_high, beta_high, delta_E, orient_energies)
    assert 0 <= acc_high <= 1
    assert not math.isnan(acc_high)
    
    # At high T, energy differences matter less
    # So acceptance should be closer to geometric factors only
    geometric_acc = fn / (n + 1) * math.exp(B_high)
    if geometric_acc < 1:
        assert abs(acc_high - geometric_acc) < 0.1


def test_parallel_logspace_calculations():
    """Test that log-space calculations work correctly in parallel scenarios."""
    # Simulate parallel calculation for multiple molecules
    n_molecules = 100
    n_orientations = 50
    
    random.seed(123)
    T = 300.0
    beta = 1.0 / (kB * T)
    fn = 0.3
    n = 50
    # FIX: Use parameters that avoid saturation
    beta = 0.05  # Lower beta to avoid saturation
    B = -5.0  # More reasonable B value
    
    results = []
    
    for mol_idx in range(n_molecules):
        # FIX: Generate energies with smaller variance to avoid saturation
        delta_E = random.gauss(0, 5)  # Reduced from 20
        orient_energies = [random.gauss(0, 5) for _ in range(n_orientations)]  # Reduced from 50
        
        # Calculate acceptance
        acc = acc_prob_insertion_logspace(fn, n, B, beta, delta_E, orient_energies)
        results.append(acc)
        
        # Verify each result is valid
        assert 0 <= acc <= 1
        assert not math.isnan(acc)
        assert not math.isinf(acc)
    
    # Check statistics of results
    assert len(results) == n_molecules
    avg_acc = sum(results) / len(results)
    assert 0 <= avg_acc <= 1
    
    # Should have some variation
    min_acc = min(results)
    max_acc = max(results)
    assert max_acc > min_acc, "Should have variation in acceptance probabilities"
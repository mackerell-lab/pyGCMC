# tests/simulation/movementInsert/logspace_advanced_tests.py
"""
Advanced log-space numerical stability tests
"""

import pytest
import math
import random
import numpy as np
import pygcmc

# Constants
kB = 0.008314463  # kJ/(mol·K)

# Import from basic tests
from .logspace_basic_tests import (
    logsumexp,
    prob_from_log,
    stable_rosenbluth_log,
    acc_prob_insertion_logspace,
    acc_prob_deletion_logspace
)

# Import system creation helpers
from .basic_insertion_helpers import (
    create_empty_system,
    create_water_molecule,
    insert_molecule,
    calculate_system_energy
)

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
# tests/simulation/movementInsert/numerical_stability_tests.py
"""
Numerical stability tests for GCMC insertion movements

Tests the numerical stability of insertion operations with extreme values,
including log-space calculations and Rosenbluth weight computations.
"""

import pytest
import math
import random
import numpy as np
import pygcmc
from typing import List, Tuple

# Constants
kB = 0.008314463  # kJ/(mol·K)
kC = 138.935456  # Coulomb constant in kJ·nm/mol/e²


def calculate_boltzmann_weights(energies: List[float], beta: float) -> List[float]:
    """Calculate Boltzmann weights with numerical stability using log-space."""
    if not energies:
        return []
    
    # Use log-space to avoid overflow/underflow
    min_energy = min(energies)
    log_weights = [-beta * (e - min_energy) for e in energies]
    
    # Find maximum for numerical stability
    max_log = max(log_weights)
    
    # Calculate weights
    weights = []
    for log_w in log_weights:
        if log_w - max_log < -700:  # Would underflow
            weights.append(0.0)
        else:
            weights.append(math.exp(log_w - max_log))
    
    # Normalize
    total = sum(weights)
    if total > 0:
        weights = [w / total for w in weights]
    
    return weights


def calculate_rosenbluth_factor(energies: List[float], beta: float) -> float:
    """Calculate Rosenbluth factor W = Σ exp(-βE_i) with numerical stability."""
    if not energies:
        return 1.0
    
    # Shift energies to avoid overflow
    min_energy = min(energies)
    
    # Calculate in log-space
    log_terms = [-beta * (e - min_energy) for e in energies]
    max_log = max(log_terms)
    
    # Sum exp(log_terms - max_log) and restore scale
    sum_exp = sum(math.exp(lt - max_log) for lt in log_terms)
    
    # Final result: exp(max_log) * sum_exp * exp(-beta * min_energy)
    return sum_exp * math.exp(max_log - beta * min_energy)


def test_extreme_energy_numerical_stability():
    """Test numerical stability with extreme energy values."""
    # Test case 1: Very high energies that would overflow
    energies = [0, 100, 200, 1000, 10000]  # kJ/mol
    T = 300.0  # K
    beta = 1.0 / (kB * T)
    
    # Calculate weights without overflow
    weights = calculate_boltzmann_weights(energies, beta)
    
    # Verify weights are normalized
    assert abs(sum(weights) - 1.0) < 1e-10, "Weights should sum to 1"
    
    # Verify highest energy has smallest weight
    assert weights[-1] <= weights[0], "Higher energy should have lower weight"
    
    # Verify no NaN or Inf
    for w in weights:
        assert not math.isnan(w), "Weight should not be NaN"
        assert not math.isinf(w), "Weight should not be Inf"
        assert 0 <= w <= 1, "Weight should be between 0 and 1"
    
    # Calculate Rosenbluth factor
    W = calculate_rosenbluth_factor(energies, beta)
    assert W > 0, "Rosenbluth factor should be positive"
    assert not math.isnan(W), "Rosenbluth factor should not be NaN"
    assert not math.isinf(W), "Rosenbluth factor should not be Inf"


def test_configurational_bias_numerical_stability():
    """Test configurational bias with many trial orientations."""
    # Simulate 100 trial orientations
    n_trials = 100
    random.seed(42)
    
    # Generate random energies with wide range
    energies = [random.uniform(-100, 100) for _ in range(n_trials)]
    
    T = 300.0
    beta = 1.0 / (kB * T)
    
    # Calculate weights
    weights = calculate_boltzmann_weights(energies, beta)
    
    # Select configuration based on weights
    cumsum = 0.0
    r = random.random()
    selected = -1
    
    for i, w in enumerate(weights):
        cumsum += w
        if r < cumsum:
            selected = i
            break
    
    assert selected >= 0, "Should select a configuration"
    assert selected < n_trials, "Selected index should be valid"
    
    # Calculate Rosenbluth factor
    W = calculate_rosenbluth_factor(energies, beta)
    
    # Calculate acceptance probability (simplified)
    # P_acc = min(1, W_new/W_old * exp(-β*ΔE) * other_factors)
    # For insertion: P_acc = min(1, f_n/(n+1) * exp(B) * W_new)
    
    f_n = 0.1  # Cavity fraction
    n = 10  # Current number of molecules
    B = -2.0  # Chemical potential term
    
    # Use log-space for acceptance calculation
    log_acc = math.log(f_n) - math.log(n + 1) + B + math.log(W)
    acc_prob = min(1.0, math.exp(log_acc)) if log_acc < 0 else 1.0
    
    assert 0 <= acc_prob <= 1, "Acceptance probability should be between 0 and 1"


def test_temperature_extremes():
    """Test numerical stability at extreme temperatures."""
    energies = [0, 10, 20, 30, 40]  # kJ/mol
    
    # Test very low temperature (should strongly favor lowest energy)
    T_low = 1.0  # K
    beta_low = 1.0 / (kB * T_low)
    weights_low = calculate_boltzmann_weights(energies, beta_low)
    
    # At low T, lowest energy should dominate
    assert weights_low[0] > 0.99, "Lowest energy should dominate at low T"
    
    # Test very high temperature (should give more uniform distribution)
    T_high = 10000.0  # K
    beta_high = 1.0 / (kB * T_high)
    weights_high = calculate_boltzmann_weights(energies, beta_high)
    
    # At high T, weights should be more uniform
    weight_ratio = max(weights_high) / min(weights_high)
    assert weight_ratio < 10, "Weights should be more uniform at high T"
    
    # Test that no numerical issues occur
    for weights in [weights_low, weights_high]:
        assert abs(sum(weights) - 1.0) < 1e-10
        for w in weights:
            assert not math.isnan(w) and not math.isinf(w)
            assert 0 <= w <= 1


def test_gcmc_insertion_acceptance_stability():
    """Test GCMC insertion acceptance calculation stability."""
    # Create test system
    system = pygcmc.MCState()
    system.info.box = [5.0, 5.0, 5.0]
    system.info.setTemperature(300.0)
    system.info.cutoff = 1.2
    
    # Create forcefield
    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 2
    ff.ljSigma = [0.3, 0.35, 0.35, 0.4]
    ff.ljEps = [0.5, 0.4, 0.4, 0.6]
    system.forcefield = ff
    
    # Test acceptance calculation with various parameters
    test_cases = [
        # (f_n, n, B, delta_E, description)
        (0.1, 0, 0.0, 0.0, "First insertion, neutral"),
        (0.01, 100, -5.0, 10.0, "Many molecules, unfavorable"),
        (0.5, 10, 2.0, -10.0, "Favorable insertion"),
        (1e-6, 1000, -10.0, 100.0, "Extreme: tiny cavity, many molecules"),
        (1.0, 0, 10.0, -100.0, "Extreme: very favorable"),
    ]
    
    T = 300.0
    beta = 1.0 / (kB * T)
    
    for f_n, n, B, delta_E, desc in test_cases:
        # Calculate acceptance probability in log-space
        if f_n > 0:
            log_acc = math.log(f_n) - math.log(n + 1) + B - beta * delta_E
            acc_prob = min(1.0, math.exp(log_acc)) if log_acc < 0 else 1.0
        else:
            acc_prob = 0.0
        
        # Verify result is valid
        assert 0 <= acc_prob <= 1, f"{desc}: Invalid acceptance probability {acc_prob}"
        assert not math.isnan(acc_prob), f"{desc}: NaN acceptance probability"
        assert not math.isinf(acc_prob), f"{desc}: Inf acceptance probability"
        
        # Additional sanity checks
        if f_n == 0:
            assert acc_prob == 0, "Zero cavity should give zero acceptance"
        if n == 0 and delta_E < 0 and B > 0:
            assert acc_prob > 0.5, "Very favorable first insertion should have high acceptance"


def test_parallel_rosenbluth_calculation():
    """Test parallel calculation of Rosenbluth weights for multiple molecules."""
    # Simulate calculating weights for multiple molecules in parallel
    n_molecules = 10
    n_orientations = 20
    
    random.seed(123)
    T = 300.0
    beta = 1.0 / (kB * T)
    
    # Generate energies for each molecule's orientations
    all_weights = []
    all_rosenbluth = []
    
    for mol_idx in range(n_molecules):
        energies = [random.gauss(0, 50) for _ in range(n_orientations)]
        
        # Calculate weights
        weights = calculate_boltzmann_weights(energies, beta)
        all_weights.append(weights)
        
        # Calculate Rosenbluth factor
        W = calculate_rosenbluth_factor(energies, beta)
        all_rosenbluth.append(W)
        
        # Verify correctness
        assert abs(sum(weights) - 1.0) < 1e-10
        assert W > 0 and not math.isnan(W) and not math.isinf(W)
    
    # Verify all calculations succeeded
    assert len(all_weights) == n_molecules
    assert len(all_rosenbluth) == n_molecules
    
    # Check for reasonable variation in Rosenbluth factors
    W_min = min(all_rosenbluth)
    W_max = max(all_rosenbluth)
    assert W_max / W_min > 1.0, "Should have some variation in Rosenbluth factors"


def test_insertion_with_pbc_wrapping():
    """Test numerical stability when coordinates wrap around periodic boundaries."""
    box_size = 5.0
    
    # Test coordinates that need wrapping
    test_coords = [
        (6.0, 3.0, 2.0),   # x > box
        (-1.0, 3.0, 2.0),  # x < 0
        (2.0, 10.0, 2.0),  # y > 2*box
        (2.0, 3.0, -5.0),  # z < -box
    ]
    
    for x, y, z in test_coords:
        # Apply PBC
        x_wrapped = x - box_size * round(x / box_size)
        y_wrapped = y - box_size * round(y / box_size)
        z_wrapped = z - box_size * round(z / box_size)
        
        # Verify wrapped coordinates are in box
        assert -box_size/2 <= x_wrapped <= box_size/2
        assert -box_size/2 <= y_wrapped <= box_size/2
        assert -box_size/2 <= z_wrapped <= box_size/2
        
        # Verify no numerical issues
        assert not math.isnan(x_wrapped)
        assert not math.isnan(y_wrapped)
        assert not math.isnan(z_wrapped)
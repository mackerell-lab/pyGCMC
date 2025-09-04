# tests/simulation/movementCPP/statistical_utils_tests_basic_funcs.py
"""Basic tests for statistical utility functions."""

import pytest
import numpy as np
from .statistical_utils import (
    poisson_diff_ok,
    ratio_CI_ok,
    effective_sample_size,
    batch_means_variance,
    no_drift,
    calibrate_mu
)

def test_poisson_diff_ok_edge_cases():
    """Test Poisson difference with edge cases."""
    # Zero counts - should return False
    assert not poisson_diff_ok(0, 0)
    assert not poisson_diff_ok(0, 10)
    assert not poisson_diff_ok(10, 0)
    
    # Equal counts - should pass
    assert poisson_diff_ok(10, 10)
    assert poisson_diff_ok(100, 100)
    
    # Small differences within Poisson noise
    assert poisson_diff_ok(100, 110)  # ~10% difference, within sqrt(n) noise
    assert poisson_diff_ok(95, 105)
    
    # Large differences - should fail
    assert not poisson_diff_ok(10, 100)
    assert not poisson_diff_ok(100, 10)
    assert not poisson_diff_ok(50, 200)
    

def test_poisson_diff_ok_statistical_power():
    """Test statistical power of Poisson test."""
    # Generate counts from Poisson distributions
    np.random.seed(42)
    
    # Same rate - should mostly pass
    rate = 50
    passes = 0
    for _ in range(100):
        k1 = np.random.poisson(rate)
        k2 = np.random.poisson(rate)
        if poisson_diff_ok(k1, k2):
            passes += 1
    
    assert passes > 90, f"Only {passes}/100 passed for same rate"
    
    # Different rates - should mostly fail
    rate1, rate2 = 50, 100
    fails = 0
    for _ in range(100):
        k1 = np.random.poisson(rate1)
        k2 = np.random.poisson(rate2)
        if not poisson_diff_ok(k1, k2):
            fails += 1
    
    assert fails > 90, f"Only {fails}/100 failed for different rates"
    

def test_ratio_ci_ok_edge_cases():
    """Test ratio confidence interval with edge cases."""
    # Zero counts
    assert not ratio_CI_ok(0, 0)
    assert not ratio_CI_ok(0, 10)
    assert not ratio_CI_ok(10, 0)
    
    # Equal counts - ratio = 1
    assert ratio_CI_ok(10, 10)
    assert ratio_CI_ok(100, 100)
    
    # Ratios near 1 (within typical CI)
    assert ratio_CI_ok(90, 100)
    assert ratio_CI_ok(100, 90)
    
    # Extreme ratios
    assert not ratio_CI_ok(10, 100)  # ratio = 0.1
    assert not ratio_CI_ok(100, 10)  # ratio = 10
    

def test_ratio_ci_ok_confidence_level():
    """Test confidence level of ratio test."""
    np.random.seed(123)
    
    # Generate paired counts with known ratio
    n_trials = 1000
    true_ratio = 1.0
    
    # Counts with ratio = 1 (independent Poisson)
    passes = 0
    for _ in range(n_trials):
        k1 = np.random.poisson(50)
        k2 = np.random.poisson(50)
        if ratio_CI_ok(k1, k2):
            passes += 1
    
    # Should pass ~95% of the time (95% CI) with wider tolerance
    assert 0.92 < passes/n_trials < 0.98, \
        f"Pass rate {passes/n_trials:.2%} outside expected 95% CI"
    

def test_effective_sample_size_independent():
    """Test ESS for independent samples."""
    np.random.seed(456)
    
    # Independent samples - ESS should be close to n
    n = 1000
    independent = np.random.normal(0, 1, n)
    ess = effective_sample_size(independent)
    
    # ESS should be at least 90% of n for independent data
    assert ess > 0.9 * n, f"ESS {ess:.0f} too low for independent data"
    assert ess <= n, f"ESS {ess:.0f} cannot exceed n={n}"
    

def test_effective_sample_size_correlated():
    """Test ESS for autocorrelated samples."""
    np.random.seed(789)
    
    # Generate autocorrelated data (AR(1) process)
    n = 1000
    rho = 0.8  # autocorrelation
    data = [0.0]
    for _ in range(n-1):
        data.append(rho * data[-1] + np.random.normal(0, 1))
    
    data = np.array(data)
    ess = effective_sample_size(data)
    
    # ESS should be much less than n for correlated data
    # Theoretical ESS for AR(1) ≈ n * (1-rho)/(1+rho)
    expected_ess = n * (1 - rho) / (1 + rho)
    
    assert ess < 0.5 * n, f"ESS {ess:.0f} too high for correlated data"
    assert ess > 0.5 * expected_ess, f"ESS {ess:.0f} too low vs theory {expected_ess:.0f}"
    

def test_effective_sample_size_edge_cases():
    """Test ESS edge cases."""
    # Constant data - no variance
    constant = np.ones(100)
    ess = effective_sample_size(constant)
    assert ess == 1.0, "ESS for constant data should be 1"
    
    # Very short series
    short = [1, 2, 3]
    ess = effective_sample_size(short)
    assert 1 <= ess <= 3, f"ESS {ess} out of bounds for n=3"
    
    # Single value
    single = [42]
    ess = effective_sample_size(single)
    assert ess == 1.0, "ESS for single value should be 1"
    

def test_batch_means_edge_cases():
    """Test batch means variance edge cases."""
    # Very short data
    short = [1, 2, 3]
    var = batch_means_variance(short, batch_size=2)
    assert var > 0, "Should return positive variance"
    
    # Constant data
    constant = np.ones(100)
    var = batch_means_variance(constant, batch_size=10)
    assert var == 0, "Should return zero variance for constant data"
    
    # Single batch
    data = [1, 2, 3, 4, 5]
    var = batch_means_variance(data, batch_size=5)
    assert var == 0, "Single batch should give zero variance"
    
# tests/simulation/movementCPP/statistical_utils_tests_advanced_funcs.py
"""Advanced tests for statistical utility functions."""

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

def test_batch_means_variance():
    """Test batch means variance estimation."""
    np.random.seed(111)
    
    # Generate data with known variance
    n = 10000
    true_var = 4.0
    data = np.random.normal(0, np.sqrt(true_var), n)
    
    # Estimate variance using different batch sizes
    for batch_size in [10, 50, 100]:
        est_var = batch_means_variance(data, batch_size)
        
        # Should be close to true variance
        assert est_var > 0, "Variance must be positive"
        assert 0.5 * true_var < est_var < 2.0 * true_var, \
            f"Estimated variance {est_var:.2f} far from true {true_var}"
    

def test_batch_means_variance_autocorrelated():
    """Test batch means with autocorrelated data."""
    np.random.seed(222)
    
    # AR(1) process
    n = 10000
    rho = 0.7
    data = [0.0]
    for _ in range(n-1):
        data.append(rho * data[-1] + np.random.normal(0, 1))
    
    data = np.array(data)
    
    # Small batches underestimate variance
    var_small = batch_means_variance(data, batch_size=5)
    
    # Large batches give better estimate
    var_large = batch_means_variance(data, batch_size=100)
    
    # Large batches should give higher (more accurate) variance estimate
    assert var_large > var_small, \
        "Large batches should give higher variance for correlated data"
    

def test_no_drift_converged():
    """Test drift detection for stable data."""
    np.random.seed(333)
    
    # Generate stable data
    n = 1000
    target_mean = 10.0
    data = np.random.normal(target_mean, 1.0, n)
    
    # Should detect no drift
    assert no_drift(data), "Should detect no drift for stable data"
    

def test_no_drift_drifting():
    """Test drift detection for drifting data."""
    np.random.seed(444)
    
    # Generate drifting data
    n = 1000
    target_mean = 10.0
    data = np.random.normal(target_mean, 1.0, n)
    
    # Add significant drift
    data += np.linspace(0, 20, n)
    
    assert not no_drift(data), "Should detect drift for drifting data"
    

def test_calibrate_mu():
    """Test chemical potential calibration."""
    # Test calibration
    kT = 2.5  # thermal energy
    mu_initial = -10.0
    n_target = 100
    n_obs = 50  # Too few particles
    
    mu_new = calibrate_mu(mu_initial, n_target, n_obs, kT)
    
    # Should increase mu to get more particles
    assert mu_new > mu_initial, "Should increase mu when n_obs < n_target"
    
    # Test opposite case
    n_obs = 200  # Too many particles
    mu_new = calibrate_mu(mu_initial, n_target, n_obs, kT)
    
    # Should decrease mu to get fewer particles
    assert mu_new < mu_initial, "Should decrease mu when n_obs > n_target"
    

def test_no_drift_window_size():
    """Test drift detection with different window sizes."""
    np.random.seed(555)
    
    # Generate data with slow drift
    n = 2000
    data = np.random.normal(10.0, 1.0, n)
    data += np.linspace(0, 5, n)  # Slow drift
    
    # Small window might miss slow drift
    assert no_drift(data, window=50), \
        "Small window might miss slow drift"
    
    # Large window should detect drift
    assert not no_drift(data, window=500), \
        "Large window should detect slow drift"
    

def test_integration_detailed_balance():
    """Integration test using multiple utilities for detailed balance."""
    np.random.seed(777)
    
    # Simulate insertion/deletion counts
    n_cycles = 1000
    insertion_rate = 50
    deletion_rate = 48  # Slightly different but within noise
    
    insertions = np.random.poisson(insertion_rate)
    deletions = np.random.poisson(deletion_rate)
    
    # All tests should agree for balanced system
    assert poisson_diff_ok(insertions, deletions), \
        "Poisson test should pass for balanced system"
    assert ratio_CI_ok(insertions, deletions), \
        "Ratio CI should pass for balanced system"
    
    # Very different rates
    insertions_biased = np.random.poisson(100)
    deletions_biased = np.random.poisson(30)
    
    assert not poisson_diff_ok(insertions_biased, deletions_biased), \
        "Poisson test should fail for biased system"
    assert not ratio_CI_ok(insertions_biased, deletions_biased), \
        "Ratio CI should fail for biased system"
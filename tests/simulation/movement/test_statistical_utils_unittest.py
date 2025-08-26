# tests/simulation/movement/test_statistical_utils_unittest.py
"""Unit tests for statistical utilities."""

import pytest
import numpy as np
from .test_statistical_utils import (
    poisson_diff_ok,
    ratio_CI_ok,
    effective_sample_size,  # Correct function name
    batch_means_variance,
    no_drift,  # Available function
    calibrate_mu  # Available function
)


class TestStatisticalUtilities:
    """Test statistical utility functions."""
    
    def test_poisson_diff_ok_edge_cases(self):
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
    
    def test_poisson_diff_ok_statistical_power(self):
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
    
    def test_ratio_ci_ok_edge_cases(self):
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
    
    def test_ratio_ci_ok_confidence_level(self):
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
    
    def test_effective_sample_size_independent(self):
        """Test ESS for independent samples."""
        np.random.seed(456)
        
        # Independent samples - ESS should be close to n
        n = 1000
        independent = np.random.normal(0, 1, n)
        ess = effective_sample_size(independent)
        
        # ESS should be at least 90% of n for independent data
        assert ess > 0.9 * n, f"ESS {ess:.0f} too low for independent data"
        assert ess <= n, f"ESS {ess:.0f} cannot exceed n={n}"
    
    def test_effective_sample_size_correlated(self):
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
    
    def test_effective_sample_size_edge_cases(self):
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
    
    def test_batch_means_variance(self):
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
    
    def test_batch_means_variance_autocorrelated(self):
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
    
    def test_no_drift_converged(self):
        """Test drift detection for stable data."""
        np.random.seed(333)
        
        # Generate stable data
        n = 1000
        target_mean = 10.0
        data = np.random.normal(target_mean, 1.0, n)
        
        # Should detect no drift
        assert no_drift(data), "Should detect no drift for stable data"
    
    def test_no_drift_drifting(self):
        """Test drift detection for drifting data."""
        np.random.seed(444)
        
        # Generate drifting data
        n = 1000
        target_mean = 10.0
        data = np.random.normal(target_mean, 1.0, n)
        
        # Add significant drift
        data += np.linspace(0, 20, n)
        
        assert not no_drift(data), "Should detect drift for drifting data"
    
    def test_calibrate_mu(self):
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
    
    def test_no_drift_window_size(self):
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
    
    def test_batch_means_edge_cases(self):
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
    
    def test_integration_detailed_balance(self):
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
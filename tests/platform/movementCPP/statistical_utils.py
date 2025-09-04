# tests/simulation/movement/test_statistical_utils.py
"""Statistical utilities for robust GCMC testing."""

import math
import numpy as np


def poisson_diff_ok(k_i, k_d, z=1.96):
    """Test if insertion/deletion counts are statistically balanced.
    
    Uses Poisson difference test: z = (k_i - k_d) / sqrt(k_i + k_d)
    Returns True if |z| <= 1.96 (95% confidence)
    """
    if k_i + k_d == 0:
        return False
    return abs(k_i - k_d) <= z * math.sqrt(k_i + k_d)


def ratio_CI_ok(k_i, k_d, z=1.96):
    """Test if insertion/deletion ratio is within confidence interval.
    
    Uses log-ratio standard error: SE ≈ sqrt(1/k_i + 1/k_d)
    Returns True if |log(R)| <= z*SE
    """
    if k_i == 0 or k_d == 0:
        return False
    se = math.sqrt(1.0/k_i + 1.0/k_d)
    return abs(math.log(k_i/k_d)) <= z * se


def no_drift(measurements, window=200, z=1.96):
    """Test if measurements show no significant drift.
    
    Compares two consecutive windows using batch means.
    Returns True if difference is within z*SE
    """
    if len(measurements) < 2*window:
        return False
    a = measurements[-2*window:-window]
    b = measurements[-window:]
    mean_a, mean_b = np.mean(a), np.mean(b)
    var_a = np.var(a, ddof=1)
    var_b = np.var(b, ddof=1)
    se = math.sqrt(var_a/len(a) + var_b/len(b))
    if se == 0:  # No variance - system stuck
        return abs(mean_a - mean_b) < 0.1
    return abs(mean_a - mean_b) <= z * se


def run_until_stable(step_fn, max_steps=20000, check_every=500, window=200):
    """Run simulation until stable or max steps reached.
    
    Checks for stability using no_drift every check_every steps.
    Returns list of measurements.
    """
    vals = []
    for t in range(1, max_steps + 1):
        vals.append(step_fn())
        if t % check_every == 0 and t >= 2*window:
            if no_drift(vals, window=window):
                break
    return vals


def calibrate_mu(mu, n_target, n_obs, kT):
    """Calibrate chemical potential to achieve target density.
    
    Uses single-step correction: μ_new = μ_old + kT * ln(n_target / n_obs)
    """
    return mu + kT * math.log(max(n_target, 1e-12) / max(n_obs, 1e-12))


def effective_sample_size(x):
    """Calculate effective sample size accounting for autocorrelation.
    
    ESS = n / (1 + 2*sum(autocorrelations))
    """
    n = len(x)
    if n < 10:
        return n
    
    x = np.array(x, dtype=float) - np.mean(x)
    c0 = np.dot(x, x) / n
    if c0 == 0:
        return 1  # Constant data has ESS = 1
    
    # Sum positive autocorrelations; negative correlations increase ESS
    sum_rho = 0.0
    for k in range(1, min(n//4, 100)):
        ck = np.dot(x[:-k], x[k:]) / (n - k)
        rho_k = ck / c0
        if abs(rho_k) < 0.05:
            break
        if rho_k > 0:
            sum_rho += rho_k
    
    ess = n / (1.0 + 2.0 * sum_rho)
    return max(1, int(ess))


def batch_means_variance(x, batch_size=None):
    """Calculate variance using batch means to account for correlation."""
    n = len(x)
    if batch_size is None:
        batch_size = max(1, int(math.sqrt(n)))
    
    n_batches = n // batch_size
    if n_batches < 2:
        # Exactly one full batch -> undefined for batch-means; return 0
        return 0.0 if n == batch_size else float(np.var(x, ddof=1))
    
    batch_means = []
    for i in range(n_batches):
        batch = x[i*batch_size:(i+1)*batch_size]
        batch_means.append(float(np.mean(batch)))
    
    # Scale by batch_size to estimate per-sample variance
    return float(np.var(batch_means, ddof=1) * batch_size)


def bootstrap_confidence_interval(data, statistic_func=np.mean, n_bootstrap=1000, confidence=0.95):
    """Calculate bootstrap confidence interval for a statistic.
    
    Args:
        data: Input data array
        statistic_func: Function to compute statistic (default: mean)
        n_bootstrap: Number of bootstrap samples
        confidence: Confidence level (default: 0.95)
    
    Returns:
        (lower, upper): Confidence interval bounds
    """
    n = len(data)
    if n < 2:
        stat = statistic_func(data)
        return (stat, stat)
    
    # Generate bootstrap samples
    bootstrap_stats = []
    rng = np.random.RandomState(42)  # Fixed seed for reproducibility
    for _ in range(n_bootstrap):
        # Resample with replacement
        indices = rng.choice(n, size=n, replace=True)
        sample = [data[i] for i in indices]
        stat = statistic_func(sample)
        bootstrap_stats.append(stat)
    
    # Calculate percentiles
    alpha = 1 - confidence
    lower_percentile = (alpha/2) * 100
    upper_percentile = (1 - alpha/2) * 100
    
    lower = np.percentile(bootstrap_stats, lower_percentile)
    upper = np.percentile(bootstrap_stats, upper_percentile)
    
    return (lower, upper)
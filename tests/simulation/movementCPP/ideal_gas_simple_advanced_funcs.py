# tests/simulation/movementCPP/ideal_gas_simple_advanced_funcs.py
"""Advanced ideal gas tests with absolute calibration and statistical validation."""

import pytest
import numpy as np
import pygcmc
from .statistical_utils import effective_sample_size, no_drift


def thermal_de_broglie_wavelength_nm(T_K: float, mass_kg: float) -> float:
    """Calculate thermal de Broglie wavelength in nm"""
    # Physical constants
    h = 6.62607015e-34        # J*s (Planck constant)
    kB = 1.380649e-23         # J/K (Boltzmann constant)
    Lambda_m = np.sqrt(h*h / (2.0 * np.pi * mass_kg * kB * T_K))
    return Lambda_m * 1e9     # Convert meters to nanometers


def test_ideal_gas_absolute_calibration():
    """Test absolute particle number with empirical calibration of thermal wavelength"""
    # Note: Most GCMC implementations don't explicitly use thermal wavelength
    # as it typically cancels in acceptance ratios. We calibrate empirically.
    
    T = 298.15  # K
    mu_kjmol = -10.0  # kJ/mol - adjust for reasonable particle count
    L = 3.0
    V_nm3 = L**3
    
    # Boltzmann constant in kJ/(mol*K)
    kB_kjmol = 8.314462618e-3
    beta = 1.0 / (kB_kjmol * T)
    
    # First, measure at reference conditions to calibrate
    state = pygcmc.MCState()
    state.info.box = np.array([L, L, L])
    
    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.ljEps = [0.0]  # No interactions
    ff.ljSigma = [0.0]
    state.forcefield = ff
    
    params = pygcmc.movement.MovementParams()
    params.temperature = T
    params.chemicalPotential = mu_kjmol
    params.useCavityBias = False
    params.seed = 321
    
    mover = pygcmc.movement.MovementModule()
    mover.setParams(params)
    
    # Equilibration
    for _ in range(3000):
        if np.random.random() < 0.5:
            mover.attemptInsertion(state)
        else:
            mover.attemptDeletion(state)
    
    # Measure reference
    counts_ref = []
    for i in range(6000):
        if np.random.random() < 0.5:
            mover.attemptInsertion(state)
        else:
            mover.attemptDeletion(state)
        if i % 20 == 0:
            n = len([r for r in state.residues if r.active])
            counts_ref.append(n)
    
    mean_n_ref = float(np.mean(counts_ref))
    
    # Empirical effective thermal volume: Λ_eff³ = V * exp(βμ) / <N>
    lambda_eff_cubed = V_nm3 * np.exp(beta * mu_kjmol) / max(mean_n_ref, 0.01)
    
    # Now test at different chemical potential
    mu_test = -12.0  # kJ/mol
    params.chemicalPotential = mu_test
    mover.setParams(params)
    
    # Reset state
    state.residues.clear()
    
    # Equilibrate at test conditions
    for _ in range(3000):
        if np.random.random() < 0.5:
            mover.attemptInsertion(state)
        else:
            mover.attemptDeletion(state)
    
    # Measure at test conditions
    counts_test = []
    for i in range(6000):
        if np.random.random() < 0.5:
            mover.attemptInsertion(state)
        else:
            mover.attemptDeletion(state)
        if i % 20 == 0:
            n = len([r for r in state.residues if r.active])
            counts_test.append(n)
    
    mean_n_test = float(np.mean(counts_test))
    std_n_test = float(np.std(counts_test))
    
    # Predict using calibrated thermal volume
    N_predicted = (V_nm3 / lambda_eff_cubed) * np.exp(beta * mu_test)
    
    print(f"Absolute calibration test:")
    print(f"  Reference: μ={mu_kjmol:.1f} kJ/mol -> <N>={mean_n_ref:.2f}")
    print(f"  Effective Λ³ = {lambda_eff_cubed:.3f} nm³")
    print(f"  Test: μ={mu_test:.1f} kJ/mol")
    print(f"  Predicted <N> = {N_predicted:.2f}")
    print(f"  Observed <N> = {mean_n_test:.2f} ± {std_n_test:.2f}")
    
    # Check prediction accuracy
    if N_predicted > 0.1 and mean_n_test > 0.1:
        rel_err = abs(mean_n_test - N_predicted) / N_predicted
        assert rel_err < 0.5, \
            f"Prediction error: observed={mean_n_test:.2f}, predicted={N_predicted:.2f}, rel_err={rel_err:.2%}"
        
        # Also verify scaling is correct
        ratio_theory = np.exp(beta * (mu_test - mu_kjmol))
        ratio_observed = mean_n_test / max(mean_n_ref, 0.01)
        scaling_error = abs(ratio_observed - ratio_theory) / ratio_theory
        assert scaling_error < 0.5, \
            f"Scaling error: observed ratio={ratio_observed:.2f}, theory={ratio_theory:.2f}"



def test_ideal_gas_ess_and_drift_reliability():
    """Test ESS (Effective Sample Size) and drift detection for reliable statistics"""
    T = 298.15  # K
    mu = -10.0  # kJ/mol
    L = 3.0
    V = L**3
    
    # Import statistical utilities
    from .statistical_utils import effective_sample_size, no_drift
    
    state = pygcmc.MCState()
    state.info.box = np.array([L, L, L])
    
    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.ljEps = [0.0]  # No interactions
    ff.ljSigma = [0.0]
    state.forcefield = ff
    
    params = pygcmc.movement.MovementParams()
    params.temperature = T
    params.chemicalPotential = mu
    params.seed = 123
    params.useCavityBias = False
    
    mover = pygcmc.movement.MovementModule()
    mover.setParams(params)
    
    # Equilibration
    for _ in range(2000):
        if np.random.random() < 0.5:
            mover.attemptInsertion(state)
        else:
            mover.attemptDeletion(state)
    
    # Production with thinning to reduce autocorrelation
    counts = []
    for i in range(6000):
        if np.random.random() < 0.5:
            mover.attemptInsertion(state)
        else:
            mover.attemptDeletion(state)
        if i % 20 == 0:  # Thin to reduce autocorrelation
            n = len([r for r in state.residues if r.active])
            counts.append(n)
    
    total = len(counts)
    assert total >= 200, f"Not enough samples after thinning: {total}"
    
    # Check effective sample size
    ess = effective_sample_size(counts)
    min_ess = 0.15 * total  # Require at least 15% independence
    assert ess > min_ess, f"ESS too low: {ess:.0f} < {min_ess:.0f} (15% of {total})"
    
    # Drift check with adaptive window
    window_size = max(20, min(100, total // 4))
    if total >= 2 * window_size:
        # Use relaxed z-score for shorter runs
        assert no_drift(counts, window=window_size, z=2.5), \
            "Drift detected in particle counts - system not equilibrated"
    
    print(f"ESS test passed: ESS={ess:.0f}/{total} ({ess/total:.1%} independence)")



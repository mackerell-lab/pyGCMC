# tests/platform/movementCPP/ideal_gas_theory_helpers.py
"""
Helper functions for theory-based GCMC validation.
These functions provide exact theoretical predictions for ideal gas systems.
"""

import numpy as np
import math


def theoretical_insertion_deletion_ratio(n_particles, volume, chemical_potential, temperature):
    """
    Calculate the theoretical insertion/deletion probability ratio for ideal gas.
    
    For an ideal gas in the grand canonical ensemble:
    p_ins / p_del = exp(β*μ) * V / (N + 1) * Λ^(-3)
    
    Where Λ is absorbed into the chemical potential definition.
    
    Args:
        n_particles: Current number of particles
        volume: System volume (nm^3)
        chemical_potential: Chemical potential (kJ/mol)
        temperature: Temperature (K)
    
    Returns:
        Theoretical ratio of insertion to deletion acceptance probability
    """
    kT = 8.314e-3 * temperature  # kJ/mol
    beta = 1.0 / kT
    
    # For ideal gas: p_ins/p_del = exp(β*μ) * V / (N+1)
    # Note: Thermal de Broglie wavelength Λ is included in μ
    if n_particles >= 0:
        ratio = math.exp(beta * chemical_potential) * volume / (n_particles + 1)
    else:
        ratio = 0.0
    
    return ratio


def validate_detailed_balance_ratio(insertion_probs, deletion_probs, n_states, 
                                   volume, chemical_potential, temperature,
                                   tolerance=0.3):
    """
    Validate that observed insertion/deletion probabilities match theory.
    
    Args:
        insertion_probs: List of (n, prob) tuples for insertion attempts
        deletion_probs: List of (n, prob) tuples for deletion attempts
        n_states: List of particle numbers for each attempt
        volume: System volume (nm^3)
        chemical_potential: Chemical potential (kJ/mol)
        temperature: Temperature (K)
        tolerance: Relative tolerance for ratio comparison
    
    Returns:
        (passed, observed_ratio, expected_ratio, relative_error)
    """
    if not insertion_probs or not deletion_probs:
        return False, 0, 0, float('inf')
    
    # Group by particle number
    ins_by_n = {}
    del_by_n = {}
    
    for i, (n, prob) in enumerate(insertion_probs):
        if n not in ins_by_n:
            ins_by_n[n] = []
        ins_by_n[n].append(prob)
    
    for i, (n, prob) in enumerate(deletion_probs):
        if n not in del_by_n:
            del_by_n[n] = []
        del_by_n[n].append(prob)
    
    # Calculate average ratio for overlapping N values
    ratios_obs = []
    ratios_exp = []
    
    for n in ins_by_n:
        if n in del_by_n and len(ins_by_n[n]) > 5 and len(del_by_n[n]) > 5:
            mean_ins = np.mean(ins_by_n[n])
            mean_del = np.mean(del_by_n[n])
            
            if mean_del > 0:
                obs_ratio = mean_ins / mean_del
                exp_ratio = theoretical_insertion_deletion_ratio(
                    n, volume, chemical_potential, temperature
                )
                
                ratios_obs.append(obs_ratio)
                ratios_exp.append(exp_ratio)
    
    if not ratios_obs:
        return False, 0, 0, float('inf')
    
    # Compare average ratios
    mean_obs = np.mean(ratios_obs)
    mean_exp = np.mean(ratios_exp)
    
    if mean_exp > 0:
        relative_error = abs(mean_obs - mean_exp) / mean_exp
        passed = relative_error < tolerance
    else:
        relative_error = float('inf')
        passed = False
    
    return passed, mean_obs, mean_exp, relative_error


def ideal_gas_mean_n(volume, chemical_potential, temperature):
    """
    Calculate expected mean particle number for ideal gas.
    
    <N> = exp(β*μ) * V / Λ^3
    
    Where Λ is the thermal de Broglie wavelength (absorbed into μ).
    """
    kT = 8.314e-3 * temperature  # kJ/mol
    beta = 1.0 / kT
    
    # For ideal gas
    mean_n = math.exp(beta * chemical_potential) * volume
    
    return mean_n


def cavity_bias_factor_check(cavity_factor_observed, n_particles, volume, probe_radius=0.15):
    """
    Check if observed cavity bias factor is reasonable.
    
    For low density: cavity_factor ≈ V_cavity / V_total
    For high density: cavity_factor << 1
    
    Returns:
        (is_reasonable, expected_range_min, expected_range_max)
    """
    # Estimate occupied volume
    particle_volume = n_particles * (4.0/3.0) * math.pi * probe_radius**3
    free_fraction = max(0, 1.0 - particle_volume / volume)
    
    # Expected cavity factor should be between 0 and free_fraction
    # Allow some margin for fluctuations
    expected_min = 0.0
    expected_max = min(1.0, free_fraction * 1.5)
    
    is_reasonable = expected_min <= cavity_factor_observed <= expected_max
    
    return is_reasonable, expected_min, expected_max
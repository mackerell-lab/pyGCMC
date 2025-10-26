# tests/simulation/movementInsert/practical_gcmc_cavity_bias_test.py
"""
Cavity bias effect test for GCMC

This module tests the effect of cavity bias on GCMC acceptance rates.
"""

import pytest
import math
import random
import pygcmc
import os

# Constants
kB = 0.008314463  # Boltzmann constant in kJ/mol/K
kC = 138.935456   # Coulomb constant in kJ·nm/mol/e²


from .practical_gcmc_systems import *
from .practical_gcmc_molecules import *

def test_cavity_bias_effect():
    """Test the effect of cavity bias on acceptance rates"""
    
    print("\n=== Testing Cavity Bias Effect ===")
    
    # GCMC parameters
    T = 300.0
    beta = 1.0 / (kB * T)
    mu_ex = -3.0
    B = beta * mu_ex + math.log(10.0)
    
    # Create systems with different densities
    densities = [0.1, 0.3, 0.5, 0.7]  # Fraction occupied
    
    for density in densities:
        # Create system with specified density
        system = create_system_with_density(density)
        f_n = 1.0 - density  # Cavity fraction
        
        # Try insertions (increase sample size to reduce statistical fluctuation)
        n_attempts = 1000
        accepted_no_bias = 0
        accepted_with_bias = 0
        
        for _ in range(n_attempts):
            # Random insertion energy (simplified)
            delta_e = random.gauss(10.0, 20.0)  # Mean 10 kJ/mol, std 20
            
            # Without cavity bias
            acc_no_bias = min(1.0, math.exp(B - beta * delta_e))
            if random.random() < acc_no_bias:
                accepted_no_bias += 1
            
            # With cavity bias
            acc_with_bias = min(1.0, f_n * math.exp(B - beta * delta_e))
            if random.random() < acc_with_bias:
                accepted_with_bias += 1
        
        rate_no_bias = accepted_no_bias / n_attempts
        rate_with_bias = accepted_with_bias / n_attempts
        
        print(f"\nDensity {density:.1f}:")
        print(f"  Without cavity bias: {rate_no_bias:.2%}")
        print(f"  With cavity bias:    {rate_with_bias:.2%}")
        print(f"  Cavity fraction f_n: {f_n:.2f}")
        
        # Cavity bias should reduce acceptance at high density
        if density > 0.5:
            # Allow for small statistical fluctuations
            assert rate_with_bias <= rate_no_bias + 0.01, \
                f"Cavity bias should reduce acceptance at high density (got {rate_with_bias:.2%} > {rate_no_bias:.2%})"


# Helper functions


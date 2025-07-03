# tests/simulation/energyOpenmm/nbfix_parameters.py
"""
NBFIX parameter tests

This module tests:
1. NBFIX parameters correctly override Lorentz-Berthelot mixing rules
2. Multiple NBFIX pairs are handled correctly
"""

import pytest
import math
import pygcmc
import os
import warnings

# Filter SWIG-related warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)

try:
    from openmm import *
    from openmm.app import *
    from openmm.unit import *
    HAS_OPENMM = True
except ImportError:
    HAS_OPENMM = False

# Constants
ANGSTROM_TO_NM = 0.1
KCAL_TO_KJ = 4.184
kC = 138.935456  # Coulomb constant in kJ·nm/mol/e^2

# Get test data directory
TEST_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(TEST_DIR)), "data")

from .nbfix_setup import *

@pytest.mark.skipif(not HAS_OPENMM, reason="OpenMM not available")
def test_nbfix_overrides_lj_combination():
    """Test that NBFIX parameters override default Lorentz-Berthelot mixing rules"""
    
    # Create system and force field
    mc_state = create_simple_ion_water_mcstate()
    forcefield = create_nbfix_forcefield()
    
    # Apply force field with NBFIX
    mc_state = setup_mcstate_with_forcefield(mc_state, forcefield)
    
    # Find SOD and CLA type indices
    sod_idx = -1
    cla_idx = -1
    for i, atom_type in enumerate(mc_state.atomTypes.atomTypes):
        if atom_type == "SOD":
            sod_idx = i
        elif atom_type == "CLA":
            cla_idx = i
    
    
    assert sod_idx >= 0 and cla_idx >= 0, "SOD and CLA types not found"
    
    # Get the interaction parameters
    n_types = mc_state.forcefield.numTotalTypes
    nbfix_idx = sod_idx * n_types + cla_idx
    
    # Expected NBFIX values from toppar_water_ions.str
    # SOD    CLA      -0.083875   3.731
    expected_eps = 0.083875 * KCAL_TO_KJ  # Convert to kJ/mol
    expected_sigma = 3.731 / math.pow(2.0, 1.0/6.0) * ANGSTROM_TO_NM  # Convert Rmin to sigma
    
    
    actual_eps = mc_state.forcefield.ljEps[nbfix_idx]
    actual_sigma = mc_state.forcefield.ljSigma[nbfix_idx]
    
    print(f"\nNBFIX test for SOD-CLA:")
    print(f"Expected epsilon: {expected_eps:.6f} kJ/mol")
    print(f"Actual epsilon: {actual_eps:.6f} kJ/mol")
    print(f"Expected sigma: {expected_sigma:.6f} nm")
    print(f"Actual sigma: {actual_sigma:.6f} nm")
    
    # Check if NBFIX values are used (within numerical tolerance)
    assert abs(actual_eps - expected_eps) < 1e-3, f"NBFIX epsilon mismatch"
    assert abs(actual_sigma - expected_sigma) < 1e-4, f"NBFIX sigma mismatch"
    
    # Now check that standard mixing would give different values
    sod_lj = forcefield.get_lj_params("SOD")
    cla_lj = forcefield.get_lj_params("CLA")
    
    # Calculate what standard mixing would give
    sod_sigma = 2.0 * sod_lj.rmin_half / math.pow(2.0, 1.0/6.0) * ANGSTROM_TO_NM
    cla_sigma = 2.0 * cla_lj.rmin_half / math.pow(2.0, 1.0/6.0) * ANGSTROM_TO_NM
    sod_eps = abs(sod_lj.epsilon) * KCAL_TO_KJ
    cla_eps = abs(cla_lj.epsilon) * KCAL_TO_KJ
    
    standard_sigma = (sod_sigma + cla_sigma) / 2.0
    standard_eps = math.sqrt(sod_eps * cla_eps)
    
    print(f"\nStandard mixing would give:")
    print(f"Standard epsilon: {standard_eps:.6f} kJ/mol")
    print(f"Standard sigma: {standard_sigma:.6f} nm")
    print(f"Difference in epsilon: {abs(actual_eps - standard_eps):.6f} kJ/mol")
    print(f"Difference in sigma: {abs(actual_sigma - standard_sigma):.6f} nm")
    
    # NBFIX should be different from standard mixing  
    # For SOD-CLA, the epsilon is very close to standard mixing, but sigma is different
    assert abs(actual_sigma - standard_sigma) > 0.001, "NBFIX sigma should differ from standard mixing"

def test_multiple_nbfix_pairs():
    """Test system with multiple NBFIX corrections"""
    
    # Create force field
    forcefield = create_nbfix_forcefield()
    
    # Check that we have multiple NBFIX pairs
    nbfix_pairs = []
    test_types = ["SOD", "CLA", "OC", "OS", "ON3"]
    
    for type1 in test_types:
        for type2 in test_types:
            if type1 <= type2:  # Avoid duplicates
                result = forcefield.get_nbfix(type1, type2)
                # get_nbfix returns a tuple (epsilon, rmin, has_nbfix)
                if len(result) >= 3 and result[2]:  # has_nbfix is True
                    epsilon = result[0]
                    rmin = result[1]
                    nbfix_pairs.append((type1, type2, epsilon, rmin))
    
    print(f"\nFound {len(nbfix_pairs)} NBFIX pairs:")
    for type1, type2, epsilon, rmin in nbfix_pairs:
        print(f"  {type1}-{type2}: epsilon={epsilon:.5f} kcal/mol, Rmin={rmin:.3f} Å")
    
    # We should have at least SOD-CLA from the water_ions file
    assert len(nbfix_pairs) > 0, "Should have at least one NBFIX pair"
    
    # Check SOD-CLA specifically (manually added)
    sod_cla_found = False
    for type1, type2, epsilon, rmin in nbfix_pairs:
        if (type1 == "SOD" and type2 == "CLA") or (type1 == "CLA" and type2 == "SOD"):
            sod_cla_found = True
            assert abs(epsilon - (-0.083875)) < 1e-4
            assert abs(rmin - 3.731) < 1e-2
            print(f"  Verified {type1}-{type2} NBFIX matches expected values")
    
    assert sod_cla_found, "SOD-CLA NBFIX pair not found"


if __name__ == "__main__":
    test_nbfix_overrides_lj_combination()
    test_multiple_nbfix_pairs()
    print("\nAll NBFIX parameter tests passed!")

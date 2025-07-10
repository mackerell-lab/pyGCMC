# tests/simulation/energyOpenmm/nbfix_pbc.py
"""
NBFIX periodic boundary condition tests

This module tests:
1. NBFIX energy calculation with atoms separated by periodic boundary
2. Correct application of minimum image convention with NBFIX parameters
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
def test_nbfix_across_pbc():
    """Test NBFIX energy calculation with atoms separated by periodic boundary
    
    This tests that PyGCMC correctly applies NBFIX parameters when atoms interact
    across periodic boundaries, which is a common source of errors.
    """
    
    # Create MCState with atoms near box boundaries
    state = pygcmc.MCState()
    
    # Small box to force PBC interactions
    state.info.box = [2.0, 2.0, 2.0]  # 2 nm box
    state.info.cutoff = 1.2  # 1.2 nm cutoff
    
    # Define atom types
    sod_idx = state.atomTypes.get_or_add_type("SOD")
    cla_idx = state.atomTypes.get_or_add_type("CLA")
    
    # Place ions near opposite box edges
    # They are 1.5 nm apart in real space, but only 0.5 nm via PBC
    atom_na = pygcmc.MCAtom()
    atom_na.x = 0.25  # Near left edge
    atom_na.y = 1.0
    atom_na.z = 1.0
    atom_na.charge = 1.0
    atom_na.type = sod_idx
    
    atom_cl = pygcmc.MCAtom()
    atom_cl.x = 1.75  # Near right edge
    atom_cl.y = 1.0
    atom_cl.z = 1.0
    atom_cl.charge = -1.0
    atom_cl.type = cla_idx
    
    state.atoms = [atom_na, atom_cl]
    state.activeAtomCount = 2
    
    # Add residues
    res_na = pygcmc.MCResidue()
    res_na.active = True
    res_na.atomStart = 0
    res_na.atomCount = 1
    res_na.type = 0
    
    res_cl = pygcmc.MCResidue()
    res_cl.active = True
    res_cl.atomStart = 1
    res_cl.atomCount = 1
    res_cl.type = 1
    
    state.residues = [res_na, res_cl]
    state.activeResidueCount = 2
    
    # Create and apply force field with NBFIX
    forcefield = create_nbfix_forcefield()
    state = setup_mcstate_with_forcefield(state, forcefield)
    
    # Calculate PyGCMC energy
    pygcmc.computeSystemEnergyPBCCutoff(state)
    pygcmc_total = 0.0
    for res in state.residues[:state.activeResidueCount]:
        pygcmc_total += res.energy_vdw + res.energy_elec
    pygcmc_corrected = pygcmc_total / 2.0
    
    # Calculate expected energy manually
    # Distance via PBC: 2.0 - 1.5 = 0.5 nm
    pbc_distance = 0.5  # nm
    
    # Get NBFIX parameters for SOD-CLA
    nbfix_result = forcefield.get_nbfix("SOD", "CLA")
    nbfix_eps = abs(nbfix_result[0]) * KCAL_TO_KJ
    nbfix_rmin = nbfix_result[1]
    nbfix_sigma = nbfix_rmin / math.pow(2.0, 1.0/6.0) * ANGSTROM_TO_NM
    
    # Calculate expected energies
    expected_vdw = 4 * nbfix_eps * (math.pow(nbfix_sigma/pbc_distance, 12) - math.pow(nbfix_sigma/pbc_distance, 6))
    expected_elec = kC * 1.0 * (-1.0) / pbc_distance
    expected_total = expected_vdw + expected_elec
    
    print(f"\n=== PBC NBFIX Test ===")
    print(f"Box size: {state.info.box}")
    print(f"Atom positions: Na+ at ({atom_na.x}, {atom_na.y}, {atom_na.z}), Cl- at ({atom_cl.x}, {atom_cl.y}, {atom_cl.z})")
    print(f"Direct distance: 1.5 nm")
    print(f"PBC distance: {pbc_distance} nm")
    print(f"\nExpected energies:")
    print(f"  VDW: {expected_vdw:.6f} kJ/mol")
    print(f"  Electrostatic: {expected_elec:.6f} kJ/mol")
    print(f"  Total: {expected_total:.6f} kJ/mol")
    print(f"\nPyGCMC energy (corrected): {pygcmc_corrected:.6f} kJ/mol")
    
    # At 0.5 nm distance, check the energy
    # The energy should be negative (attractive) due to opposite charges
    assert pygcmc_corrected < 0, "Energy should be attractive for opposite charges at 0.5 nm"
    
    # Check that the energy is in the right ballpark
    # Allow larger tolerance due to potential cutoff effects
    rel_error = abs(pygcmc_corrected - expected_total) / abs(expected_total) if abs(expected_total) > 1e-10 else 0
    assert rel_error < 0.05, f"PBC energy mismatch: PyGCMC={pygcmc_corrected}, Expected={expected_total}"



if __name__ == "__main__":
    test_nbfix_across_pbc()
    print("\nAll NBFIX PBC tests passed!")

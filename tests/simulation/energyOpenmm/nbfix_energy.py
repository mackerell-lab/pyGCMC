# tests/simulation/energyOpenmm/nbfix_energy.py
"""
NBFIX energy tests

This module tests:
1. Ion-water interaction energies with NBFIX
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
def test_ion_water_nbfix_energy():
    """Test ion-water interaction energies with NBFIX"""
    
    # Create system and force field
    mc_state = create_simple_ion_water_mcstate()
    forcefield = create_nbfix_forcefield()
    
    # Apply force field with NBFIX
    mc_state = setup_mcstate_with_forcefield(mc_state, forcefield)
    
    # Calculate energy with PyGCMC
    pygcmc.computeSystemEnergyPBCCutoff(mc_state)
    
    # Get residue energies
    total_vdw = 0.0
    total_elec = 0.0
    for i in range(mc_state.activeResidueCount):
        total_vdw += mc_state.residues[i].energy_vdw
        total_elec += mc_state.residues[i].energy_elec
    total_energy = total_vdw + total_elec
    
    print(f"\nPyGCMC Energy Components:")
    print(f"VDW energy: {total_vdw:.6f} kJ/mol")
    print(f"Electrostatic energy: {total_elec:.6f} kJ/mol")
    print(f"Total energy: {total_energy:.6f} kJ/mol")
    
    # Detailed residue breakdown
    residue_names = ["SOD", "TIP3", "CLA"]
    for i, res in enumerate(mc_state.residues[:mc_state.activeResidueCount]):
        print(f"\nResidue {i} ({residue_names[i]}):")
        print(f"  VDW: {res.energy_vdw:.6f} kJ/mol")
        print(f"  Elec: {res.energy_elec:.6f} kJ/mol")
        print(f"  Total: {res.energy_vdw + res.energy_elec:.6f} kJ/mol")
    
    # Check that energies are reasonable
    assert abs(total_elec) > 0.1, "Should have non-zero electrostatic energy"
    assert total_energy < 0, "Ion-water interaction should be favorable (negative)"



if __name__ == "__main__":
    test_ion_water_nbfix_energy()
    print("\nAll NBFIX energy tests passed!")

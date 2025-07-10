# tests/simulation/energyOpenmm/nbfix_openmm.py
"""
NBFIX OpenMM comparison tests

This module tests:
1. Energy consistency between PyGCMC and OpenMM for NBFIX systems
2. VDW and electrostatic component comparisons
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
from .nbfix_openmm_helpers import *

@pytest.mark.skipif(not HAS_OPENMM, reason="OpenMM not available")  
def test_nbfix_energy_vs_openmm():
    """Compare NBFIX energy calculations between PyGCMC and OpenMM"""
    
    # Create system and force field
    mc_state = create_simple_ion_water_mcstate()
    forcefield = create_nbfix_forcefield()
    
    # Apply force field with NBFIX
    mc_state = setup_mcstate_with_forcefield(mc_state, forcefield)
    pygcmc.computeSystemEnergyPBCCutoff(mc_state)
    
    pygcmc_vdw = 0.0
    pygcmc_elec = 0.0
    for i in range(mc_state.activeResidueCount):
        pygcmc_vdw += mc_state.residues[i].energy_vdw
        pygcmc_elec += mc_state.residues[i].energy_elec
    pygcmc_total = pygcmc_vdw + pygcmc_elec
    
    # Set up OpenMM system
    omm_system, positions = setup_openmm_system_with_nbfix(mc_state, forcefield)
    
    # Calculate OpenMM energy
    integrator = VerletIntegrator(0.001 * picoseconds)
    platform = Platform.getPlatformByName('Reference')
    context = Context(omm_system, integrator, platform)
    context.setPositions(positions)
    
    state = context.getState(getEnergy=True)
    omm_energy = state.getPotentialEnergy().value_in_unit(kilojoules_per_mole)
    
    # IMPORTANT: PyGCMC's computeSystemEnergy* functions count each interaction twice
    # (both i-j and j-i), while OpenMM counts each pair once. We divide by 2 to compensate.
    # This is a known behavior of PyGCMC's system energy calculation.
    # NOTE: computeMovementEnergy* functions do NOT double count.
    pygcmc_corrected = pygcmc_total / 2.0
    
    print(f"\nEnergy comparison:")
    print(f"PyGCMC total (before correction): {pygcmc_total:.6f} kJ/mol")
    print(f"PyGCMC total (after /2 correction): {pygcmc_corrected:.6f} kJ/mol") 
    print(f"OpenMM total: {omm_energy:.6f} kJ/mol")
    
    # Allow 1% relative error due to implementation differences
    rel_error = abs(pygcmc_corrected - omm_energy) / abs(omm_energy) if abs(omm_energy) > 1e-10 else 0
    assert rel_error < 0.01, f"Energy mismatch: PyGCMC={pygcmc_corrected}, OpenMM={omm_energy}"
    
    # Test moving an atom and recalculating
    print("\n--- Testing energy change after moving Cl- ---")
    
    # Move chloride ion farther away
    mc_state.atoms[4].z = 2.5  # Move from 2.0 to 2.5 nm
    
    # Recalculate PyGCMC energy
    pygcmc.computeSystemEnergyPBCCutoff(mc_state)
    pygcmc_total_new = 0.0
    for i in range(mc_state.activeResidueCount):
        pygcmc_total_new += mc_state.residues[i].energy_vdw + mc_state.residues[i].energy_elec
    
    # Update OpenMM positions
    new_positions = list(positions)
    new_positions[4] = Vec3(1.5, 1.5, 2.5) * nanometers
    context.setPositions(new_positions)
    
    state_new = context.getState(getEnergy=True)
    omm_energy_new = state_new.getPotentialEnergy().value_in_unit(kilojoules_per_mole)
    
    # Compare energy changes - need to correct for factor of 2
    pygcmc_delta = (pygcmc_total_new - pygcmc_total) / 2.0
    omm_delta = omm_energy_new - omm_energy
    
    print(f"\nEnergy change after moving Cl-:")
    print(f"PyGCMC delta: {pygcmc_delta:.6f} kJ/mol")
    print(f"OpenMM delta: {omm_delta:.6f} kJ/mol")
    print(f"Difference in deltas: {abs(pygcmc_delta - omm_delta):.6f} kJ/mol")
    
    # Energy changes should be consistent
    assert abs(pygcmc_delta - omm_delta) < 0.1, "Energy changes should be consistent"
    
    # Moving opposite charges apart should increase energy
    assert pygcmc_delta > 0, "Moving Na+ and Cl- apart should increase energy"
    
    # Clean up OpenMM context
    del context

@pytest.mark.skipif(not HAS_OPENMM, reason="OpenMM not available")
def test_nbfix_energy_components_separately():
    """Compare VDW and electrostatic energy components separately between PyGCMC and OpenMM
    
    This test provides more detailed validation by comparing energy components individually,
    which helps identify whether discrepancies come from VDW or electrostatic calculations.
    """
    
    # Create system and force field
    mc_state = create_simple_ion_water_mcstate()
    forcefield = create_nbfix_forcefield()
    
    # Apply force field with NBFIX
    mc_state = setup_mcstate_with_forcefield(mc_state, forcefield)
    pygcmc.computeSystemEnergyPBCCutoff(mc_state)
    
    # Get PyGCMC energy components
    pygcmc_vdw = 0.0
    pygcmc_elec = 0.0
    for i in range(mc_state.activeResidueCount):
        pygcmc_vdw += mc_state.residues[i].energy_vdw
        pygcmc_elec += mc_state.residues[i].energy_elec
    
    # IMPORTANT: PyGCMC's computeSystemEnergy* functions count each interaction twice
    # This is a known behavior that requires division by 2 for comparison with OpenMM
    pygcmc_vdw_corrected = pygcmc_vdw / 2.0
    pygcmc_elec_corrected = pygcmc_elec / 2.0
    
    # Set up OpenMM system with separate forces
    omm_system, positions = setup_openmm_system_with_nbfix(mc_state, forcefield, separate_forces=True)
    
    # Calculate OpenMM energy components
    integrator = VerletIntegrator(0.001 * picoseconds)
    platform = Platform.getPlatformByName('Reference')
    context = Context(omm_system, integrator, platform)
    context.setPositions(positions)
    
    # Get VDW energy (force group 1)
    state_vdw = context.getState(getEnergy=True, groups={1})
    omm_vdw = state_vdw.getPotentialEnergy().value_in_unit(kilojoules_per_mole)
    
    # Get electrostatic energy (force group 2)
    state_elec = context.getState(getEnergy=True, groups={2})
    omm_elec = state_elec.getPotentialEnergy().value_in_unit(kilojoules_per_mole)
    
    # Get total energy
    state_total = context.getState(getEnergy=True)
    omm_total = state_total.getPotentialEnergy().value_in_unit(kilojoules_per_mole)
    
    print(f"\n=== Detailed Energy Component Comparison ===")
    print(f"VDW Energy:")
    print(f"  PyGCMC (raw): {pygcmc_vdw:.6f} kJ/mol")
    print(f"  PyGCMC (corrected /2): {pygcmc_vdw_corrected:.6f} kJ/mol")
    print(f"  OpenMM: {omm_vdw:.6f} kJ/mol")
    print(f"  Relative error: {abs(pygcmc_vdw_corrected - omm_vdw) / abs(omm_vdw) * 100:.2f}%")
    
    print(f"\nElectrostatic Energy:")
    print(f"  PyGCMC (raw): {pygcmc_elec:.6f} kJ/mol")
    print(f"  PyGCMC (corrected /2): {pygcmc_elec_corrected:.6f} kJ/mol")
    print(f"  OpenMM: {omm_elec:.6f} kJ/mol")
    print(f"  Relative error: {abs(pygcmc_elec_corrected - omm_elec) / abs(omm_elec) * 100:.2f}%")
    
    print(f"\nTotal Energy:")
    print(f"  PyGCMC (corrected): {pygcmc_vdw_corrected + pygcmc_elec_corrected:.6f} kJ/mol")
    print(f"  OpenMM: {omm_total:.6f} kJ/mol")
    
    # Compare components separately with 1% tolerance
    rel_tol = 0.01
    
    # Check VDW energy
    if abs(omm_vdw) > 1e-6:  # Only check relative error if VDW is non-negligible
        vdw_rel_error = abs(pygcmc_vdw_corrected - omm_vdw) / abs(omm_vdw) if abs(omm_vdw) > 1e-10 else 0
        assert vdw_rel_error < rel_tol, f"VDW energy mismatch: PyGCMC={pygcmc_vdw_corrected}, OpenMM={omm_vdw}"
    else:
        assert abs(pygcmc_vdw_corrected - omm_vdw) < 1e-3, f"VDW energy should be near zero"
    
    # Check electrostatic energy
    elec_rel_error = abs(pygcmc_elec_corrected - omm_elec) / abs(omm_elec) if abs(omm_elec) > 1e-10 else 0
    assert elec_rel_error < rel_tol, f"Electrostatic energy mismatch: PyGCMC={pygcmc_elec_corrected}, OpenMM={omm_elec}"
    
    # Check total energy consistency
    omm_sum = omm_vdw + omm_elec
    assert abs(omm_total - omm_sum) < 1e-6, "OpenMM total should equal sum of components"
    
    # Clean up
    del context


if __name__ == "__main__":
    if HAS_OPENMM:
        test_nbfix_energy_vs_openmm()
        test_nbfix_energy_components_separately()
        print("\nAll NBFIX OpenMM comparison tests passed!")
    else:
        print("OpenMM not available - skipping tests")

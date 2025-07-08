"""
Test full nonbonded PME interactions (electrostatic + LJ)
"""

import pytest
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pygcmc
from pygcmc import initializePMEParameters, computeSystemEnergyPME, computeSystemEnergyPMEFixed

from .full_nonbonded_helpers import create_test_system_with_lj, calculate_openmm_full_energy, OPENMM_AVAILABLE

@pytest.mark.skipif(not OPENMM_AVAILABLE, reason="OpenMM not available")
def test_full_nonbonded_pme_lj():
    """Test full nonbonded energy (PME + LJ) comparison with OpenMM"""
    
    # Create test system
    state = create_test_system_with_lj()
    
    print(f"\nTest system with LJ and charges:")
    print(f"  Atoms: {state.activeAtomCount}")
    print(f"  Box: {state.info.box[0]} nm")
    print(f"  Cutoff: {state.info.cutoff} nm")
    print(f"  LJ parameters: eps={state.forcefield.ljEps[0]} kJ/mol, sigma={state.forcefield.ljSigma[0]} nm")
    
    # PME parameters
    alpha = 3.2
    mesh_size = [32, 32, 32]
    spline_order = 4
    
    # Initialize PyGCMC PME
    pygcmc.setPMEParameters(alpha, mesh_size, spline_order)
    initializePMEParameters(state.info.cutoff, state.info.box, alpha, mesh_size, spline_order)
    
    # Calculate PyGCMC energy using fixed function
    elec_total, vdw_total, pme_dict = computeSystemEnergyPMEFixed(state)
    pygcmc_elec = pme_dict['total']
    pygcmc_lj = vdw_total
    
    pygcmc_total = pygcmc_elec + pygcmc_lj
    
    # Calculate OpenMM energy
    openmm_total, openmm_elec, openmm_lj = calculate_openmm_full_energy(state, alpha)
    
    if openmm_total is not None:
        print(f"\nEnergy comparison:")
        print(f"  PyGCMC:")
        print(f"    Electrostatic: {pygcmc_elec:.6f} kJ/mol")
        print(f"    LJ: {pygcmc_lj:.6f} kJ/mol")
        print(f"    Total: {pygcmc_total:.6f} kJ/mol")
        print(f"  OpenMM:")
        print(f"    Electrostatic: {openmm_elec:.6f} kJ/mol")
        print(f"    LJ: {openmm_lj:.6f} kJ/mol")
        print(f"    Total: {openmm_total:.6f} kJ/mol")
        
        # Calculate differences
        elec_diff = abs(pygcmc_elec - openmm_elec)
        elec_rel = elec_diff / abs(openmm_elec) if openmm_elec != 0 else 0
        
        lj_diff = abs(pygcmc_lj - openmm_lj)
        lj_rel = lj_diff / abs(openmm_lj) if openmm_lj != 0 else 0
        
        total_diff = abs(pygcmc_total - openmm_total)
        total_rel = total_diff / abs(openmm_total) if openmm_total != 0 else 0
        
        print(f"\nRelative differences:")
        print(f"  Electrostatic: {elec_rel*100:.3f}%")
        print(f"  LJ: {lj_rel*100:.3f}%")
        print(f"  Total: {total_rel*100:.3f}%")
        
        # Assertions
        # PME implementations can differ slightly between PyGCMC and OpenMM
        assert elec_rel < 0.05, f"Electrostatic energies differ by {elec_rel*100:.3f}% (> 5%)"
        assert lj_rel < 0.01, f"LJ energies differ by {lj_rel*100:.3f}% (> 1%)"
        assert total_rel < 0.05, f"Total energies differ by {total_rel*100:.3f}% (> 5%)"
        
        print("\n✓ Full nonbonded energy calculation matches OpenMM!")
    else:
        print("\nOpenMM not available for comparison")



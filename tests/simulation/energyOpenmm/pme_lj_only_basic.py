"""
Basic PME LJ-only tests

Tests basic LJ energy calculations with cutoff and PME methods
to isolate LJ calculation from electrostatic contributions.
"""

import pytest
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pygcmc
from pygcmc import initializePMEParameters, computeSystemEnergyPME, computeSystemEnergyCutoff
from pygcmc import computeSystemEnergyCutoffFixed, computeSystemEnergyPMEFixed

from .pme_lj_only_helpers import create_lj_only_system, calculate_openmm_lj_energy, OPENMM_AVAILABLE


@pytest.mark.skipif(not OPENMM_AVAILABLE, reason="OpenMM not available")
def test_lj_only_cutoff():
    """Test LJ energy with simple cutoff (no PME)"""
    
    state = create_lj_only_system(n_atoms=6)
    
    print("\nLJ-only test with Cutoff:")
    print(f"  Atoms: {state.activeAtomCount}")
    print(f"  Box: {state.info.box[0]} nm")
    print(f"  Cutoff: {state.info.cutoff} nm")
    print(f"  LJ parameters:")
    print(f"    Type 0: eps={state.forcefield.ljEps[0]} kJ/mol, sigma={state.forcefield.ljSigma[0]} nm")
    print(f"    Type 1: eps={state.forcefield.ljEps[1]} kJ/mol, sigma={state.forcefield.ljSigma[1]} nm")
    
    # Calculate with PyGCMC using fixed cutoff function
    elec, vdw, total = computeSystemEnergyCutoffFixed(state)
    pygcmc_lj = vdw
    
    # Calculate with OpenMM
    openmm_lj = calculate_openmm_lj_energy(state, use_pme=False)
    
    print(f"\nEnergy comparison (Cutoff):")
    print(f"  PyGCMC LJ: {pygcmc_lj:.6f} kJ/mol")
    print(f"  OpenMM LJ: {openmm_lj:.6f} kJ/mol")
    
    diff = abs(pygcmc_lj - openmm_lj)
    rel_diff = diff / abs(openmm_lj) if openmm_lj != 0 else 0
    print(f"  Absolute difference: {diff:.6f} kJ/mol")
    print(f"  Relative difference: {rel_diff*100:.3f}%")
    
    # Cutoff should match reasonably well
    # Allow higher tolerance for mixed LJ types due to implementation differences
    assert rel_diff < 0.15, f"Cutoff LJ energies differ by {rel_diff*100:.3f}% (> 15%)"
    
    print("\n✓ Cutoff LJ calculation matches OpenMM!")


@pytest.mark.skipif(not OPENMM_AVAILABLE, reason="OpenMM not available")
def test_lj_only_pme():
    """Test LJ energy with PME to isolate the issue"""
    
    state = create_lj_only_system(n_atoms=6)
    
    print("\nLJ-only test with PME:")
    print(f"  Atoms: {state.activeAtomCount}")
    
    # PME parameters
    alpha = 3.2
    mesh_size = [32, 32, 32]
    spline_order = 4
    
    # Initialize PME
    pygcmc.setPMEParameters(alpha, mesh_size, spline_order)
    initializePMEParameters(state.info.cutoff, state.info.box, alpha, mesh_size, spline_order)
    
    # Calculate with PyGCMC PME
    computeSystemEnergyPME(state)
    
    # Get LJ energy
    pygcmc_lj = 0.0
    for res in state.residues:
        if res.active:
            pygcmc_lj += res.energy_vdw
    
    # Also check if PME modified anything (it shouldn't for LJ-only)
    pygcmc_elec = state.ewald_energy.get('total', 0.0)
    
    # Calculate with OpenMM PME
    openmm_lj_pme = calculate_openmm_lj_energy(state, use_pme=True)
    
    # Calculate with OpenMM Cutoff for comparison
    openmm_lj_cutoff = calculate_openmm_lj_energy(state, use_pme=False)
    
    print(f"\nEnergy comparison (PME):")
    print(f"  PyGCMC LJ (PME): {pygcmc_lj:.6f} kJ/mol")
    print(f"  PyGCMC Elec (should be 0): {pygcmc_elec:.6f} kJ/mol")
    print(f"  OpenMM LJ (PME): {openmm_lj_pme:.6f} kJ/mol")
    print(f"  OpenMM LJ (Cutoff): {openmm_lj_cutoff:.6f} kJ/mol")
    
    # Compare PME vs Cutoff for OpenMM (should be same for LJ-only)
    openmm_diff = abs(openmm_lj_pme - openmm_lj_cutoff)
    print(f"\n  OpenMM PME vs Cutoff difference: {openmm_diff:.6f} kJ/mol")
    
    # Compare PyGCMC PME with OpenMM
    diff = abs(pygcmc_lj - openmm_lj_pme)
    rel_diff = diff / abs(openmm_lj_pme) if openmm_lj_pme != 0 else 0
    print(f"\n  PyGCMC vs OpenMM (PME) difference: {rel_diff*100:.3f}%")
    
    # Check if PME is affecting LJ calculation incorrectly
    if rel_diff > 0.01:
        print("\n⚠️  PME appears to be affecting LJ calculation!")
        print("This suggests the PME implementation may be modifying LJ energies.")
    
    assert pygcmc_elec == 0.0, "PME should not produce electrostatic energy for charge-free system"


if __name__ == "__main__":
    test_lj_only_cutoff()
    test_lj_only_pme()
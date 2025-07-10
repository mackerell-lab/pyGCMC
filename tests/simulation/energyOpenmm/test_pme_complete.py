"""
Test PME Complete functions that include intramolecular interactions

These tests verify that computeSystemEnergyPMEComplete and computeSystemEnergyCutoffComplete
properly match OpenMM by including intramolecular LJ interactions.
"""

import pytest
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pygcmc
from pygcmc import computeSystemEnergyPMEComplete, computeSystemEnergyCutoffComplete
from pygcmc import initializePMEParameters

from .pme_lj_only_helpers import create_lj_only_system_with_molecules, calculate_openmm_lj_energy, OPENMM_AVAILABLE
from .pme_medium_complexity import create_medium_complexity_system, calculate_openmm_energy_medium


@pytest.mark.skipif(not OPENMM_AVAILABLE, reason="OpenMM not available")
def test_complete_cutoff_lj_only():
    """Test that Complete cutoff matches OpenMM exactly for LJ-only systems"""
    
    # Create system with multi-atom molecules
    state = create_lj_only_system_with_molecules(n_molecules=3, atoms_per_mol=2)
    
    print("\nComplete Cutoff LJ-only test:")
    print(f"  System: {state.activeAtomCount} atoms in {len([r for r in state.residues if r.active])} residues")
    print(f"  Box: {state.info.box[0]} nm")
    print(f"  Cutoff: {state.info.cutoff} nm")
    
    # Calculate with PyGCMC Complete cutoff
    elec, vdw, total = computeSystemEnergyCutoffComplete(state)
    
    # Calculate with OpenMM
    openmm_lj = calculate_openmm_lj_energy(state, use_pme=False)
    
    print(f"\nEnergy comparison:")
    print(f"  PyGCMC Complete: {vdw:.6f} kJ/mol")
    print(f"  OpenMM:          {openmm_lj:.6f} kJ/mol")
    
    rel_diff = abs(vdw - openmm_lj) / abs(openmm_lj) if openmm_lj != 0 else 0
    print(f"  Relative difference: {rel_diff*100:.3f}%")
    
    # Should match very closely now
    assert rel_diff < 0.01, f"Complete cutoff LJ differs by {rel_diff*100:.3f}% (> 1%)"
    assert elec == 0.0, "Should have no electrostatic energy for charge-free system"


@pytest.mark.skipif(not OPENMM_AVAILABLE, reason="OpenMM not available")
def test_complete_pme_lj_only():
    """Test that Complete PME matches OpenMM for LJ-only systems"""
    
    # Create system with multi-atom molecules
    state = create_lj_only_system_with_molecules(n_molecules=3, atoms_per_mol=2)
    
    print("\nComplete PME LJ-only test:")
    print(f"  System: {state.activeAtomCount} atoms in {len([r for r in state.residues if r.active])} residues")
    
    # Initialize PME
    alpha = 3.2
    mesh_size = [32, 32, 32]
    spline_order = 4
    
    pygcmc.setPMEParameters(alpha, mesh_size, spline_order)
    initializePMEParameters(state.info.cutoff, state.info.box, alpha, mesh_size, spline_order)
    
    # Calculate with PyGCMC Complete PME
    elec, vdw, total = computeSystemEnergyPMEComplete(state)
    
    # Calculate with OpenMM PME
    openmm_lj = calculate_openmm_lj_energy(state, use_pme=True)
    
    print(f"\nEnergy comparison:")
    print(f"  PyGCMC Complete PME: {vdw:.6f} kJ/mol")
    print(f"  OpenMM PME:          {openmm_lj:.6f} kJ/mol")
    
    rel_diff = abs(vdw - openmm_lj) / abs(openmm_lj) if openmm_lj != 0 else 0
    print(f"  Relative difference: {rel_diff*100:.3f}%")
    
    # Should match very closely
    assert rel_diff < 0.01, f"Complete PME LJ differs by {rel_diff*100:.3f}% (> 1%)"
    assert elec == 0.0, "Should have no electrostatic energy for charge-free system"


@pytest.mark.skipif(not OPENMM_AVAILABLE, reason="OpenMM not available")
def test_complete_pme_full_system():
    """Test Complete PME with full electrostatics and LJ"""
    
    # Create a medium complexity system
    state = create_medium_complexity_system()
    
    print("\nComplete PME full system test:")
    print(f"  System: {state.activeAtomCount} atoms")
    print(f"  Box: {state.info.box[0]} nm")
    
    # Initialize PME
    alpha = 5.6
    mesh_size = [32, 32, 32]
    spline_order = 4
    
    pygcmc.setPMEParameters(alpha, mesh_size, spline_order)
    initializePMEParameters(state.info.cutoff, state.info.box, alpha, mesh_size, spline_order)
    
    # Calculate with PyGCMC Complete PME
    elec, vdw, total = computeSystemEnergyPMEComplete(state)
    
    print(f"\nPyGCMC Complete PME energies:")
    print(f"  Electrostatic: {elec:.6f} kJ/mol")
    print(f"  VdW:           {vdw:.6f} kJ/mol")
    print(f"  Total:         {total:.6f} kJ/mol")
    
    # Calculate with OpenMM
    openmm_total = calculate_openmm_energy_medium(state, alpha)
    
    print(f"\nOpenMM PME total energy: {openmm_total:.6f} kJ/mol")
    
    # Check total energy
    total_diff = abs(total - openmm_total) / abs(openmm_total) if openmm_total != 0 else 0
    
    print(f"\nTotal energy relative difference: {total_diff*100:.3f}%")
    
    # Complete should match OpenMM total energy well
    assert total_diff < 0.05, f"Total energy differs by {total_diff*100:.3f}% (> 5%)"


@pytest.mark.skipif(not OPENMM_AVAILABLE, reason="OpenMM not available")
def test_complete_vs_fixed_difference():
    """Test that Complete includes intramolecular interactions that Fixed excludes"""
    
    # Create system with multi-atom molecules
    state = create_lj_only_system_with_molecules(n_molecules=2, atoms_per_mol=3)
    
    print("\nComplete vs Fixed comparison:")
    print(f"  System: {state.activeAtomCount} atoms in {len([r for r in state.residues if r.active])} residues")
    
    # Calculate with Fixed function (excludes intramolecular)
    from pygcmc import computeSystemEnergyCutoffFixed
    computeSystemEnergyCutoffFixed(state)
    vdw_fixed = sum(res.energy_vdw for res in state.residues if res.active) / 2.0
    
    # Calculate with Complete function (includes intramolecular)
    elec, vdw_complete, total = computeSystemEnergyCutoffComplete(state)
    
    print(f"\nEnergy comparison:")
    print(f"  Fixed (excludes intra):    {vdw_fixed:.6f} kJ/mol")
    print(f"  Complete (includes intra): {vdw_complete:.6f} kJ/mol")
    
    # The difference is the intramolecular contribution
    intra_contrib = vdw_complete - vdw_fixed
    print(f"  Intramolecular contribution: {intra_contrib:.6f} kJ/mol")
    
    # Complete should be different from Fixed (intramolecular can be repulsive or attractive)
    assert abs(intra_contrib) > 1e-6, "Complete should differ from Fixed by intramolecular contribution"
    
    # For close atoms, intramolecular can be positive (repulsive)
    print(f"  Note: Intramolecular is {'repulsive' if intra_contrib > 0 else 'attractive'}")


@pytest.mark.skipif(not OPENMM_AVAILABLE, reason="OpenMM not available")
def test_complete_consistency():
    """Test that Complete PME and Complete Cutoff give similar VdW energies"""
    
    # Create a simple system
    state = create_lj_only_system_with_molecules(n_molecules=4, atoms_per_mol=1)
    
    print("\nComplete consistency test:")
    print(f"  System: {state.activeAtomCount} atoms")
    print(f"  Force field types: {state.forcefield.numTotalTypes}")
    print(f"  LJ params length: eps={len(state.forcefield.ljEps)}, sigma={len(state.forcefield.ljSigma)}")
    
    # Calculate with Complete Cutoff
    elec_cut, vdw_cut, total_cut = computeSystemEnergyCutoffComplete(state)
    print(f"  Cutoff result: elec={elec_cut}, vdw={vdw_cut}, total={total_cut}")
    
    # Initialize PME
    alpha = 3.2
    mesh_size = [32, 32, 32]
    spline_order = 4
    
    pygcmc.setPMEParameters(alpha, mesh_size, spline_order)
    initializePMEParameters(state.info.cutoff, state.info.box, alpha, mesh_size, spline_order)
    
    # Calculate with Complete PME
    elec_pme, vdw_pme, total_pme = computeSystemEnergyPMEComplete(state)
    print(f"  PME result: elec={elec_pme}, vdw={vdw_pme}, total={total_pme}")
    
    print(f"\nVdW energy comparison:")
    print(f"  Complete Cutoff: {vdw_cut:.6f} kJ/mol")
    print(f"  Complete PME:    {vdw_pme:.6f} kJ/mol")
    
    # VdW should be identical (PME doesn't affect LJ)
    vdw_diff = abs(vdw_pme - vdw_cut)
    print(f"  Absolute difference: {vdw_diff:.6f} kJ/mol")
    
    assert vdw_diff < 1e-6, f"VdW energies should be identical, differ by {vdw_diff}"


if __name__ == "__main__":
    test_complete_cutoff_lj_only()
    test_complete_pme_lj_only()
    test_complete_pme_full_system()
    test_complete_vs_fixed_difference()
    test_complete_consistency()
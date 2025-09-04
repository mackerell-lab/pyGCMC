#!/usr/bin/env python3
"""Test the getTotalEnergyComponents function"""

import pytest
import pygcmc
import numpy as np

def test_get_total_energy_components_simple():
    """Test getTotalEnergyComponents with a simple two-atom system"""
    # Create a simple test system with two atoms
    state = pygcmc.MCState()
    state.info.cutoff = 12.0
    state.info.box = [30.0, 30.0, 30.0]
    
    # Add two atoms with charges and LJ parameters
    atom1 = pygcmc.MCAtom()
    atom1.x, atom1.y, atom1.z = 0.0, 0.0, 0.0
    atom1.charge = 1.0
    atom1.type = 0
    
    atom2 = pygcmc.MCAtom()
    atom2.x, atom2.y, atom2.z = 3.0, 0.0, 0.0
    atom2.charge = -1.0
    atom2.type = 0
    
    # Set the atoms and count
    state.atoms = [atom1, atom2]
    state.activeAtomCount = 2
    
    # Add residues
    res1 = pygcmc.MCResidue()
    res1.atomStart = 0
    res1.atomCount = 1
    res1.active = True
    
    res2 = pygcmc.MCResidue()
    res2.atomStart = 1
    res2.atomCount = 1
    res2.active = True
    
    state.residues = [res1, res2]
    state.activeResidueCount = 2
    
    # Set force field parameters
    state.forcefield.numTotalTypes = 1
    state.forcefield.ljSigma = [3.0]
    state.forcefield.ljEps = [0.1]
    
    # Compute energy
    pygcmc.computeSystemEnergy(state)
    
    # Get energy components
    elec, vdw = pygcmc.getTotalEnergyComponents(state)
    
    # Check that we get sensible values
    assert elec < 0, "Electrostatic energy should be negative for opposite charges"
    assert vdw == 0.0, "VdW energy should be 0 at this distance with sigma=3.0"
    
    # Check that the total matches the sum of residue energies
    total_elec_from_residues = sum(res.energy_elec for res in state.residues if res.active) / 2.0
    total_vdw_from_residues = sum(res.energy_vdw for res in state.residues if res.active) / 2.0
    
    assert abs(elec - total_elec_from_residues) < 1e-6, "Electrostatic energy mismatch"
    assert abs(vdw - total_vdw_from_residues) < 1e-6, "VdW energy mismatch"

def test_get_total_energy_components_with_vdw():
    """Test getTotalEnergyComponents with significant VdW interactions"""
    # Create a system where atoms are close enough for VdW interactions
    state = pygcmc.MCState()
    state.info.cutoff = 12.0
    state.info.box = [30.0, 30.0, 30.0]
    
    # Add two atoms closer together
    atom1 = pygcmc.MCAtom()
    atom1.x, atom1.y, atom1.z = 0.0, 0.0, 0.0
    atom1.charge = 0.5
    atom1.type = 0
    
    atom2 = pygcmc.MCAtom()
    atom2.x, atom2.y, atom2.z = 2.5, 0.0, 0.0  # Closer for VdW interaction
    atom2.charge = -0.5
    atom2.type = 0
    
    state.atoms = [atom1, atom2]
    state.activeAtomCount = 2
    
    # Add residues
    res1 = pygcmc.MCResidue()
    res1.atomStart = 0
    res1.atomCount = 1
    res1.active = True
    
    res2 = pygcmc.MCResidue()
    res2.atomStart = 1
    res2.atomCount = 1
    res2.active = True
    
    state.residues = [res1, res2]
    state.activeResidueCount = 2
    
    # Set force field parameters with significant epsilon
    state.forcefield.numTotalTypes = 1
    state.forcefield.ljSigma = [2.0]  # Smaller sigma
    state.forcefield.ljEps = [1.0]    # Larger epsilon
    
    # Compute energy
    pygcmc.computeSystemEnergy(state)
    
    # Get energy components
    elec, vdw = pygcmc.getTotalEnergyComponents(state)
    
    # Check that we get sensible values
    assert elec < 0, "Electrostatic energy should be negative for opposite charges"
    assert vdw != 0.0, "VdW energy should be non-zero at this distance"
    
    # At r=2.5 with sigma=2.0, we're at r/sigma = 1.25
    # This is past the minimum (r_min = 2^(1/6) * sigma ≈ 2.24)
    # So VdW should be slightly negative (attractive)
    assert vdw < 0, "VdW energy should be negative (attractive) at this distance"

def test_get_total_energy_components_multiple_residues():
    """Test getTotalEnergyComponents with multiple residues"""
    # Create a more complex system
    state = pygcmc.MCState()
    state.info.cutoff = 12.0
    state.info.box = [30.0, 30.0, 30.0]
    
    # Add 4 atoms in a square arrangement
    positions = [(0, 0, 0), (3, 0, 0), (3, 3, 0), (0, 3, 0)]
    charges = [1.0, -1.0, 1.0, -1.0]
    
    atoms = []
    for i, (pos, charge) in enumerate(zip(positions, charges)):
        atom = pygcmc.MCAtom()
        atom.x, atom.y, atom.z = pos
        atom.charge = charge
        atom.type = 0
        atoms.append(atom)
    
    state.atoms = atoms
    state.activeAtomCount = 4
    
    # Add residues (2 atoms per residue)
    residues = []
    for i in range(2):
        res = pygcmc.MCResidue()
        res.atomStart = i * 2
        res.atomCount = 2
        res.active = True
        residues.append(res)
    
    state.residues = residues
    state.activeResidueCount = 2
    
    # Set force field parameters
    state.forcefield.numTotalTypes = 1
    state.forcefield.ljSigma = [2.5]
    state.forcefield.ljEps = [0.5]
    
    # Compute energy
    pygcmc.computeSystemEnergy(state)
    
    # Get energy components
    elec, vdw = pygcmc.getTotalEnergyComponents(state)
    
    # Check that energies are non-zero
    assert elec != 0.0, "Electrostatic energy should be non-zero"
    assert vdw != 0.0, "VdW energy should be non-zero"
    
    # Verify consistency with residue energies
    total_elec_from_residues = sum(res.energy_elec for res in state.residues if res.active) / 2.0
    total_vdw_from_residues = sum(res.energy_vdw for res in state.residues if res.active) / 2.0
    
    assert abs(elec - total_elec_from_residues) < 1e-6, "Electrostatic energy mismatch"
    assert abs(vdw - total_vdw_from_residues) < 1e-6, "VdW energy mismatch"

if __name__ == "__main__":
    test_get_total_energy_components_simple()
    test_get_total_energy_components_with_vdw()
    test_get_total_energy_components_multiple_residues()
    print("All tests passed!")
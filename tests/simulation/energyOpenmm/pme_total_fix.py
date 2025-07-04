"""
Test PME Total field fix

Verify that PME total energy calculation is correct after bug fix
"""

import pytest
import numpy as np
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pygcmc
from pygcmc import MCState, MCAtom, MCResidue, MCForceField
from pygcmc import initializePMEParameters, computeSystemEnergyPME


def test_pme_total_single_residue():
    """Test PME total for single residue system (should be correct)"""
    state = MCState()
    state.info.box = [5.0, 5.0, 5.0]
    state.info.cutoff = 2.0
    
    ff = MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.ljEps = [0.0]  # No LJ
    ff.ljSigma = [0.3]
    state.forcefield = ff
    
    # 2 atoms in 1 residue
    positions = [[2.0, 2.5, 2.5], [3.0, 2.5, 2.5]]
    charges = [1.0, -1.0]
    
    atoms = []
    for i in range(2):
        atom = MCAtom()
        atom.x, atom.y, atom.z = positions[i]
        atom.charge = charges[i]
        atom.type = 0
        atoms.append(atom)
    
    state.atoms = atoms
    state.activeAtomCount = 2
    
    # Single residue
    res = MCResidue()
    res.active = True
    res.fixed = False
    res.atomStart = 0
    res.atomCount = 2
    res.type = 0
    
    state.residues = [res]
    state.activeResidueCount = 1
    
    # Initialize and compute
    initializePMEParameters(state.info.cutoff, state.info.box, 2.5)
    computeSystemEnergyPME(state)
    
    # Check that total equals sum of components
    total = state.ewald_energy.get('total', 0.0)
    expected = (state.ewald_energy.get('real_space', 0.0) + 
                state.ewald_energy.get('reciprocal', 0.0) + 
                state.ewald_energy.get('self', 0.0))
    
    assert abs(total - expected) < 1e-6, f"PME total {total} != sum {expected}"


def test_pme_total_multiple_residues():
    """Test PME total for multiple residue system (bug was here)"""
    state = MCState()
    state.info.box = [5.0, 5.0, 5.0]
    state.info.cutoff = 2.0
    
    ff = MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.ljEps = [0.0]  # No LJ
    ff.ljSigma = [0.3]
    state.forcefield = ff
    
    # 2 atoms in 2 residues
    positions = [[2.0, 2.5, 2.5], [3.0, 2.5, 2.5]]
    charges = [1.0, -1.0]
    
    atoms = []
    for i in range(2):
        atom = MCAtom()
        atom.x, atom.y, atom.z = positions[i]
        atom.charge = charges[i]
        atom.type = 0
        atoms.append(atom)
    
    state.atoms = atoms
    state.activeAtomCount = 2
    
    # Two residues
    residues = []
    for i in range(2):
        res = MCResidue()
        res.active = True
        res.fixed = False
        res.atomStart = i
        res.atomCount = 1
        res.type = 0
        residues.append(res)
    
    state.residues = residues
    state.activeResidueCount = 2
    
    # Initialize and compute
    initializePMEParameters(state.info.cutoff, state.info.box, 2.5)
    computeSystemEnergyPME(state)
    
    # Check that total equals sum of components
    total = state.ewald_energy.get('total', 0.0)
    expected = (state.ewald_energy.get('real_space', 0.0) + 
                state.ewald_energy.get('reciprocal', 0.0) + 
                state.ewald_energy.get('self', 0.0))
    
    assert abs(total - expected) < 1e-6, f"PME total {total} != sum {expected}"
    
    # Verify no huge offset (bug would give ~38000 kJ/mol offset)
    assert abs(total) < 1000, f"PME total {total} has unreasonable magnitude"


def test_pme_total_vs_ewald():
    """Test that PME and Ewald give consistent total energies"""
    from pygcmc import initializeEwaldParameters, computeSystemEnergyEwald
    
    state = MCState()
    state.info.box = [5.0, 5.0, 5.0]
    state.info.cutoff = 2.0
    
    ff = MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.ljEps = [0.0]
    ff.ljSigma = [0.3]
    state.forcefield = ff
    
    # Create a small system
    positions = [[2.0, 2.0, 2.5], [3.0, 2.0, 2.5], [2.5, 3.0, 2.5]]
    charges = [1.0, -1.0, 0.5]
    
    atoms = []
    for i in range(3):
        atom = MCAtom()
        atom.x, atom.y, atom.z = positions[i]
        atom.charge = charges[i]
        atom.type = 0
        atoms.append(atom)
    
    state.atoms = atoms
    state.activeAtomCount = 3
    
    residues = []
    for i in range(3):
        res = MCResidue()
        res.active = True
        res.fixed = False
        res.atomStart = i
        res.atomCount = 1
        res.type = 0
        residues.append(res)
    
    state.residues = residues
    state.activeResidueCount = 3
    
    # Compute PME
    initializePMEParameters(state.info.cutoff, state.info.box, 2.5)
    computeSystemEnergyPME(state)
    pme_total = state.ewald_energy.get('total', 0.0)
    
    # Compute Ewald
    initializeEwaldParameters(state.info.cutoff, state.info.box, 2.5)
    result = computeSystemEnergyEwald(state)
    ewald_total = result[0]  # First element of tuple is total electrostatic
    
    # They should be very close (within 1% for electrostatic part)
    assert abs(pme_total - ewald_total) / abs(ewald_total) < 0.01, \
        f"PME total {pme_total} differs from Ewald {ewald_total} by more than 1%"


if __name__ == "__main__":
    test_pme_total_single_residue()
    print("✓ Single residue test passed")
    
    test_pme_total_multiple_residues()
    print("✓ Multiple residue test passed")
    
    test_pme_total_vs_ewald()
    print("✓ PME vs Ewald consistency test passed")
    
    print("\nAll tests passed!")
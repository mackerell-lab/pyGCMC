"""
Test PME Total field fix

Verify that PME total energy calculation is correct after bug fix
"""

import pytest
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pygcmc
from pygcmc import MCState, MCAtom, MCResidue, MCForceField
from pygcmc import initializePMEParameters, computeSystemEnergyPME, initializeEwaldParameters, computeSystemEnergyEwald


@pytest.fixture(autouse=True)
def clear_pme_state():
    """Clear PME engine state before each test to prevent contamination"""
    try:
        pygcmc.clearPMEEngine()
    except:
        pass
    yield
    try:
        pygcmc.clearPMEEngine()
    except:
        pass


def test_pme_total_single_residue():
    """Test PME total for single residue system (should be correct)"""
    # Create state
    state = MCState()
    state.info.box = [5.0, 5.0, 5.0]
    state.info.cutoff = 2.0
    
    ff = MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.ljEps = [0.0]
    ff.ljSigma = [0.3]
    state.forcefield = ff
    
    # Create atoms
    atoms = []
    atom1 = MCAtom()
    atom1.x, atom1.y, atom1.z = 2.0, 2.5, 2.5
    atom1.charge = 1.0
    atom1.type = 0
    atoms.append(atom1)
    
    atom2 = MCAtom()
    atom2.x, atom2.y, atom2.z = 3.0, 2.5, 2.5
    atom2.charge = -1.0
    atom2.type = 0
    atoms.append(atom2)
    
    state.atoms = atoms
    state.activeAtomCount = len(atoms)
    
    # Create residue
    residues = []
    res = MCResidue()
    res.active = True
    res.fixed = False
    res.atomStart = 0
    res.atomCount = 2
    res.type = 0
    residues.append(res)
    
    state.residues = residues
    state.activeResidueCount = len(residues)
    
    # Run PME
    initializePMEParameters(state.info.cutoff, state.info.box, 2.5)
    computeSystemEnergyPME(state)
    
    # Check that total equals sum of components
    total = state.ewald_energy.get('total', 0.0)
    expected = state.ewald_energy.get('real_space', 0.0) + state.ewald_energy.get('reciprocal', 0.0) + state.ewald_energy.get('self', 0.0)
    
    assert abs(total - expected) < 1e-6, f"PME total {total} != sum {expected}"


def test_pme_total_multiple_residues():
    """Test PME total for multiple residue system (bug was here)"""
    # Create state
    state = MCState()
    state.info.box = [5.0, 5.0, 5.0]
    state.info.cutoff = 2.0
    
    ff = MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.ljEps = [0.0]
    ff.ljSigma = [0.3]
    state.forcefield = ff
    
    # Create atoms
    atoms = []
    atom1 = MCAtom()
    atom1.x, atom1.y, atom1.z = 2.0, 2.5, 2.5
    atom1.charge = 1.0
    atom1.type = 0
    atoms.append(atom1)
    
    atom2 = MCAtom()
    atom2.x, atom2.y, atom2.z = 3.0, 2.5, 2.5
    atom2.charge = -1.0
    atom2.type = 0
    atoms.append(atom2)
    
    state.atoms = atoms
    state.activeAtomCount = len(atoms)
    
    # Create residues (each atom in separate residue)
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
    state.activeResidueCount = len(residues)
    
    # Run PME
    initializePMEParameters(state.info.cutoff, state.info.box, 2.5)
    computeSystemEnergyPME(state)
    
    # Check that total equals sum of components
    total = state.ewald_energy.get('total', 0.0)
    expected = state.ewald_energy.get('real_space', 0.0) + state.ewald_energy.get('reciprocal', 0.0) + state.ewald_energy.get('self', 0.0)
    
    assert abs(total - expected) < 1e-6, f"PME total {total} != sum {expected}"
    
    # Verify no huge offset (bug would give ~38000 kJ/mol offset)
    assert abs(total) < 1000, f"PME total {total} has unreasonable magnitude"


def test_pme_total_vs_ewald():
    """Test that PME and Ewald give consistent total energies"""
    # Create state for PME
    state_pme = MCState()
    state_pme.info.box = [5.0, 5.0, 5.0]
    state_pme.info.cutoff = 2.0
    
    ff = MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.ljEps = [0.0]
    ff.ljSigma = [0.3]
    state_pme.forcefield = ff
    
    # Create atoms
    positions = [[2.0, 2.0, 2.5], [3.0, 2.0, 2.5], [2.5, 3.0, 2.5]]
    charges = [1.0, -1.0, 0.5]
    
    atoms_pme = []
    for i in range(3):
        atom = MCAtom()
        atom.x, atom.y, atom.z = positions[i]
        atom.charge = charges[i]
        atom.type = 0
        atoms_pme.append(atom)
    
    state_pme.atoms = atoms_pme
    state_pme.activeAtomCount = len(atoms_pme)
    
    # Create residues (each atom in separate residue)
    residues_pme = []
    for i in range(3):
        res = MCResidue()
        res.active = True
        res.fixed = False
        res.atomStart = i
        res.atomCount = 1
        res.type = 0
        residues_pme.append(res)
    
    state_pme.residues = residues_pme
    state_pme.activeResidueCount = len(residues_pme)
    
    # Clear PME state and run PME
    try:
        pygcmc.clearPMEEngine()
    except:
        pass
    
    initializePMEParameters(state_pme.info.cutoff, state_pme.info.box, 2.5)
    computeSystemEnergyPME(state_pme)
    pme_total = state_pme.ewald_energy.get('total', 0.0)
    
    # Create independent state for Ewald
    state_ewald = MCState()
    state_ewald.info.box = [5.0, 5.0, 5.0]
    state_ewald.info.cutoff = 2.0
    state_ewald.forcefield = ff
    
    atoms_ewald = []
    for i in range(3):
        atom = MCAtom()
        atom.x, atom.y, atom.z = positions[i]
        atom.charge = charges[i]
        atom.type = 0
        atoms_ewald.append(atom)
    
    state_ewald.atoms = atoms_ewald
    state_ewald.activeAtomCount = len(atoms_ewald)
    
    residues_ewald = []
    for i in range(3):
        res = MCResidue()
        res.active = True
        res.fixed = False
        res.atomStart = i
        res.atomCount = 1
        res.type = 0
        residues_ewald.append(res)
    
    state_ewald.residues = residues_ewald
    state_ewald.activeResidueCount = len(residues_ewald)
    
    # Run Ewald
    initializeEwaldParameters(state_ewald.info.cutoff, state_ewald.info.box, 2.5)
    ewald_result = computeSystemEnergyEwald(state_ewald)
    ewald_total = ewald_result[0]
    
    # They should be very close (within 1% for electrostatic part)
    if abs(ewald_total) > 1e-6:  # Avoid division by zero
        relative_diff = abs(pme_total - ewald_total) / abs(ewald_total)
        assert relative_diff < 0.01, \
            f"PME total {pme_total} differs from Ewald {ewald_total} by {relative_diff*100:.2f}%"
    else:
        assert abs(pme_total - ewald_total) < 1e-6, \
            f"PME total {pme_total} differs from Ewald {ewald_total}"


# The test_pme_total_vs_ewald function is now defined above


if __name__ == "__main__":
    test_pme_total_single_residue()
    print("✓ Single residue test passed")
    
    test_pme_total_multiple_residues()
    print("✓ Multiple residue test passed")
    
    test_pme_total_vs_ewald()
    print("✓ PME vs Ewald consistency test passed")
    
    print("\nAll tests passed!")
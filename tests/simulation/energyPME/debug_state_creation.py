"""
Debug tests for MCState creation and atom manipulation
"""
import pytest
import pygcmc
from pygcmc import MCState, MCAtom, MCResidue, MCForceField

def test_state_creation():
    """Test MCState creation and basic operations"""
    state = MCState()
    
    # Check initial state
    assert state.activeAtomCount == 0
    assert state.activeResidueCount == 0
    assert len(state.atoms) == 0
    assert len(state.residues) == 0
    
    # Initialize arrays
    state.atoms = []
    state.residues = []
    assert len(state.atoms) == 0
    assert len(state.residues) == 0

def test_atom_addition_methods():
    """Test different methods of adding atoms to MCState"""
    state = MCState()
    
    # Create atom
    atom = MCAtom()
    atom.x = 1.5
    atom.y = 1.5
    atom.z = 1.5
    atom.charge = 1.0
    atom.type = 0
    
    # Verify atom properties
    assert atom.x == 1.5
    assert atom.y == 1.5
    assert atom.z == 1.5
    assert atom.charge == 1.0
    assert atom.type == 0
    
    # The correct way to add atoms: create a list and assign it
    atoms = [atom]
    state.atoms = atoms
    state.activeAtomCount = len(atoms)
    
    # Verify the atom was added
    assert len(state.atoms) == 1
    assert state.activeAtomCount == 1
    assert state.atoms[0].x == 1.5
    assert state.atoms[0].charge == 1.0
    
    # Test adding multiple atoms
    atom2 = MCAtom()
    atom2.x = 2.5
    atom2.y = 2.5
    atom2.z = 2.5
    atom2.charge = -1.0
    atom2.type = 1
    
    atoms = [atom, atom2]
    state.atoms = atoms
    state.activeAtomCount = len(atoms)
    
    assert len(state.atoms) == 2
    assert state.activeAtomCount == 2
    assert state.atoms[1].x == 2.5
    assert state.atoms[1].charge == -1.0

def test_state_with_residues():
    """Test MCState with atoms and residues"""
    state = MCState()
    state.atoms = []
    state.residues = []
    
    # Create and add atoms
    atoms = []
    for i in range(3):
        atom = MCAtom()
        atom.x = atom.y = atom.z = float(i)
        atom.charge = 1.0 if i % 2 == 0 else -1.0
        atom.type = i % 2
        atoms.append(atom)
    
    state.atoms = atoms
    state.activeAtomCount = len(atoms)
    
    # Create residue
    res = MCResidue()
    res.atomStart = 0
    res.atomCount = 3
    res.active = True
    res.type = 0
    
    state.residues = [res]
    state.activeResidueCount = 1
    
    # Verify state
    assert state.activeAtomCount == 3
    assert state.activeResidueCount == 1
    assert len(state.atoms) == 3
    assert len(state.residues) == 1
    assert state.residues[0].atomCount == 3
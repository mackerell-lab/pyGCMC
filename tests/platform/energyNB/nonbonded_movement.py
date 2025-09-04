# tests/simulation/energyNB/nonbonded_movement.py
"""Non-bonded movement and validation tests."""

import pytest
import pygcmc
import math


def test_three_movement_molecules():
    """Test nonbonded energy calculation for three movement molecules.
    
    Setup:
    - Three movement residues (all type 0)
    - Each residue has one atom with zero charge (only vdw interaction)
    - All residues are active
    - Force field parameters:
        eps = 1.0, sigma = 1.0
        r = 1.0 (distance between atoms)
    
    Expected:
    - First movement residue interacts with two other residues
    - Each interaction at r = sigma gives V = 0
    - Total energy for first residue = 0
    """
    # Create a state object
    state = pygcmc.MCState()
    
    # 1. Set up force field
    state.forcefield.numTotalTypes = 1  # Only one type (movement type 0)
    state.forcefield.numMovementTypes = 1  # One movement type
    
    # Initialize force field parameters
    # For movement type 0:
    #   - interaction with type 0: index = 0 * 1 + 0 = 0
    state.forcefield.ljEps = [1.0]     # eps = 1.0 for movement-movement interaction
    state.forcefield.ljSigma = [1.0]   # sigma = 1.0 for movement-movement interaction
    
    # Set up movement atom types
    state.movementAtomTypes = [0]  # Type 0 is a movement type
    state.numMovementAtomTypes = 1
    
    # 2. Set up atoms (place them at unit distance from each other)
    atoms = []
    positions = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0)]  # Form a right triangle with unit distances
    
    for x, y, z in positions:
        atom = pygcmc.MCAtom()
        atom.x = x
        atom.y = y
        atom.z = z
        atom.charge = 0.0  # No electrostatic interaction
        atom.type = 0  # All atoms are movement type
        atoms.append(atom)
    
    state.atoms = atoms
    state.activeAtomCount = 3
    
    # 3. Set up three movement residues
    residues = []
    for i in range(3):
        res = pygcmc.MCResidue()
        res.active = True
        res.type = 0  # All are movement type
        res.atomStart = i
        res.atomCount = 1
        residues.append(res)
    
    state.residues = residues
    state.activeResidueCount = 3
    
    # 4. Set up movement residue info
    movement_info = pygcmc.MCMovementResidueInfo()
    movement_info.startIndex = 0
    movement_info.activeCount = 3  # All three residues are active
    movement_info.totalCount = 3
    
    state.movementResidues = [movement_info]
    
    # Calculate energy
    pygcmc.computeMovementEnergy(state)
    energy = state.residues[0].energy_vdw + state.residues[0].energy_elec
    
    # Check result
    assert abs(energy) < 1e-5, "Expected zero energy at r = sigma"


def test_invalid_forcefield_params():
    """Test handling of invalid force field parameters."""
    state = pygcmc.MCState()
    
    # Set up force field with incorrect parameter array size
    state.forcefield.numTotalTypes = 2
    state.forcefield.numMovementTypes = 1
    state.forcefield.ljEps = [1.0]  # Should be size 2
    state.forcefield.ljSigma = [1.0, 1.0]
    
    # Set up minimal state
    atom = pygcmc.MCAtom()
    state.atoms = [atom]
    res = pygcmc.MCResidue()
    state.residues = [res]
    movement_info = pygcmc.MCMovementResidueInfo()
    state.movementResidues = [movement_info]
    
    # Expect runtime error due to invalid parameter array size
    with pytest.raises(RuntimeError):
        pygcmc.computeMovementEnergy(state)


def test_inactive_residue():
    """Test energy calculation with inactive residues.
    
    Setup:
    - Two residues, but second one is inactive
    - Should result in zero energy since no interaction is computed
    """
    state = pygcmc.MCState()
    
    # 1. Set up force field
    state.forcefield.numTotalTypes = 2
    state.forcefield.numMovementTypes = 1
    state.forcefield.ljEps = [1.0, 1.0, 1.0, 1.0]    # Complete 2x2 matrix for eps
    state.forcefield.ljSigma = [1.0, 1.0, 1.0, 1.0]  # Complete 2x2 matrix for sigma
    
    # 2. Set up atoms
    atom1 = pygcmc.MCAtom()
    atom1.x = 0.0
    atom1.y = 0.0
    atom1.z = 0.0
    atom1.charge = 1.0
    atom1.type = 0
    
    atom2 = pygcmc.MCAtom()
    atom2.x = 1.0
    atom2.y = 0.0
    atom2.z = 0.0
    atom2.charge = -1.0
    atom2.type = 1
    
    state.atoms = [atom1, atom2]
    state.activeAtomCount = 2
    
    # 3. Set up residues (second one inactive)
    movement_res = pygcmc.MCResidue()
    movement_res.active = True
    movement_res.type = 0
    movement_res.atomStart = 0
    movement_res.atomCount = 1
    
    fixed_res = pygcmc.MCResidue()
    fixed_res.active = False  # Inactive
    fixed_res.type = 1
    fixed_res.atomStart = 1
    fixed_res.atomCount = 1
    
    state.residues = [movement_res, fixed_res]
    state.activeResidueCount = 1  # Only one active residue
    
    # Set up movementAtomTypes
    state.movementAtomTypes = [0]  # Type 0 is a movement type
    state.numMovementAtomTypes = 1
    
    # 4. Set up movement residue info
    movement_info = pygcmc.MCMovementResidueInfo()
    movement_info.startIndex = 0
    movement_info.activeCount = 1
    movement_info.totalCount = 1
    movement_info.resName = "MOV"  # Add residue name
    
    state.movementResidues = [movement_info]
    
    # Calculate energy
    pygcmc.computeMovementEnergy(state)
    energy = state.residues[0].energy_vdw + state.residues[0].energy_elec
    
    # Check result
    assert abs(energy) < 1e-6, "Expected zero energy with inactive residue"
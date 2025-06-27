# tests/simulation/energyNB/nonbonded_system.py
"""Non-bonded system-wide tests."""

import pytest
import pygcmc
import math


def test_all_residues_nonbonded():
    """Test nonbonded energy calculation for all residues in the system.
    
    Setup:
    - Multiple residues with various interactions
    - Mix of movement and fixed residues
    - Both VDW and electrostatic interactions
    """
    state = pygcmc.MCState()
    
    # 1. Set up force field
    state.forcefield.numTotalTypes = 3  # Three types
    state.forcefield.numMovementTypes = 3
    
    # Set up force field parameters (3x3 matrix flattened)
    state.forcefield.ljEps = [1.0] * 9    # All interactions eps=1.0
    state.forcefield.ljSigma = [1.0] * 9  # All interactions sigma=1.0
    
    # 2. Set up three atoms forming an equilateral triangle
    atoms = []
    # Place atoms at the vertices of an equilateral triangle
    positions = [
        (0.0, 0.0, 0.0),           # Origin
        (1.0, 0.0, 0.0),           # At 1nm on x-axis
        (0.5, 0.866, 0.0)          # Complete the equilateral triangle
    ]
    charges = [0.5, -0.5, 0.5]     # Alternating charges
    
    for i, (x, y, z) in enumerate(positions):
        atom = pygcmc.MCAtom()
        atom.x = x
        atom.y = y
        atom.z = z
        atom.charge = charges[i]
        atom.type = i
        atoms.append(atom)
    
    state.atoms = atoms
    state.activeAtomCount = 3
    
    # 3. Set up three residues
    residues = []
    for i in range(3):
        res = pygcmc.MCResidue()
        res.active = True
        res.type = i
        res.atomStart = i
        res.atomCount = 1
        residues.append(res)
    
    state.residues = residues
    state.activeResidueCount = 3
    
    # Calculate energy
    pygcmc.computeSystemEnergy(state)
    
    # Verify each residue has energy
    for i in range(3):
        energy = state.residues[i].energy_vdw + state.residues[i].energy_elec
        assert not math.isnan(energy), f"Residue {i} has NaN energy"
        assert not math.isinf(energy), f"Residue {i} has infinite energy"
    
    # Calculate total energy
    total_vdw = sum(res.energy_vdw for res in state.residues)
    total_elec = sum(res.energy_elec for res in state.residues)
    total_energy = total_vdw + total_elec
    
    # Total energy needs to be divided by 2 because each interaction was calculated twice
    system_energy = total_energy / 2.0
    
    # Verify the reasonableness of total energy
    assert system_energy < 1e6, "Total system energy too large"
    assert system_energy > -1e6, "Total system energy too negative"
    
    # Verify that energy magnitude is reasonable
    # For our setup:
    # - LJ energy is 0 at r=sigma
    # - Electrostatic energy = k_c * q1*q2/r
    # For three residues, there are three pairs of interactions:
    # 1. residue 0-1: q1=0.5, q2=-0.5, r=1.0
    # 2. residue 1-2: q1=-0.5, q2=0.5, r=1.0
    # 3. residue 2-0: q1=0.5, q2=0.5, r=1.0
    # Each pair of interactions is calculated once and added to the two related residues
    COULOMB = 138.935458  # kJ·nm/mol/e²
    
    # Calculate energy for each interaction pair
    e_01 = COULOMB * 0.5 * (-0.5) / 1.0  # residue 0-1 interaction
    e_12 = COULOMB * (-0.5) * 0.5 / 1.0  # residue 1-2 interaction
    e_20 = COULOMB * 0.5 * 0.5 / 1.0     # residue 2-0 interaction
    
    # Total system energy is the sum of all interactions (each interaction has already been added twice in the code)
    expected_system_energy = (e_01 + e_12 + e_20)
    
    rel_tol = 0.1  # 10% relative error tolerance
    assert abs((system_energy - expected_system_energy) / expected_system_energy) < rel_tol, \
           f"System energy {system_energy} differs too much from expected {expected_system_energy}"


def test_all_residues_inactive():
    """Test all residues nonbonded energy with inactive residues."""
    state = pygcmc.MCState()
    
    # Basic setup
    state.forcefield.numTotalTypes = 2
    state.forcefield.numMovementTypes = 2
    state.forcefield.ljEps = [1.0] * 4
    state.forcefield.ljSigma = [1.0] * 4
    
    # Set up two atoms
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
    
    # Set up two inactive residues
    res1 = pygcmc.MCResidue()
    res1.active = False
    res1.type = 0
    res1.atomStart = 0
    res1.atomCount = 1
    
    res2 = pygcmc.MCResidue()
    res2.active = False
    res2.type = 1
    res2.atomStart = 1
    res2.atomCount = 1
    
    state.residues = [res1, res2]
    state.activeResidueCount = 0
    
    # Calculate energy
    pygcmc.computeSystemEnergy(state)
    
    # Verify all energies are 0
    for res in state.residues:
        assert abs(res.energy_vdw) < 1e-6, "Inactive residue has non-zero VDW energy"
        assert abs(res.energy_elec) < 1e-6, "Inactive residue has non-zero electrostatic energy"
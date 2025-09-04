# tests/simulation/energyNB/nonbonded_cutoff.py
"""Non-bonded cutoff tests."""

import pytest
import pygcmc
import math


def test_cutoff_nonperiodic():
    """Test nonbonded energy calculation with distance cutoff.
    
    Setup:
    - Multiple residues at various distances
    - Only interactions within cutoff distance are considered
    - No periodic boundary conditions
    """
    state = pygcmc.MCState()
    
    # 1. Set up force field
    state.forcefield.numTotalTypes = 3  # Three types
    state.forcefield.numMovementTypes = 3
    
    # Set up force field parameters (3x3 matrix flattened)
    state.forcefield.ljEps = [1.0] * 9    # All interactions eps=1.0
    state.forcefield.ljSigma = [1.0] * 9  # All interactions sigma=1.0
    
    # Set cutoff distance
    state.info.cutoff = 1.5  # nm
    
    # 2. Set up three atoms arranged sequentially on the x-axis
    atoms = []
    positions = [
        (0.0, 0.0, 0.0),  # Origin
        (1.0, 0.0, 0.0),  # x=1.0 nm, within cutoff distance
        (2.0, 0.0, 0.0)   # x=2.0 nm, beyond cutoff distance
    ]
    charges = [0.5, -0.5, 0.5]  # Alternating charges
    
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
    pygcmc.computeSystemEnergyCutoff(state)
    
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
    # Only residue 0-1 interaction is within cutoff distance:
    # - q1=0.5, q2=-0.5, r=1.0
    # - The energy of each interaction is fully added to both residues
    COULOMB = 138.935458  # kJ·nm/mol/e²
    expected_energy = 2.0 * COULOMB * 0.5 * (-0.5) / 1.0  # Interaction energy between one pair of atoms, multiplied by 2 because it's added twice in the code

    rel_tol = 0.1  # 10% relative error tolerance
    assert abs((system_energy - expected_energy) / expected_energy) < rel_tol, \
           f"System energy {system_energy} differs too much from expected {expected_energy}"
    
    # Verify that residue 2 energy is 0 (because its distance to other residues exceeds cutoff)
    res2_energy = state.residues[2].energy_vdw + state.residues[2].energy_elec
    # residue 2 should have energy because it receives energy from interaction with residue 1
    expected_res2_energy = COULOMB * 0.5 * (-0.5) / 1.0  # residue 0-1 interaction
    assert abs((res2_energy - expected_res2_energy) / expected_res2_energy) < rel_tol, \
           f"Residue 2 energy {res2_energy} differs too much from expected {expected_res2_energy}"
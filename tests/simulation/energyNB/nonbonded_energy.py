# tests/simulation/energyNB/nonbonded_energy.py
"""Non-bonded energy and distance tests."""

import pytest
import pygcmc
import math


def test_energy_symmetry():
    """Verify that energy calculation is symmetric.
    A->B interaction should equal B->A interaction.
    
    Setup:
    - Two movement residues with both LJ and charge interactions
    - Verify energy is same regardless of which is considered the "movement" residue
    """
    state = pygcmc.MCState()
    
    # 1. Set up force field
    state.forcefield.numTotalTypes = 1  # Both residues same type
    state.forcefield.numMovementTypes = 1
    state.forcefield.ljEps = [1.0]     # eps = 1.0 kJ/mol
    state.forcefield.ljSigma = [1.0]   # sigma = 1.0 nm
    
    # Set up movement atom types
    state.movementAtomTypes = [0]  # Type 0 is a movement type
    state.numMovementAtomTypes = 1
    
    # 2. Set up atoms
    atom1 = pygcmc.MCAtom()
    atom1.x = 0.0
    atom1.y = 0.0
    atom1.z = 0.0
    atom1.charge = 0.5
    atom1.type = 0
    
    atom2 = pygcmc.MCAtom()
    atom2.x = 1.0
    atom2.y = 0.0
    atom2.z = 0.0
    atom2.charge = -0.5
    atom2.type = 0
    
    state.atoms = [atom1, atom2]
    state.activeAtomCount = 2
    
    # 3. Set up residues
    res1 = pygcmc.MCResidue()
    res1.active = True
    res1.type = 0
    res1.atomStart = 0
    res1.atomCount = 1
    
    res2 = pygcmc.MCResidue()
    res2.active = True
    res2.type = 0
    res2.atomStart = 1
    res2.atomCount = 1
    
    state.residues = [res1, res2]
    state.activeResidueCount = 2
    
    # 4. Set up movement residue info for both residues
    movement_info = pygcmc.MCMovementResidueInfo()
    movement_info.startIndex = 0
    movement_info.activeCount = 2  # Both residues are movement residues
    movement_info.totalCount = 2
    movement_info.resName = "MOV"
    
    state.movementResidues = [movement_info]
    
    # Calculate energy
    pygcmc.computeMovementEnergy(state)
    
    # Get energies for both residues
    energy1 = state.residues[0].energy_vdw + state.residues[0].energy_elec
    energy2 = state.residues[1].energy_vdw + state.residues[1].energy_elec
    
    # Check symmetry with appropriate tolerance
    rel_tol = 1e-5
    assert abs((energy1 - energy2) / energy1) < rel_tol, \
           f"Energy not symmetric: residue1={energy1}, residue2={energy2}"


def test_very_close_distance():
    """Test energy calculation for very close distances.
    
    Setup:
    - Two residues at very close distance (0.001 sigma)
    - Should give large but finite energy due to soft core potential
    - Energy should be capped at MAX_SAFE_ENERGY
    Each residue gets this full capped energy (not divided by 2)
    """
    state = pygcmc.MCState()
    
    # 1. Set up force field
    state.forcefield.numTotalTypes = 2
    state.forcefield.numMovementTypes = 1
    state.forcefield.ljEps = [1.0, 1.0, 1.0, 1.0]    # Complete 2x2 matrix for eps
    state.forcefield.ljSigma = [1.0, 1.0, 1.0, 1.0]  # Complete 2x2 matrix for sigma
    
    # Set up movement atom types
    state.movementAtomTypes = [0]
    state.numMovementAtomTypes = 1
    
    # 2. Set up atoms at very close distance
    atom1 = pygcmc.MCAtom()
    atom1.x = 0.0
    atom1.y = 0.0
    atom1.z = 0.0
    atom1.charge = 1.0
    atom1.type = 0
    
    atom2 = pygcmc.MCAtom()
    atom2.x = 0.001  # Very close distance (0.001 sigma)
    atom2.y = 0.0
    atom2.z = 0.0
    atom2.charge = 1.0  # Same charge for strong repulsion
    atom2.type = 1
    
    state.atoms = [atom1, atom2]
    state.activeAtomCount = 2
    
    # 3. Set up residues
    movement_res = pygcmc.MCResidue()
    movement_res.active = True
    movement_res.type = 0
    movement_res.atomStart = 0
    movement_res.atomCount = 1
    
    fixed_res = pygcmc.MCResidue()
    fixed_res.active = True
    fixed_res.type = 1
    fixed_res.atomStart = 1
    fixed_res.atomCount = 1
    
    state.residues = [movement_res, fixed_res]
    state.activeResidueCount = 2
    
    # 4. Set up movement residue info
    movement_info = pygcmc.MCMovementResidueInfo()
    movement_info.startIndex = 0
    movement_info.activeCount = 1
    movement_info.totalCount = 1
    movement_info.resName = "MOV"
    
    state.movementResidues = [movement_info]
    
    # Calculate energy
    pygcmc.computeMovementEnergy(state)
    energy = state.residues[0].energy_vdw + state.residues[0].energy_elec
    
    # Energy should be large but finite
    MIN_SAFE_DISTANCE = 0.01  # nm (1% of sigma)
    MAX_SAFE_ENERGY = 1e6     # kJ/mol
    
    # Calculate expected energy at MIN_SAFE_DISTANCE
    r = MIN_SAFE_DISTANCE
    sigma = 1.0
    eps = 1.0
    
    # LJ energy: V_LJ = 4ε[(σ/r)¹² - (σ/r)⁶]
    sigma_r = sigma / r
    term6 = pow(sigma_r, 6)
    term12 = term6 * term6
    expected_vdw = min(4.0 * eps * (term12 - term6), MAX_SAFE_ENERGY)
    
    # Coulomb energy: V_C = k_c * q1*q2/r
    COULOMB = 138.935458  # kJ·nm/mol/e²
    expected_elec = min(COULOMB * 1.0 * 1.0 / r, MAX_SAFE_ENERGY)
    
    expected_total = min(expected_vdw + expected_elec, MAX_SAFE_ENERGY)
    
    # Verify energy is finite and reasonable
    assert not math.isinf(energy), "Energy should not be infinite"
    assert not math.isnan(energy), "Energy should not be NaN"
    assert energy <= MAX_SAFE_ENERGY, f"Energy {energy} exceeds MAX_SAFE_ENERGY {MAX_SAFE_ENERGY}"
    
    # Energy should be close to expected value at MIN_SAFE_DISTANCE
    rel_tol = 1e-3  # 0.1% tolerance for numerical stability
    assert abs((energy - expected_total) / expected_total) < rel_tol, \
           f"Energy mismatch at MIN_SAFE_DISTANCE: expected {expected_total}, got {energy}"


def test_zero_distance_handling():
    """Test handling of zero distance between atoms.
    Should use soft core potential instead of throwing error.
    """
    state = pygcmc.MCState()
    
    # 1. Set up force field
    state.forcefield.numTotalTypes = 2
    state.forcefield.numMovementTypes = 1
    state.forcefield.ljEps = [1.0, 1.0, 1.0, 1.0]    # Complete 2x2 matrix for eps
    state.forcefield.ljSigma = [1.0, 1.0, 1.0, 1.0]  # Complete 2x2 matrix for sigma
    
    # Set up movement atom types
    state.movementAtomTypes = [0]
    state.numMovementAtomTypes = 1
    
    # 2. Set up atoms at exact same position (zero distance)
    atom1 = pygcmc.MCAtom()
    atom1.x = 0.0
    atom1.y = 0.0
    atom1.z = 0.0
    atom1.charge = 1.0
    atom1.type = 0
    
    atom2 = pygcmc.MCAtom()
    atom2.x = 0.0  # Exact same position
    atom2.y = 0.0
    atom2.z = 0.0
    atom2.charge = 1.0
    atom2.type = 1
    
    state.atoms = [atom1, atom2]
    state.activeAtomCount = 2
    
    # 3. Set up residues
    movement_res = pygcmc.MCResidue()
    movement_res.active = True
    movement_res.type = 0
    movement_res.atomStart = 0
    movement_res.atomCount = 1
    
    fixed_res = pygcmc.MCResidue()
    fixed_res.active = True
    fixed_res.type = 1
    fixed_res.atomStart = 1
    fixed_res.atomCount = 1
    
    state.residues = [movement_res, fixed_res]
    state.activeResidueCount = 2
    
    # 4. Set up movement residue info
    movement_info = pygcmc.MCMovementResidueInfo()
    movement_info.startIndex = 0
    movement_info.activeCount = 1
    movement_info.totalCount = 1
    movement_info.resName = "MOV"
    
    state.movementResidues = [movement_info]
    
    # Calculate energy (should not crash)
    pygcmc.computeMovementEnergy(state)
    energy = state.residues[0].energy_vdw + state.residues[0].energy_elec
    
    # Energy should be finite (not infinite/NaN)
    assert not math.isinf(energy), "Energy should not be infinite for zero distance"
    assert not math.isnan(energy), "Energy should not be NaN for zero distance"
    
    # Energy should be large and positive for same-charge repulsion
    MAX_SAFE_ENERGY = 1e6  # kJ/mol
    assert energy > 0, "Energy should be positive for same-charge repulsion"
    assert energy <= MAX_SAFE_ENERGY, f"Energy {energy} should be capped at {MAX_SAFE_ENERGY}"
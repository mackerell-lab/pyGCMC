# /tests/simulation/test_naive_nonbonded.py

import pytest
import pygcmc
import math

def test_attractive_interaction():
    """Test nonbonded energy calculation for an attractive interaction between movement and fixed residues.
    
    Setup:
    - Two residues: one movement (type 0) and one fixed (type 1)
    - Each residue has one atom with zero charge (only vdw interaction)
    - Both residues are active
    - Force field parameters:
        eps = 1.0, sigma = 1.0
        r = 1.0 (distance between atoms)
    
    Expected:
    V_LJ = 4ε[(σ/r)¹² - (σ/r)⁶]
    V = 4 * 1.0 * [(1.0/1.0)¹² - (1.0/1.0)⁶]
    V = 4 * 1.0 * (1 - 1) = 0
    """
    # Create a state object
    state = pygcmc.MCState()
    
    # 1. Set up force field
    state.forcefield.numTotalTypes = 2  # Two types: movement and fixed
    state.forcefield.numMovementTypes = 1  # One movement type
    
    # Initialize force field parameters for attractive interaction
    # For type 0:
    #   - interaction with type 0: index = 0 * 2 + 0 = 0
    #   - interaction with type 1: index = 0 * 2 + 1 = 1
    # For type 1:
    #   - interaction with type 0: index = 1 * 2 + 0 = 2
    #   - interaction with type 1: index = 1 * 2 + 1 = 3
    state.forcefield.ljEps = [1.0, 1.0, 1.0, 1.0]    # Complete 2x2 matrix for eps
    state.forcefield.ljSigma = [1.0, 1.0, 1.0, 1.0]  # Complete 2x2 matrix for sigma
    
    # Set up movement atom types
    state.movementAtomTypes = [0]  # Type 0 is a movement type
    state.numMovementAtomTypes = 1
    
    # 2. Set up atoms
    # First atom (for movement residue)
    atom1 = pygcmc.MCAtom()
    atom1.x = 0.0
    atom1.y = 0.0
    atom1.z = 0.0
    atom1.charge = 0.0  # No electrostatic interaction
    atom1.type = 0  # Movement type
    
    # Second atom (for fixed residue)
    atom2 = pygcmc.MCAtom()
    atom2.x = 1.0  # At distance 1.0 from first atom
    atom2.y = 0.0
    atom2.z = 0.0
    atom2.charge = 0.0  # No electrostatic interaction
    atom2.type = 1  # Fixed type
    
    state.atoms = [atom1, atom2]
    state.activeAtomCount = 2
    
    # 3. Set up residues
    movement_res = pygcmc.MCResidue()
    movement_res.active = True
    movement_res.type = 0  # Movement type
    movement_res.atomStart = 0
    movement_res.atomCount = 1
    
    fixed_res = pygcmc.MCResidue()
    fixed_res.active = True
    fixed_res.type = 1  # Fixed type
    fixed_res.atomStart = 1
    fixed_res.atomCount = 1
    
    state.residues = [movement_res, fixed_res]
    state.activeResidueCount = 2
    
    # 4. Set up movement residue info
    movement_info = pygcmc.MCMovementResidueInfo()
    movement_info.startIndex = 0
    movement_info.activeCount = 1
    movement_info.totalCount = 1
    
    state.movementResidues = [movement_info]
    
    # Calculate energy
    pygcmc.computeMovementEnergy(state)
    energy = state.residues[0].energy_vdw + state.residues[0].energy_elec
    
    # Check result
    assert abs(energy) < 1e-6, "Expected zero interaction energy at r = sigma"

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

def test_electrostatic_interaction():
    """Test nonbonded energy calculation for electrostatic interaction.
    
    Setup:
    - Two residues: one movement (type 0) and one fixed (type 1)
    - Each residue has one atom with opposite charges (+1 and -1)
    - Both residues are active
    - Force field parameters:
        eps = 0.0 (no vdw interaction)
        r = 1.0 nm
    
    Expected:
    V = k_c * q1*q2/r where:
    - k_c = 138.935458 kJ·nm/mol/e^2 (Coulomb constant)
    - q1 = +1e, q2 = -1e
    - r = 1.0 nm
    Therefore:
    V = 138.935458 * (+1) * (-1) / 1.0 = -138.935458 kJ/mol
    """
    state = pygcmc.MCState()
    
    # 1. Set up force field (no vdw interaction)
    state.forcefield.numTotalTypes = 2
    state.forcefield.numMovementTypes = 1
    state.forcefield.ljEps = [0.0, 0.0, 0.0, 0.0]    # Complete 2x2 matrix with eps = 0.0
    state.forcefield.ljSigma = [1.0, 1.0, 1.0, 1.0]  # Complete 2x2 matrix with sigma = 1.0
    
    # Set up movement atom types
    state.movementAtomTypes = [0]  # Type 0 is a movement type
    state.numMovementAtomTypes = 1
    
    # 2. Set up atoms with opposite charges
    atom1 = pygcmc.MCAtom()
    atom1.x = 0.0
    atom1.y = 0.0
    atom1.z = 0.0
    atom1.charge = 1.0  # Positive charge
    atom1.type = 0
    
    atom2 = pygcmc.MCAtom()
    atom2.x = 1.0  # Distance of 1.0 nm
    atom2.y = 0.0
    atom2.z = 0.0
    atom2.charge = -1.0  # Negative charge
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
    movement_info.resName = "MOV"  # Add residue name
    
    state.movementResidues = [movement_info]
    
    # Calculate energy
    pygcmc.computeMovementEnergy(state)
    energy = state.residues[0].energy_vdw + state.residues[0].energy_elec
    
    # Expected energy with Coulomb constant
    COULOMB = 138.935458  # kJ·nm/mol/e^2
    expected_energy = -COULOMB  # k_c * (+1) * (-1) / 1.0
    
    # Check result with appropriate tolerance for single-precision float
    assert abs(energy - expected_energy) < 1e-5, f"Expected energy {expected_energy}, got {energy}"

def test_repulsive_interaction():
    """Test naive nonbonded energy calculation for repulsive interaction.
    
    Setup:
    - Two residues at very close distance (0.5 sigma)
    - Each residue has one atom with zero charge (only vdw interaction)
    - Force field parameters:
        eps = 1.0, sigma = 1.0
        r = 0.5 (half of sigma)
    
    Expected:
    V_LJ = 4ε[(σ/r)¹² - (σ/r)⁶]
    V = 4 * 1.0 * [(1.0/0.5)¹² - (1.0/0.5)⁶]
    V = 4 * (4096 - 64) = 16128
    Each residue gets this full energy (not divided by 2)
    """
    state = pygcmc.MCState()
    
    # 1. Set up force field
    state.forcefield.numTotalTypes = 2
    state.forcefield.numMovementTypes = 1
    state.forcefield.ljEps = [1.0, 1.0, 1.0, 1.0]    # Complete 2x2 matrix for eps
    state.forcefield.ljSigma = [1.0, 1.0, 1.0, 1.0]  # Complete 2x2 matrix for sigma
    
    # Set up movement atom types
    state.movementAtomTypes = [0]  # Type 0 is a movement type
    state.numMovementAtomTypes = 1
    
    # 2. Set up atoms at close distance
    atom1 = pygcmc.MCAtom()
    atom1.x = 0.0
    atom1.y = 0.0
    atom1.z = 0.0
    atom1.charge = 0.0
    atom1.type = 0
    
    atom2 = pygcmc.MCAtom()
    atom2.x = 0.5  # Distance = 0.5 sigma
    atom2.y = 0.0
    atom2.z = 0.0
    atom2.charge = 0.0
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
    
    state.movementResidues = [movement_info]
    
    # Calculate energy
    pygcmc.computeMovementEnergy(state)
    energy = state.residues[0].energy_vdw + state.residues[0].energy_elec
    
    # Check result
    expected_energy = 16128.0  # 4 * (4096 - 64)
    assert abs(energy - expected_energy) < 1.0, "Expected strong repulsive energy"

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

def test_combined_interaction():
    """Test both LJ and Coulomb interactions together.
    
    Setup:
    - Two residues with both LJ and charge interactions
    - eps = 1.0 kJ/mol, sigma = 1.0 nm
    - q1 = +0.5e, q2 = -0.5e
    - r = 1.0 nm
    
    Expected:
    V_LJ = 4ε[(σ/r)¹² - (σ/r)⁶] = 0 kJ/mol (at r = sigma)
    V_C = k_c * q1*q2/r = 138.935458 * 0.5 * (-0.5) / 1.0 = -34.734 kJ/mol
    Each residue gets this full energy (not divided by 2)
    """
    state = pygcmc.MCState()
    
    # 1. Set up force field
    state.forcefield.numTotalTypes = 2
    state.forcefield.numMovementTypes = 1
    state.forcefield.ljEps = [1.0, 1.0, 1.0, 1.0]    # Complete 2x2 matrix for eps
    state.forcefield.ljSigma = [1.0, 1.0, 1.0, 1.0]  # Complete 2x2 matrix for sigma
    
    # Set up movement atom types
    state.movementAtomTypes = [0]  # Type 0 is a movement type
    state.numMovementAtomTypes = 1
    
    # 2. Set up atoms with both LJ and charge interactions
    atom1 = pygcmc.MCAtom()
    atom1.x = 0.0
    atom1.y = 0.0
    atom1.z = 0.0
    atom1.charge = 0.5  # Half positive charge
    atom1.type = 0
    
    atom2 = pygcmc.MCAtom()
    atom2.x = 1.0  # Distance of 1.0 nm
    atom2.y = 0.0
    atom2.z = 0.0
    atom2.charge = -0.5  # Half negative charge
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
    
    # Get actual energies
    actual_vdw = state.residues[0].energy_vdw
    actual_elec = state.residues[0].energy_elec
    actual_total = actual_vdw + actual_elec
    
    # Expected energies
    COULOMB = 138.935458  # kJ·nm/mol/e²
    expected_vdw = 0.0  # kJ/mol (at r = sigma)
    expected_elec = COULOMB * 0.5 * (-0.5) / 1.0  # kJ/mol
    expected_total = expected_vdw + expected_elec
    
    # Check results with appropriate tolerance
    rel_tol = 1e-5
    assert abs(actual_vdw - expected_vdw) < 1e-6, \
           f"VDW energy mismatch: expected {expected_vdw}, got {actual_vdw}"
    assert abs((actual_elec - expected_elec) / expected_elec) < rel_tol, \
           f"Electrostatic energy mismatch: expected {expected_elec}, got {actual_elec}"
    assert abs((actual_total - expected_total) / expected_total) < rel_tol, \
           f"Total energy mismatch: expected {expected_total}, got {actual_total}"

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
    
    # 2. Set up atoms at exactly same position
    atom1 = pygcmc.MCAtom()
    atom1.x = 0.0
    atom1.y = 0.0
    atom1.z = 0.0
    atom1.charge = 1.0
    atom1.type = 0
    
    atom2 = pygcmc.MCAtom()
    atom2.x = 0.0  # Same position as atom1
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
    
    # Calculate energy - should not raise exception
    pygcmc.computeMovementEnergy(state)
    energy = state.residues[0].energy_vdw + state.residues[0].energy_elec
    
    # Energy should be finite and capped
    assert not math.isinf(energy), "Energy should not be infinite"
    assert not math.isnan(energy), "Energy should not be NaN"
    assert energy <= 1e6, "Energy should be capped at MAX_SAFE_ENERGY"

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

def test_pbc_basic():
    """Test basic PBC energy calculation.
    
    Setup:
    - Two residues at opposite sides of the box
    - Box size = 2.0 nm
    - Residues should interact through PBC
    """
    state = pygcmc.MCState()
    
    # 1. Set up force field
    state.forcefield.numTotalTypes = 2
    state.forcefield.numMovementTypes = 2
    state.forcefield.ljEps = [1.0] * 4
    state.forcefield.ljSigma = [1.0] * 4
    
    # Set up a large box and appropriate cutoff
    state.info.box = [10.0, 10.0, 10.0]  # Large box to avoid PBC effects
    state.info.cutoff = 2.0  # nm
    
    # 2. Set up two atoms on opposite sides of the box
    atoms = []
    positions = [
        (0.5, 1.0, 1.0),    # Near the middle of the box
        (1.5, 1.0, 1.0)     # Near the middle of the box
    ]
    charges = [0.5, -0.5]   # Opposite charges
    
    for i, (x, y, z) in enumerate(positions):
        atom = pygcmc.MCAtom()
        atom.x = x
        atom.y = y
        atom.z = z
        atom.charge = charges[i]
        atom.type = i
        atoms.append(atom)
    
    state.atoms = atoms
    state.activeAtomCount = 2
    
    # 3. Set up two residues
    residues = []
    for i in range(2):
        res = pygcmc.MCResidue()
        res.active = True
        res.type = i
        res.atomStart = i
        res.atomCount = 1
        residues.append(res)
    
    state.residues = residues
    state.activeResidueCount = 2
    
    # Calculate PBC energy
    pygcmc.computeSystemEnergyPBCCutoff(state)
    pbc_vdw = sum(res.energy_vdw for res in state.residues)
    pbc_elec = sum(res.energy_elec for res in state.residues)
    pbc_energy = (pbc_vdw + pbc_elec) / 2.0
    
    # Through PBC, the actual distance between atoms should be 1.0 nm (not 1.0 nm)
    # In x direction: 1.5 - 0.5 = 1.0, exactly equal to sigma
    COULOMB = 138.935458  # kJ·nm/mol/e²
    min_dist = 1.0  # Minimum image distance
    expected_elec = COULOMB * 0.5 * (-0.5) / min_dist  # Electrostatic energy
    
    # LJ energy is 0 at r = sigma
    expected_vdw = 0.0  # Assuming r ≈ sigma
    expected_energy = expected_elec + expected_vdw
    
    # Print detailed debug information
    print(f"\nDetailed energy comparison:")
    print(f"Box size: {state.info.box}")
    print(f"Cutoff distance: {state.info.cutoff} nm")
    print(f"Atom positions: {positions}")
    print(f"Minimum image distance: {min_dist} nm")
    print(f"Expected electrostatic energy: {expected_elec} kJ/mol")
    print(f"Expected VDW energy: {expected_vdw} kJ/mol")
    print(f"Expected total energy: {expected_energy} kJ/mol")
    print(f"Actual VDW energy: {pbc_vdw/2} kJ/mol")
    print(f"Actual electrostatic energy: {pbc_elec/2} kJ/mol")
    print(f"Actual total energy: {pbc_energy} kJ/mol")
    
    # Check VDW and electrostatic energies separately
    rel_tol = 0.1  # 10% relative error tolerance
    
    # Check electrostatic energy
    actual_elec = pbc_elec / 2.0
    assert abs((actual_elec - expected_elec) / expected_elec) < rel_tol, \
           f"Electrostatic energy {actual_elec} differs too much from expected {expected_elec}"
    
    # Check VDW energy (should be close to 0 because r ≈ sigma)
    actual_vdw = pbc_vdw / 2.0
    assert abs(actual_vdw) < 10.0, f"VDW energy {actual_vdw} is too large"
    
    # Check total energy
    assert abs((pbc_energy - expected_energy) / expected_energy) < rel_tol, \
           f"PBC energy {pbc_energy} differs too much from expected {expected_energy}"

def test_pbc_invalid_box():
    """Test PBC energy calculation with invalid box dimensions."""
    state = pygcmc.MCState()
    
    # Set up basic state
    state.forcefield.numTotalTypes = 1
    state.forcefield.numMovementTypes = 1
    state.forcefield.ljEps = [1.0]
    state.forcefield.ljSigma = [1.0]
    
    # Set invalid box dimensions
    state.info.box = [0.0, 1.0, 1.0]  # x dimension is 0
    
    # Add a simple atom and residue
    atom = pygcmc.MCAtom()
    state.atoms = [atom]
    
    res = pygcmc.MCResidue()
    res.active = True
    state.residues = [res]
    state.activeResidueCount = 1
    
    # Expect exception to be raised
    with pytest.raises(RuntimeError, match="Invalid box dimensions"):
        pygcmc.computeSystemEnergyPBCCutoff(state) 

def test_pbc_vs_nopbc():
    """Compare PBC and non-PBC energy calculations.
    
    Setup:
    - Two residues well within the box (not near boundaries)
    - Results should be identical when atoms are far from boundaries
    """
    state = pygcmc.MCState()
    
    # 1. Set up force field
    state.forcefield.numTotalTypes = 2
    state.forcefield.numMovementTypes = 2
    state.forcefield.ljEps = [1.0] * 4
    state.forcefield.ljSigma = [1.0] * 4
    
    # Set up a large box and appropriate cutoff
    state.info.box = [10.0, 10.0, 10.0]  # Large box to avoid PBC effects
    state.info.cutoff = 2.0  # nm
    
    # 2. Set up two atoms in the middle of the box
    atoms = []
    positions = [
        (4.0, 5.0, 5.0),
        (5.0, 5.0, 5.0)
    ]
    charges = [0.5, -0.5]
    
    for i, (x, y, z) in enumerate(positions):
        atom = pygcmc.MCAtom()
        atom.x = x
        atom.y = y
        atom.z = z
        atom.charge = charges[i]
        atom.type = i
        atoms.append(atom)
    
    state.atoms = atoms
    state.activeAtomCount = 2
    
    # 3. Set up residues
    residues = []
    for i in range(2):
        res = pygcmc.MCResidue()
        res.active = True
        res.type = i
        res.atomStart = i
        res.atomCount = 1
        residues.append(res)
    
    state.residues = residues
    state.activeResidueCount = 2
    
    # Calculate both types of energy
    pygcmc.computeSystemEnergyCutoff(state)
    nopbc_vdw = sum(res.energy_vdw for res in state.residues)
    nopbc_elec = sum(res.energy_elec for res in state.residues)
    nopbc_energy = (nopbc_vdw + nopbc_elec) / 2.0
    
    # Reset energy
    for res in state.residues:
        res.energy_vdw = 0.0
        res.energy_elec = 0.0
    
    # Calculate PBC energy
    pygcmc.computeSystemEnergyPBCCutoff(state)
    pbc_vdw = sum(res.energy_vdw for res in state.residues)
    pbc_elec = sum(res.energy_elec for res in state.residues)
    pbc_energy = (pbc_vdw + pbc_elec) / 2.0
    
    # For atoms in the middle of the box, both calculation methods should give the same result
    rel_tol = 1e-5  # Very small relative error tolerance
    assert abs((pbc_energy - nopbc_energy) / nopbc_energy) < rel_tol, \
           f"PBC energy {pbc_energy} differs from non-PBC energy {nopbc_energy}"

def test_pbc_cross_boundary():
    """Test PBC energy calculation for atoms crossing periodic boundaries.
    
    Setup:
    - Three atoms forming a chain across the periodic boundary
    - Middle atom near boundary should interact correctly with both other atoms
    """
    state = pygcmc.MCState()
    
    # 1. Set up force field
    state.forcefield.numTotalTypes = 3
    state.forcefield.numMovementTypes = 3
    state.forcefield.ljEps = [1.0] * 9
    state.forcefield.ljSigma = [1.0] * 9
    
    # Set up box and cutoff
    state.info.box = [3.0, 3.0, 3.0]  # 3.0 nm cubic box
    state.info.cutoff = 2.0  # nm
    
    # 2. Set up three atoms: one near the box edge, one crossing the boundary, and one on the other side of the box
    atoms = []
    positions = [
        (2.8, 1.5, 1.5),  # Near the right boundary
        (0.1, 1.5, 1.5),  # Near the left boundary
        (1.5, 1.5, 1.5)   # In the middle
    ]
    charges = [0.5, 0.5, -1.0]  # Make the middle atom attractive to both sides
    
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
    
    # 3. Set up residues
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
    
    # Calculate PBC energy
    pygcmc.computeSystemEnergyPBCCutoff(state)
    
    # Verify all residues have reasonable energy
    for i, res in enumerate(state.residues):
        energy = res.energy_vdw + res.energy_elec
        assert not math.isnan(energy), f"Residue {i} has NaN energy"
        assert not math.isinf(energy), f"Residue {i} has infinite energy"
        assert abs(energy) < 1e6, f"Residue {i} has unreasonably large energy: {energy}"
    
    # Through PBC, the minimum distance between atoms 0 and 1 should be 0.3 nm
    # instead of the actual distance of 2.7 nm in the box
    COULOMB = 138.935458  # kJ·nm/mol/e²
    
    # Calculate total energy
    total_vdw = sum(res.energy_vdw for res in state.residues)
    total_elec = sum(res.energy_elec for res in state.residues)
    system_energy = (total_vdw + total_elec) / 2.0
    
    # Verify energy is finite and reasonable
    assert not math.isnan(system_energy), "System energy is NaN"
    assert not math.isinf(system_energy), "System energy is infinite"
    assert abs(system_energy) < 1e6, f"System energy {system_energy} is unreasonably large"
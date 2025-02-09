import pytest
import pygcmc
import math

def test_attractive_interaction():
    """Test naive nonbonded energy calculation for an attractive interaction.
    
    Setup:
    - Two residues: one movement (type 0) and one fixed (type 1)
    - Each residue has one atom with zero charge (only vdw interaction)
    - Both residues are active
    - Force field parameters:
        eps = 1.0, sigma = 1.0
        r = 1.0 (fixed in naive implementation)
    
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
    # For movement type 0:
    #   - interaction with type 0: index = 0 * 2 + 0 = 0
    #   - interaction with type 1: index = 0 * 2 + 1 = 1
    state.forcefield.ljEps = [1.0, 1.0]    # eps = 1.0 for both interactions
    state.forcefield.ljSigma = [1.0, 1.0]  # sigma = 1.0 for both interactions
    
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
    pygcmc.computeNaiveNonbondedEnergy(state)
    energy = state.residues[0].energy_vdw + state.residues[0].energy_elec
    
    # Check result
    assert abs(energy) < 1e-6, "Expected zero interaction energy at r = sigma"

def test_three_movement_molecules():
    """Test naive nonbonded energy calculation for three movement molecules.
    
    Setup:
    - Three movement residues (all type 0)
    - Each residue has one atom with zero charge (only vdw interaction)
    - All residues are active
    - Force field parameters:
        eps = 1.0, sigma = 1.0
        r = 1.0 (fixed in naive implementation)
    
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
    pygcmc.computeNaiveNonbondedEnergy(state)
    energy = state.residues[0].energy_vdw + state.residues[0].energy_elec
    
    # Check result
    assert abs(energy) < 1e-5, "Expected zero energy at r = sigma"

def test_electrostatic_interaction():
    """Test naive nonbonded energy calculation for electrostatic interaction.
    
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
    Each residue gets this full energy (not divided by 2)
    """
    state = pygcmc.MCState()
    
    # 1. Set up force field (no vdw interaction)
    state.forcefield.numTotalTypes = 2
    state.forcefield.numMovementTypes = 1
    state.forcefield.ljEps = [0.0, 0.0]    # eps = 0.0 to disable vdw
    state.forcefield.ljSigma = [1.0, 1.0]  # sigma doesn't matter when eps = 0
    
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
    pygcmc.computeNaiveNonbondedEnergy(state)
    energy = state.residues[0].energy_vdw + state.residues[0].energy_elec
    
    # Expected energy with Coulomb constant
    COULOMB = 138.935458  # kJ·nm/mol/e^2
    expected_energy = -COULOMB  # k_c * (+1) * (-1) / 1.0
    
    # Check result with appropriate tolerance for single-precision float
    rel_tol = 1e-5  # 0.001% relative tolerance
    abs_diff = abs(energy - expected_energy)
    rel_diff = abs_diff / abs(expected_energy)
    assert rel_diff < rel_tol, f"Expected electrostatic energy of {expected_energy} kJ/mol, got {energy} kJ/mol (relative error: {rel_diff})"

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
    state.forcefield.ljEps = [1.0, 1.0]
    state.forcefield.ljSigma = [1.0, 1.0]
    
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
    pygcmc.computeNaiveNonbondedEnergy(state)
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
        pygcmc.computeNaiveNonbondedEnergy(state)

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
    state.forcefield.ljEps = [1.0, 1.0]
    state.forcefield.ljSigma = [1.0, 1.0]
    
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
    pygcmc.computeNaiveNonbondedEnergy(state)
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
    state.forcefield.ljEps = [1.0, 1.0]    # eps = 1.0 kJ/mol
    state.forcefield.ljSigma = [1.0, 1.0]  # sigma = 1.0 nm
    
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
    pygcmc.computeNaiveNonbondedEnergy(state)
    
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
    pygcmc.computeNaiveNonbondedEnergy(state)
    
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
    state.forcefield.ljEps = [1.0, 1.0]
    state.forcefield.ljSigma = [1.0, 1.0]
    
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
    pygcmc.computeNaiveNonbondedEnergy(state)
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
    state.forcefield.ljEps = [1.0, 1.0]
    state.forcefield.ljSigma = [1.0, 1.0]
    
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
    pygcmc.computeNaiveNonbondedEnergy(state)
    energy = state.residues[0].energy_vdw + state.residues[0].energy_elec
    
    # Energy should be finite and capped
    assert not math.isinf(energy), "Energy should not be infinite"
    assert not math.isnan(energy), "Energy should not be NaN"
    assert energy <= 1e6, "Energy should be capped at MAX_SAFE_ENERGY"

def test_all_residues_nonbonded():
    """Test nonbonded energy calculation for all residues.
    
    Setup:
    - Three active residues in a triangular configuration
    - Each residue has one atom with both charge and LJ interactions
    - All residues should have correct pairwise interactions
    
    Note:
    - Each pairwise interaction is calculated once and added to both residues
    - Total system energy should be sum of all residue energies divided by 2
      (because each interaction is counted twice in the sum)
    """
    state = pygcmc.MCState()
    
    # 1. 设置力场
    state.forcefield.numTotalTypes = 3  # 三种类型
    state.forcefield.numMovementTypes = 3
    
    # 设置力场参数 (3x3矩阵展平)
    state.forcefield.ljEps = [1.0] * 9    # 所有相互作用eps=1.0
    state.forcefield.ljSigma = [1.0] * 9  # 所有相互作用sigma=1.0
    
    # 2. 设置三个原子，形成等边三角形
    atoms = []
    # 在等边三角形顶点上放置原子
    positions = [
        (0.0, 0.0, 0.0),           # 原点
        (1.0, 0.0, 0.0),           # x轴上1nm处
        (0.5, 0.866, 0.0)          # 完成等边三角形
    ]
    charges = [0.5, -0.5, 0.5]     # 交替电荷
    
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
    
    # 3. 设置三个residues
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
    
    # 计算能量
    pygcmc.computeAllNonbondedEnergy(state)
    
    # 验证每个residue都有能量
    for i in range(3):
        energy = state.residues[i].energy_vdw + state.residues[i].energy_elec
        assert not math.isnan(energy), f"Residue {i} has NaN energy"
        assert not math.isinf(energy), f"Residue {i} has infinite energy"
    
    # 计算总能量
    total_vdw = sum(res.energy_vdw for res in state.residues)
    total_elec = sum(res.energy_elec for res in state.residues)
    total_energy = total_vdw + total_elec
    
    # 总能量除以2（因为每个相互作用被计算了两次）
    system_energy = total_energy / 2.0
    
    # 验证总能量合理性
    assert system_energy < 1e6, "Total system energy too large"
    assert system_energy > -1e6, "Total system energy too negative"
    
    # 验证能量的数量级合理性
    # 对于我们的设置：
    # - LJ能量在r=sigma时为0
    # - 静电能量 = k_c * q1*q2/r
    # 对于三个residue，有三对相互作用：
    # 1. residue 0-1: q1=0.5, q2=-0.5, r=1.0
    # 2. residue 1-2: q1=-0.5, q2=0.5, r=1.0
    # 3. residue 2-0: q1=0.5, q2=0.5, r=1.0
    # 每对相互作用计算一次，加到两个相关的residue上
    COULOMB = 138.935458  # kJ·nm/mol/e²
    
    # 计算每对相互作用的能量
    e_01 = COULOMB * 0.5 * (-0.5) / 1.0  # residue 0-1 相互作用
    e_12 = COULOMB * (-0.5) * 0.5 / 1.0  # residue 1-2 相互作用
    e_20 = COULOMB * 0.5 * 0.5 / 1.0     # residue 2-0 相互作用
    
    # 总系统能量是所有相互作用的和
    expected_system_energy = e_01 + e_12 + e_20
    
    rel_tol = 0.1  # 10%的相对误差容忍度
    assert abs((system_energy - expected_system_energy) / expected_system_energy) < rel_tol, \
           f"System energy {system_energy} differs too much from expected {expected_system_energy}"

def test_all_residues_inactive():
    """Test all residues nonbonded energy with inactive residues."""
    state = pygcmc.MCState()
    
    # 基本设置
    state.forcefield.numTotalTypes = 2
    state.forcefield.numMovementTypes = 2
    state.forcefield.ljEps = [1.0] * 4
    state.forcefield.ljSigma = [1.0] * 4
    
    # 设置两个原子
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
    
    # 设置两个inactive residues
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
    
    # 计算能量
    pygcmc.computeAllNonbondedEnergy(state)
    
    # 验证所有能量都是0
    for res in state.residues:
        assert abs(res.energy_vdw) < 1e-6, "Inactive residue has non-zero VDW energy"
        assert abs(res.energy_elec) < 1e-6, "Inactive residue has non-zero electrostatic energy" 
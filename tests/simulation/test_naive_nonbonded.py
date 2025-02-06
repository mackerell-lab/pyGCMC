import pytest
import pygcmc

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
    V = eps * [(sigma/r)^12 - 2*(sigma/r)^6]
    V = 1.0 * (1.0 - 2.0) = -1.0
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
    assert abs(energy - (-1.0)) < 1e-6, "Expected attractive interaction energy of -1.0"

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
    - Each interaction contributes V = -1.0 (vdw only)
    - Total energy for first residue = 2 * (-1.0) = -2.0
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
    # First residue interacts with two other residues at unit distance
    # Each interaction contributes -1.0 energy
    expected_energy = -2.0  # 2 * (-1.0)
    assert abs(energy - expected_energy) < 1e-5, f"Expected energy of first residue to be {expected_energy}"

def test_electrostatic_interaction():
    """Test naive nonbonded energy calculation for electrostatic interaction.
    
    Setup:
    - Two residues: one movement (type 0) and one fixed (type 1)
    - Each residue has one atom with opposite charges (+1 and -1)
    - Both residues are active
    - Force field parameters:
        eps = 0.0 (no vdw interaction)
        r = 1.0 (fixed in naive implementation)
    
    Expected:
    V = q1*q2/r = (+1)*(-1)/1.0 = -1.0
    """
    state = pygcmc.MCState()
    
    # 1. Set up force field (no vdw interaction)
    state.forcefield.numTotalTypes = 2
    state.forcefield.numMovementTypes = 1
    state.forcefield.ljEps = [0.0, 0.0]    # eps = 0.0 to disable vdw
    state.forcefield.ljSigma = [1.0, 1.0]  # sigma doesn't matter when eps = 0
    
    # 2. Set up atoms with opposite charges
    atom1 = pygcmc.MCAtom()
    atom1.x = 0.0
    atom1.y = 0.0
    atom1.z = 0.0
    atom1.charge = 1.0  # Positive charge
    atom1.type = 0
    
    atom2 = pygcmc.MCAtom()
    atom2.x = 1.0
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
    
    state.movementResidues = [movement_info]
    
    # Calculate energy
    pygcmc.computeNaiveNonbondedEnergy(state)
    energy = state.residues[0].energy_vdw + state.residues[0].energy_elec
    
    # Check result
    assert abs(energy - (-1.0)) < 1e-6, "Expected electrostatic energy of -1.0"

def test_repulsive_interaction():
    """Test naive nonbonded energy calculation for repulsive interaction.
    
    Setup:
    - Two residues at very close distance (0.5 sigma)
    - Each residue has one atom with zero charge (only vdw interaction)
    - Force field parameters:
        eps = 1.0, sigma = 1.0
        r = 0.5 (half of sigma)
    
    Expected:
    V = eps * [(sigma/r)^12 - 2*(sigma/r)^6]
    V = 1.0 * [(1.0/0.5)^12 - 2*(1.0/0.5)^6]
    V = 1.0 * (4096 - 128) = 3968
    """
    state = pygcmc.MCState()
    
    # 1. Set up force field
    state.forcefield.numTotalTypes = 2
    state.forcefield.numMovementTypes = 1
    state.forcefield.ljEps = [1.0, 1.0]
    state.forcefield.ljSigma = [1.0, 1.0]
    
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
    expected_energy = 3968.0  # 1.0 * (4096 - 128)
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
    assert abs(energy) < 1e-6, "Expected zero energy with inactive residue" 
import pytest
import pygcmc

def test_attractive_interaction():
    """Test naive nonbonded energy calculation for an attractive interaction.
    
    Setup:
    - Two residues: one movement (type 0) and one fixed (type 1)
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
    state.forcefield.maxTypes = 2  # Two types: movement and fixed
    state.forcefield.numMovementTypes = 1  # One movement type
    
    # Initialize force field parameters for attractive interaction
    # For movement type 0:
    #   - interaction with type 0: index = 0 * 2 + 0 = 0
    #   - interaction with type 1: index = 0 * 2 + 1 = 1
    state.forcefield.ljEps = [1.0, 1.0]    # eps = 1.0 for both interactions
    state.forcefield.ljSigma = [1.0, 1.0]  # sigma = 1.0 for both interactions
    
    # 2. Set up residues
    movement_res = pygcmc.MCResidue()
    movement_res.active = True
    movement_res.type = 0  # Movement type
    
    fixed_res = pygcmc.MCResidue()
    fixed_res.active = True
    fixed_res.type = 1  # Fixed type
    
    state.residues = [movement_res, fixed_res]
    
    # 3. Set up movement residue info
    movement_info = pygcmc.MCMovementResidueInfo()
    movement_info.startIndex = 0
    movement_info.activeCount = 1
    movement_info.totalCount = 1
    
    state.movementResidues = [movement_info]
    
    # Calculate energy
    energy = pygcmc.computeNaiveNonbondedEnergy(state)
    
    # Check result
    assert abs(energy - (-1.0)) < 1e-6, "Expected attractive interaction energy of -1.0"

def test_three_movement_molecules():
    """Test naive nonbonded energy calculation for three movement molecules.
    
    Setup:
    - Three movement residues (all type 0)
    - All residues are active
    - Force field parameters:
        eps = 1.0, sigma = 1.0
        r = 1.0 (fixed in naive implementation)
    
    Expected:
    - Each pair of movement residues contributes V = -1.0
    - Total number of pairs = 3 choose 2 = 3 pairs
    - Total energy = 3 * (-1.0) = -3.0
    """
    # Create a state object
    state = pygcmc.MCState()
    
    # 1. Set up force field
    state.forcefield.maxTypes = 1  # Only one type (movement type 0)
    state.forcefield.numMovementTypes = 1  # One movement type
    
    # Initialize force field parameters
    # For movement type 0:
    #   - interaction with type 0: index = 0 * 1 + 0 = 0
    state.forcefield.ljEps = [1.0]     # eps = 1.0 for movement-movement interaction
    state.forcefield.ljSigma = [1.0]   # sigma = 1.0 for movement-movement interaction
    
    # 2. Set up three movement residues
    residues = []
    for _ in range(3):
        res = pygcmc.MCResidue()
        res.active = True
        res.type = 0  # All are movement type
        residues.append(res)
    
    state.residues = residues
    
    # 3. Set up movement residue info
    movement_info = pygcmc.MCMovementResidueInfo()
    movement_info.startIndex = 0
    movement_info.activeCount = 3  # All three residues are active
    movement_info.totalCount = 3
    
    state.movementResidues = [movement_info]
    
    # Calculate energy
    energy = pygcmc.computeNaiveNonbondedEnergy(state)
    
    # Check result
    # We expect 3 pairs of interactions:
    # - residue 0 with residue 1
    # - residue 0 with residue 2
    # - residue 1 with residue 2
    # Each pair contributes -1.0 to the total energy
    expected_energy = -3.0  # 3 pairs * (-1.0)
    assert abs(energy - expected_energy) < 1e-6, f"Expected total energy of {expected_energy}" 
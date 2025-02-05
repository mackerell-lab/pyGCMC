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
    state.forcefield.ljEps = [1.0]    # eps = 1.0 for the interaction
    state.forcefield.ljSigma = [1.0]  # sigma = 1.0 for the interaction
    
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
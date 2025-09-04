"""
Pytest configuration and shared fixtures for movementCPP tests.

This file is automatically loaded by pytest and provides fixtures
available to all test files in this directory and subdirectories.
"""

import pytest
import pygcmc


@pytest.fixture
def setup_system_with_params():
    """Setup test system with standard parameters.
    
    Returns a tuple of (state, params).
    """
    state = pygcmc.MCState()
    state.info.box = [5.0, 5.0, 5.0]
    
    # Setup force field with known interactions
    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.ljEps = [0.5]  # kJ/mol
    ff.ljSigma = [0.3]  # nm
    state.forcefield = ff
    
    params = pygcmc.movement.MovementParams()
    params.temperature = 298.15  # K
    params.chemicalPotential = -15.7  # kJ/mol
    params.seed = 42
    
    return state, params


@pytest.fixture
def setup_system():
    """Setup test system - returns only state for compatibility.
    
    Returns a state with 5x5x5 nm box and standard force field parameters.
    """
    state = pygcmc.MCState()
    state.info.box = [5.0, 5.0, 5.0]
    
    # Setup force field with known interactions
    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 2
    ff.numMovementTypes = 2
    ff.ljEps = [0.65, 0.0, 0.0, 0.0]  # kJ/mol
    ff.ljSigma = [0.3165, 0.0, 0.0, 0.0]  # nm
    state.forcefield = ff
    
    return state


@pytest.fixture
def clean_system():
    """Create a fresh system for each test."""
    state = pygcmc.MCState()
    state.info.box = [5.0, 5.0, 5.0]
    
    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 2
    ff.numMovementTypes = 2
    ff.ljEps = [0.65, 0.0, 0.0, 0.0]
    ff.ljSigma = [0.3165, 0.0, 0.0, 0.0]
    state.forcefield = ff
    
    return state


@pytest.fixture
def fresh_system():
    """Create a fresh system for statistics tests."""
    state = pygcmc.MCState()
    state.info.box = [5.0, 5.0, 5.0]
    
    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 2
    ff.numMovementTypes = 2
    ff.ljEps = [0.65, 0.0, 0.0, 0.0]
    ff.ljSigma = [0.3165, 0.0, 0.0, 0.0]
    state.forcefield = ff
    
    return state


@pytest.fixture
def create_state():
    """Create a test MCState with given box size."""
    def _create(box_nm=5.0):
        state = pygcmc.MCState()
        if isinstance(box_nm, (list, tuple)):
            state.info.box = list(box_nm)
        else:
            state.info.box = [box_nm, box_nm, box_nm]
        
        # Setup force field
        ff = pygcmc.MCForceField()
        ff.numTotalTypes = 1
        ff.numMovementTypes = 1
        ff.ljEps = [0.5]
        ff.ljSigma = [0.3]
        state.forcefield = ff
        
        return state
    return _create
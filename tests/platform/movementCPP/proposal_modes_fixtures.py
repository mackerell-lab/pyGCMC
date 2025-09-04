# tests/simulation/movementCPP/proposal_modes_fixtures.py
"""Fixtures for proposal mode tests."""

import pytest
import pygcmc
import numpy as np


@pytest.fixture
def setup_system():
    """Setup a basic system for testing - returns only state for compatibility."""
    state = pygcmc.MCState()
    state.info.box = np.array([5.0, 5.0, 5.0])
    
    # Setup force field
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
    state.info.box = np.array([5.0, 5.0, 5.0])
    
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
    state.info.box = np.array([5.0, 5.0, 5.0])
    
    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 2
    ff.numMovementTypes = 2
    ff.ljEps = [0.65, 0.0, 0.0, 0.0]
    ff.ljSigma = [0.3165, 0.0, 0.0, 0.0]
    state.forcefield = ff
    
    return state
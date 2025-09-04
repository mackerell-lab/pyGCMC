# tests/simulation/movementCPP/multi_insertion_fixtures.py
"""Multi-insertion CBMC test fixtures and utilities."""

import pytest
import pygcmc
import numpy as np
import concurrent.futures
import time
import math
    
# Fixture function
def setup_multi_insertion_system():
    """Setup test system for multi-insertion."""
    state = pygcmc.MCState()
    state.info.box = np.array([6.0, 6.0, 6.0])
    
    # Setup force field
    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.ljEps = [0.5]
    ff.ljSigma = [0.3]
    state.forcefield = ff
    
    params = pygcmc.movement.MovementParams()
    params.temperature = 298.15
    params.chemicalPotential = -15.0
    params.seed = 42
    
    return state, params
    

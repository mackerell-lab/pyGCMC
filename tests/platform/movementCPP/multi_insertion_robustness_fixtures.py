# tests/simulation/movementCPP/multi_insertion_robustness_fixtures.py
"""Multi-insertion CBMC robustness test fixtures and utilities."""

import pytest
import pygcmc
import numpy as np


def basic_robustness_system():
    """Basic system for robustness testing."""
    state = pygcmc.MCState()
    state.info.box = np.array([5.0, 5.0, 5.0])
    state.info.cutoff = 1.2  # nm
    
    # Simple force field
    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 2  # Water O and H
    ff.numMovementTypes = 2
    ff.ljEps = [0.6, 0.3,
                0.3, 0.1]
    ff.ljSigma = [0.315, 0.2,
                  0.2, 0.1]
    state.forcefield = ff
    
    params = pygcmc.movement.MovementParams()
    params.temperature = 298.15
    params.chemicalPotential = -10.0
    params.seed = 99999
    params.useMultiInsertionCBMC = True
    params.maxParallelInsertions = 2
    params.numConfigTrials = 5
    params.minRegionSeparationNm = 1.5
    
    return state, params


def make_test_state(box_size=5.0, cutoff=1.2):
    """Create a test state with specified box size."""
    state = pygcmc.MCState()
    state.info.box = np.array([box_size, box_size, box_size])
    state.info.cutoff = cutoff
    
    # Simple force field
    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.ljEps = [0.5]
    ff.ljSigma = [0.3]
    state.forcefield = ff
    
    return state
# tests/simulation/movementCPP/multi_insertion_improved_fixtures.py
"""Multi-insertion CBMC improved test fixtures and utilities."""

import pytest
import pygcmc
import numpy as np

# Import statistical utilities
from .statistical_utils import (
    poisson_diff_ok,
    ratio_CI_ok,
    no_drift,
    run_until_stable,
    effective_sample_size,
    batch_means_variance
)


def has_multi_insertion():
    """Check if multi-insertion CBMC is available."""
    try:
        mover = pygcmc.movement.MovementModule()
        return hasattr(mover, 'attemptMultiInsertionCBMC')
    except:
        return False


def setup_improved_system():
    """Setup test system for multi-insertion."""
    state = pygcmc.MCState()
    state.info.box = np.array([6.0, 6.0, 6.0])
    
    # Setup force field
    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 2  # Support O and H atoms
    ff.numMovementTypes = 2
    # Simple 2x2 LJ matrix for O and H atoms
    ff.ljEps = [0.5, 0.45,
                0.45, 0.2]
    ff.ljSigma = [0.30, 0.275,
                  0.275, 0.20]
    state.forcefield = ff
    
    params = pygcmc.movement.MovementParams()
    params.temperature = 298.15
    params.chemicalPotential = -5.0  # Less negative for better acceptance
    params.useCavityBias = False  # Disable for reproducibility
    params.seed = 42  # Fixed seed for reproducibility
    
    # Multi-insertion specific
    params.useMultiInsertionCBMC = True
    params.maxParallelInsertions = 4
    params.numConfigTrials = 10
    params.minRegionSeparationNm = 1.5  # nm
    params.multiDisplacementFraction = 0.5
    params.multiUseRegionVolume = False
    
    return state, params
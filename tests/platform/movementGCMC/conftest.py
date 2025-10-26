# tests/simulation/movementGCMC/conftest.py
"""
Fixtures for GCMC movement tests
"""
import pytest
import sys
import os
from acceptance_log_utils import (  # noqa: F401
    acceptance_statistics,
    compute_detailed_balance_ratio,
    count_by_move_and_species,
    count_by_species,
    filter_by_move,
    filter_by_species,
    match_insert_delete_pairs,
    read_jsonl,
)

# Add build directory to path
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(__file__)))), 'build'))
import pygcmc
import numpy as np


@pytest.fixture
def setup_gcmc_state():
    """Create a basic MCState for GCMC testing"""
    state = pygcmc.MCState()
    
    # Initialize info structure
    state.info = pygcmc.MCInfo()
    state.info.box = [3.0, 3.0, 3.0]  # 3x3x3 nm box
    state.info.setTemperature(300.0)  # 300K
    state.info.cutoff = 1.2  # 1.2 nm cutoff
    state.info.max_residues = 1000
    state.info.max_atoms = 3000
    state.info.volume = 27.0  # 3^3
    
    # Setup force field
    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 10
    ff.numMovementTypes = 1
    ff.maxTypes = 10
    
    # Initialize LJ parameters (10x10 matrix)
    n_types = ff.numTotalTypes
    ff.ljSigma = [0.0] * (n_types * n_types)
    ff.ljEps = [0.0] * (n_types * n_types)
    
    # TIP3P water parameters
    # O-O (type 0)
    ff.ljSigma[0] = 0.3151  # nm
    ff.ljEps[0] = 0.6364    # kJ/mol
    
    state.forcefield = ff
    
    # Initialize empty residue and atom lists
    state.residues = []
    state.atoms = []
    state.activeResidueCount = 0
    state.activeAtomCount = 0
    
    return state


@pytest.fixture
def gcmc_params():
    """Create GCMC movement parameters"""
    params = pygcmc.movement.MovementParams()
    params.temperature = 300.0  # K
    params.chemicalPotential = -15.7  # kJ/mol (water at 300K)
    params.maxTranslation = 0.1  # nm
    params.maxRotation = 0.2  # radians
    params.useCavityBias = True
    params.cavityGridSpacing = 0.2  # nm
    params.probeRadius = 0.14  # nm
    params.seed = 12345  # For reproducibility
    return params


@pytest.fixture
def water_template():
    """Create a water molecule template (TIP3P)"""
    # Note: This requires FragmentTemplate to be exposed in Python bindings
    # For now, return a dictionary representation
    template = {
        'name': 'WAT',
        'atoms': [
            {'name': 'O', 'type': 0, 'charge': -0.834, 'mass': 15.999,
             'x': 0.0, 'y': 0.0, 'z': 0.0},
            {'name': 'H1', 'type': 1, 'charge': 0.417, 'mass': 1.008,
             'x': 0.0957, 'y': 0.0, 'z': 0.0},
            {'name': 'H2', 'type': 1, 'charge': 0.417, 'mass': 1.008,
             'x': -0.024, 'y': 0.0927, 'z': 0.0}
        ],
        'chemical_potential': -15.7,  # kJ/mol
        'use_cavity_bias': True,
        'use_config_bias': False
    }
    return template


@pytest.fixture
def gcmc_mover(gcmc_params):
    """Create a configured GCMC movement module with proper seeding"""
    # Use MovementModule(params) for deterministic seeding
    mover = pygcmc.movement.MovementModule(gcmc_params)
    mover.resetStatistics()
    return mover


@pytest.fixture
def small_state():
    """Create a small box state for testing"""
    state = pygcmc.MCState()
    
    state.info = pygcmc.MCInfo()
    state.info.box = [2.0, 2.0, 2.0]  # 2x2x2 nm box
    state.info.setTemperature(300.0)
    state.info.cutoff = 0.9  # Smaller cutoff for small box
    state.info.max_residues = 500
    state.info.max_atoms = 1500
    state.info.volume = 8.0
    
    # Minimal force field setup
    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 10
    ff.numMovementTypes = 1
    ff.maxTypes = 10
    ff.ljSigma = [0.0] * 100
    ff.ljEps = [0.0] * 100
    ff.ljSigma[0] = 0.3151  # O-O
    ff.ljEps[0] = 0.6364
    
    state.forcefield = ff
    state.residues = []
    state.atoms = []
    state.activeResidueCount = 0
    state.activeAtomCount = 0
    
    return state


@pytest.fixture
def populated_state(setup_gcmc_state, gcmc_mover):
    """Create a state with some water molecules already inserted"""
    state = setup_gcmc_state
    
    # Attempt to insert a few water molecules
    for _ in range(10):
        gcmc_mover.attemptInsertion(state)
    
    return state


# Helper functions for tests
def calculate_acceptance_rate(mover):
    """Calculate overall acceptance rate from statistics"""
    stats = mover.getStatistics()
    if hasattr(stats, 'insert'):
        total_attempts = (stats.insert.attempts + stats.delete.attempts + 
                         stats.translate.attempts + stats.rotate.attempts)
        total_accepts = (stats.insert.accepts + stats.delete.accepts +
                        stats.translate.accepts + stats.rotate.accepts)
        if total_attempts > 0:
            return total_accepts / total_attempts
    return 0.0


def run_gcmc_steps(state, mover, n_steps=1000, seed=None):
    """Run a specified number of GCMC steps with optional seeding"""
    if seed is not None:
        np.random.seed(seed)
    
    results = []
    for _ in range(n_steps):
        # Choose move type randomly
        r = np.random.random()
        if r < 0.25:
            result = mover.attemptInsertion(state)
        elif r < 0.5 and state.activeResidueCount > 0:
            result = mover.attemptDeletion(state)
        elif r < 0.75 and state.activeResidueCount > 0:
            result = mover.attemptTranslation(state)
        elif state.activeResidueCount > 0:
            result = mover.attemptRotation(state)
        else:
            continue
        results.append(result)
    return results

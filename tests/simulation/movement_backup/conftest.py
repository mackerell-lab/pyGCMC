# tests/simulation/movement/conftest.py
"""Shared fixtures and utilities for movement tests."""

import pytest
import pygcmc
import numpy as np
import os
import random
import time


@pytest.fixture(autouse=True)
def deterministic_seed():
    """Set deterministic random seed for reproducibility.
    
    Each test gets a unique but deterministic seed based on:
    - Base seed from environment or default
    - Worker ID (for parallel workers)
    """
    import re
    base = int(os.getenv("PYTEST_SEED", "12345"))
    worker = os.getenv("PYTEST_XDIST_WORKER", "gw0")
    m = re.match(r"gw(\d+)", worker)
    wid = int(m.group(1)) if m else 0
    seed = (base + wid) % (2**32)
    
    # Set Python and NumPy seeds
    random.seed(seed)
    np.random.seed(seed)
    
    # Return seed for tests that need to set movement params
    yield seed


@pytest.fixture
def simple_system():
    """Create a simple test system with one atom type."""
    state = pygcmc.MCState()
    state.info.box = np.array([5.0, 5.0, 5.0])
    
    # Setup minimal force field
    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.ljEps = [0.5]
    ff.ljSigma = [0.3]
    state.forcefield = ff
    
    return state


@pytest.fixture
def movement_params(deterministic_seed):
    """Create default movement parameters with deterministic seed."""
    params = pygcmc.movement.MovementParams()
    params.temperature = 298.15
    params.chemicalPotential = -15.7
    params.seed = deterministic_seed  # Use fixture seed for reproducibility
    return params


@pytest.fixture
def configured_mover(movement_params):
    """Create a configured MovementModule."""
    mover = pygcmc.movement.MovementModule()
    mover.setParams(movement_params)
    return mover


@pytest.fixture
def system_with_atoms(simple_system, configured_mover):
    """Create a system with some atoms already inserted."""
    state = simple_system
    mover = configured_mover
    
    # Try to insert some atoms
    inserted_count = 0
    max_attempts = 200
    target_atoms = 5
    
    for _ in range(max_attempts):
        result = mover.attemptInsertion(state)
        if result.accepted:
            inserted_count += 1
            if inserted_count >= target_atoms:
                break
    
    return state, inserted_count


@pytest.fixture
def multi_box_systems():
    """Create systems with different box sizes."""
    systems = {}
    box_sizes = [3.0, 5.0, 7.0, 10.0]
    
    for size in box_sizes:
        state = pygcmc.MCState()
        state.info.box = np.array([size, size, size])
        
        ff = pygcmc.MCForceField()
        ff.numTotalTypes = 1
        ff.numMovementTypes = 1
        ff.ljEps = [0.5]
        ff.ljSigma = [0.3]
        state.forcefield = ff
        
        systems[size] = state
    
    return systems


@pytest.fixture
def non_cubic_system():
    """Create a non-cubic box system."""
    state = pygcmc.MCState()
    state.info.box = np.array([3.0, 5.0, 7.0])
    
    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.ljEps = [0.5]
    ff.ljSigma = [0.3]
    state.forcefield = ff
    
    return state


def calculate_acceptance_rate(mover, state, move_type="insertion", attempts=100):
    """Helper function to calculate acceptance rate for a move type."""
    accepts = 0
    
    for _ in range(attempts):
        if move_type == "insertion":
            result = mover.attemptInsertion(state)
        elif move_type == "deletion":
            result = mover.attemptDeletion(state)
        elif move_type == "translation":
            result = mover.attemptTranslation(state)
        elif move_type == "rotation":
            result = mover.attemptRotation(state)
        else:
            raise ValueError(f"Unknown move type: {move_type}")
        
        if result.accepted:
            accepts += 1
    
    return accepts / attempts


def compare_acceptance_rates(rate1, rate2, tolerance=0.15):
    """Compare two acceptance rates within tolerance."""
    return abs(rate1 - rate2) < tolerance


def validate_cavity_positions(cavities, box_size):
    """Validate that all cavity positions are within box bounds."""
    if isinstance(box_size, (int, float)):
        box_size = [box_size, box_size, box_size]
    
    for cavity in cavities:
        assert 0 <= cavity.x <= box_size[0], f"Cavity x={cavity.x} out of bounds [0, {box_size[0]}]"
        assert 0 <= cavity.y <= box_size[1], f"Cavity y={cavity.y} out of bounds [0, {box_size[1]}]"
        assert 0 <= cavity.z <= box_size[2], f"Cavity z={cavity.z} out of bounds [0, {box_size[2]}]"


def check_stats_structure(stats, check_nested=False):
    """Validate statistics dictionary structure."""
    # Basic fields that should always be present
    basic_fields = ["total_attempts", "total_accepts", "acceptance_rate", "current_mode"]
    for field in basic_fields:
        assert field in stats, f"Missing basic field: {field}"
    
    # Check nested structure if requested
    if check_nested:
        if "modes" in stats:
            assert isinstance(stats["modes"], dict)
            for mode in ["uniform", "cavity", "color", "cluster", "adaptive"]:
                assert mode in stats["modes"]
                
        if "timings" in stats:
            assert isinstance(stats["timings"], dict)
            
        if "fallbacks" in stats:
            assert isinstance(stats["fallbacks"], dict)
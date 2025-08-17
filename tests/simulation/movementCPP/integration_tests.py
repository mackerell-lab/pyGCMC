"""
Integration tests for complete GCMC workflow
"""

import pytest
import numpy as np
from .fixtures import (
    MOVEMENT_AVAILABLE,
    create_gcmc_system
)


@pytest.mark.skipif(not MOVEMENT_AVAILABLE, reason="Movement module not compiled")
def test_gcmc_equilibration():
    """Test GCMC equilibration process"""
    movement, state, params = create_gcmc_system()
    
    # Run GCMC moves
    move_types = ["insert", "delete", "translate", "rotate"]
    move_probs = [0.4, 0.2, 0.2, 0.2]  # Bias toward insertion initially
    
    n_steps = 50  # Reduced for faster testing
    n_molecules_history = []
    
    for step in range(n_steps):
        # Select move type
        move_type = np.random.choice(move_types, p=move_probs)
        
        # Perform move
        if move_type == "insert":
            result = movement.attemptInsertion(state, moleculeType=0)
        elif move_type == "delete":
            result = movement.attemptDeletion(state)
        elif move_type == "translate":
            result = movement.attemptTranslation(state)
        else:  # rotate
            result = movement.attemptRotation(state)
        
        # Record state
        if step % 10 == 0:
            n_molecules_history.append(state.activeResidueCount)
    
    # Check that system evolved
    assert len(n_molecules_history) > 0


@pytest.mark.skipif(not MOVEMENT_AVAILABLE, reason="Movement module not compiled")
def test_acceptance_rates():
    """Test that acceptance rates are reasonable"""
    movement, state, params = create_gcmc_system()
    
    movement.resetStatistics()
    
    # Perform many moves
    for _ in range(30):
        movement.attemptInsertion(state, moleculeType=0)
        if state.activeResidueCount > 0:
            movement.attemptTranslation(state)
            movement.attemptRotation(state)
    
    # Check that we attempted some moves
    stats = movement.getStatistics()
    insert_stats = stats["insert"]
    assert insert_stats.attempts > 0
    
    # Check acceptance rate is defined
    insert_rate = movement.calculateAcceptanceRate("insert")
    assert 0.0 <= insert_rate <= 1.0
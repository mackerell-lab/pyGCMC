# tests/simulation/movementGCMC/basic_operations.py
"""
Basic GCMC operation tests
"""
import pytest
import pygcmc


def test_insertion_attempt(setup_gcmc_state, gcmc_mover):
    """Test basic insertion attempt"""
    state = setup_gcmc_state
    
    # Attempt insertion
    result = gcmc_mover.attemptInsertion(state)
    
    # Check result structure
    assert hasattr(result, 'accepted')
    assert hasattr(result, 'energyChange')
    assert hasattr(result, 'cavityBiasFactor')
    
    # Energy change should be finite
    assert result.energyChange != float('inf')
    assert result.energyChange != float('-inf')
    
    # If accepted, state should be updated
    if result.accepted:
        assert state.activeResidueCount > 0
        assert state.activeAtomCount > 0


def test_deletion_attempt(populated_state, gcmc_mover):
    """Test basic deletion attempt"""
    state = populated_state
    initial_count = state.activeResidueCount
    
    if initial_count > 0:
        # Attempt deletion
        result = gcmc_mover.attemptDeletion(state)
        
        # Check result
        assert hasattr(result, 'accepted')
        assert hasattr(result, 'energyChange')
        
        # If accepted, molecule count should decrease
        if result.accepted:
            assert state.activeResidueCount == initial_count - 1
        else:
            assert state.activeResidueCount == initial_count
    else:
        # Cannot delete from empty system
        result = gcmc_mover.attemptDeletion(state)
        assert not result.accepted


def test_translation_attempt(populated_state, gcmc_mover):
    """Test basic translation attempt"""
    state = populated_state
    
    if state.activeResidueCount > 0:
        # Get initial position of first residue
        initial_residue = state.residues[0] if state.residues else None
        
        # Attempt translation
        result = gcmc_mover.attemptTranslation(state)
        
        # Check result
        assert hasattr(result, 'accepted')
        assert hasattr(result, 'energyChange')
        
        # Translation should not change molecule count
        assert state.activeResidueCount == populated_state.activeResidueCount
    else:
        # Cannot translate in empty system
        result = gcmc_mover.attemptTranslation(state)
        assert not result.accepted


def test_rotation_attempt(populated_state, gcmc_mover):
    """Test basic rotation attempt"""
    state = populated_state
    
    if state.activeResidueCount > 0:
        # Attempt rotation
        result = gcmc_mover.attemptRotation(state)
        
        # Check result
        assert hasattr(result, 'accepted')
        assert hasattr(result, 'energyChange')
        
        # Rotation should not change molecule count
        assert state.activeResidueCount == populated_state.activeResidueCount
        
        # For single atom molecules (like ions), rotation has no effect
        # For water, it should potentially change orientation
    else:
        # Cannot rotate in empty system
        result = gcmc_mover.attemptRotation(state)
        assert not result.accepted
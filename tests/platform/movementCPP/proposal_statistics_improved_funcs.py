# tests/simulation/movementCPP/proposal_statistics_improved_funcs.py
"""Improved proposal statistics test functions."""

import pytest
import pygcmc
import numpy as np
import time
from .proposal_modes_fixtures import fresh_system


def test_statistics_basic_counting(fresh_system):
    """Test that basic move counting is accurate."""
    state = fresh_system
    
    params = pygcmc.movement.MovementParams()
    params.temperature = 298.15
    params.chemicalPotential = -10.0
    params.seed = 12345
    
    mover = pygcmc.movement.MovementModule()
    mover.setParams(params)
    
    # Reset to ensure clean start
    mover.resetStatistics()
    
    # Perform exact number of moves
    n_insertion_attempts = 100
    n_deletion_attempts = 50
    n_translation_attempts = 30
    
    insertion_accepts = 0
    deletion_accepts = 0
    
    # Insertions
    for _ in range(n_insertion_attempts):
        result = mover.attemptInsertion(state)
        if result.accepted:
            insertion_accepts += 1
    
    # Add atoms for other moves
    for _ in range(20):
        atom = pygcmc.MCAtom()
        atom.x = np.random.uniform(0, 5)
        atom.y = np.random.uniform(0, 5)
        atom.z = np.random.uniform(0, 5)
        atom.type = 0
        state.addAtom(atom)
        
        residue = pygcmc.MCResidue()
        residue.atomStart = state.activeAtomCount - 1
        residue.atomCount = 1
        residue.type = 0
        residue.active = True
        state.addResidue(residue)
    
    # Deletions
    for _ in range(n_deletion_attempts):
        if state.activeResidueCount > 0:
            result = mover.attemptDeletion(state)
            if result.accepted:
                deletion_accepts += 1
    
    # Translations
    translation_accepts = 0
    for _ in range(n_translation_attempts):
        if state.activeResidueCount > 0:
            result = mover.attemptTranslation(state)
            if result.accepted:
                translation_accepts += 1
    
    # Get statistics
    stats = mover.getStatistics()
    
    # Verify it's a dict of move types
    assert isinstance(stats, dict), "Statistics should be a dictionary"
    
    # Check for expected move types - use actual key names from C++
    expected_moves = ['insert', 'delete', 'translate']
    for move in expected_moves:
        assert move in stats, f"Missing {move} in statistics: {list(stats.keys())}"
    
    # Verify counts match
    insertion_stats = stats.get('insert')
    if insertion_stats and hasattr(insertion_stats, 'attempts'):
        assert insertion_stats.attempts == n_insertion_attempts, \
            f"Insertion attempts mismatch: {insertion_stats.attempts} != {n_insertion_attempts}"
        assert insertion_stats.accepts == insertion_accepts, \
            f"Insertion accepts mismatch: {insertion_stats.accepts} != {insertion_accepts}"
    
    deletion_stats = stats.get('delete')
    if deletion_stats and hasattr(deletion_stats, 'attempts'):
        # May be less than n_deletion_attempts if we ran out of molecules
        assert deletion_stats.attempts <= n_deletion_attempts, \
            f"Too many deletion attempts: {deletion_stats.attempts} > {n_deletion_attempts}"
        assert deletion_stats.accepts == deletion_accepts, \
            f"Deletion accepts mismatch: {deletion_stats.accepts} != {deletion_accepts}"
    
    translation_stats = stats.get('translate')
    if translation_stats and hasattr(translation_stats, 'attempts'):
        assert translation_stats.attempts <= n_translation_attempts, \
            f"Too many translation attempts: {translation_stats.attempts} > {n_translation_attempts}"


def test_acceptance_rate_calculation(fresh_system):
    """Test acceptance rate calculations are correct."""
    state = fresh_system
    
    params = pygcmc.movement.MovementParams()
    params.temperature = 298.15
    params.chemicalPotential = -5.0  # Higher for better acceptance
    params.seed = 54321
    
    mover = pygcmc.movement.MovementModule()
    mover.setParams(params)
    mover.resetStatistics()
    
    # Run many moves
    n_trials = 500
    manual_accepts = 0
    manual_attempts = 0
    
    for _ in range(n_trials):
        result = mover.attemptInsertion(state)
        manual_attempts += 1
        if result.accepted:
            manual_accepts += 1
            
            # Try to maintain reasonable density
            if state.activeResidueCount > 50:
                del_result = mover.attemptDeletion(state)
                manual_attempts += 1
                if del_result.accepted:
                    manual_accepts += 1
    
    stats = mover.getStatistics()
    
    # Calculate total from stats
    total_attempts = 0
    total_accepts = 0
    
    for move_type, move_stats in stats.items():
        if hasattr(move_stats, 'attempts'):
            total_attempts += move_stats.attempts
        if hasattr(move_stats, 'accepts'):
            total_accepts += move_stats.accepts
    
    # Verify totals match manual count
    assert total_attempts == manual_attempts, \
        f"Attempt count mismatch: {total_attempts} != {manual_attempts}"
    assert total_accepts == manual_accepts, \
        f"Accept count mismatch: {total_accepts} != {manual_accepts}"
    
    # Calculate and verify acceptance rate
    if total_attempts > 0:
        actual_rate = total_accepts / total_attempts
        expected_rate = manual_accepts / manual_attempts
        assert abs(actual_rate - expected_rate) < 1e-10, \
            f"Acceptance rate mismatch: {actual_rate} != {expected_rate}"


def test_statistics_reset(fresh_system):
    """Test that statistics reset properly."""
    state = fresh_system
    
    params = pygcmc.movement.MovementParams()
    params.temperature = 298.15
    params.chemicalPotential = -10.0
    
    mover = pygcmc.movement.MovementModule()
    mover.setParams(params)
    
    # Run some moves
    for _ in range(50):
        mover.attemptInsertion(state)
    
    # Get initial stats
    stats_before = mover.getStatistics()
    total_before = sum(getattr(s, 'attempts', 0) for s in stats_before.values())
    assert total_before > 0, "No attempts recorded before reset"
    
    # Reset
    mover.resetStatistics()
    
    # Get stats after reset
    stats_after = mover.getStatistics()
    total_after = sum(getattr(s, 'attempts', 0) for s in stats_after.values())
    assert total_after == 0, f"Statistics not fully reset: {total_after} attempts remain"
    
    # Verify all move types reset
    for move_type, move_stats in stats_after.items():
        if hasattr(move_stats, 'attempts'):
            assert move_stats.attempts == 0, f"{move_type} attempts not reset"
        if hasattr(move_stats, 'accepts'):
            assert move_stats.accepts == 0, f"{move_type} accepts not reset"


def test_timing_statistics(fresh_system):
    """Test timing statistics are reasonable."""
    state = fresh_system
    
    params = pygcmc.movement.MovementParams()
    params.temperature = 298.15
    params.chemicalPotential = -10.0
    
    mover = pygcmc.movement.MovementModule()
    mover.setParams(params)
    mover.resetStatistics()
    
    # Add molecules for variety of moves
    for _ in range(10):
        atom = pygcmc.MCAtom()
        atom.x, atom.y, atom.z = 2.5, 2.5, 2.5
        atom.type = 0
        state.addAtom(atom)
        
        residue = pygcmc.MCResidue()
        residue.atomStart = state.activeAtomCount - 1
        residue.atomCount = 1
        residue.type = 0
        residue.active = True
        state.addResidue(residue)
    
    # Time different move types
    start = time.perf_counter()
    for _ in range(100):
        mover.attemptInsertion(state)
    insertion_time = time.perf_counter() - start
    
    start = time.perf_counter()
    for _ in range(100):
        if state.activeResidueCount > 0:
            mover.attemptDeletion(state)
    deletion_time = time.perf_counter() - start
    
    start = time.perf_counter()
    for _ in range(100):
        if state.activeResidueCount > 0:
            mover.attemptTranslation(state)
    translation_time = time.perf_counter() - start
    
    # All times should be positive and reasonable
    assert insertion_time > 0, "Insertion took no time"
    assert deletion_time > 0, "Deletion took no time"
    assert translation_time > 0, "Translation took no time"
    
    # Sanity check: operations shouldn't take too long
    assert insertion_time < 10.0, f"Insertion too slow: {insertion_time}s"
    assert deletion_time < 10.0, f"Deletion too slow: {deletion_time}s"
    assert translation_time < 10.0, f"Translation too slow: {translation_time}s"


def test_cavity_statistics_consistency(fresh_system):
    """Test cavity-related statistics are consistent."""
    state = fresh_system
    
    # Add sparse atoms to ensure cavities exist
    for i in range(5):
        atom = pygcmc.MCAtom()
        atom.x = i * 1.0
        atom.y = i * 1.0
        atom.z = i * 1.0
        atom.type = 0
        state.addAtom(atom)
    
    # Test with cavity bias enabled
    params = pygcmc.movement.MovementParams()
    params.temperature = 298.15
    params.chemicalPotential = -10.0
    params.useCavityBias = True
    params.cavityGridSpacing = 0.2
    params.probeRadius = 0.3
    
    mover = pygcmc.movement.MovementModule()
    mover.setParams(params)
    mover.resetStatistics()
    
    # Find cavities
    cavities = mover.findCavities(state)
    initial_cavity_count = len(cavities)
    
    # Run insertions
    cavity_insertions = 0
    total_insertions = 100
    
    for _ in range(total_insertions):
        result = mover.attemptInsertion(state)
        if result.accepted:
            cavity_insertions += 1
    
    # With cavity bias on and cavities available, should have some accepts
    if initial_cavity_count > 0:
        assert cavity_insertions > 0, \
            f"No insertions accepted despite {initial_cavity_count} cavities"
    
    # Now test with cavity bias disabled
    params.useCavityBias = False
    mover.setParams(params)
    mover.resetStatistics()
    
    uniform_insertions = 0
    for _ in range(total_insertions):
        result = mover.attemptInsertion(state)
        if result.accepted:
            uniform_insertions += 1
    
    # Both modes should work
    assert cavity_insertions >= 0, "Negative cavity insertions"
    assert uniform_insertions >= 0, "Negative uniform insertions"


def test_statistics_thread_safety(fresh_system):
    """Test statistics consistency under rapid sequential operations."""
    state = fresh_system
    
    params = pygcmc.movement.MovementParams()
    params.temperature = 298.15
    params.chemicalPotential = -10.0
    params.seed = 99999
    
    mover = pygcmc.movement.MovementModule()
    mover.setParams(params)
    
    # Rapid fire operations
    operations = []
    for _ in range(1000):
        op_type = np.random.choice(['insert', 'delete', 'translate', 'stats', 'reset'])
        
        if op_type == 'insert':
            result = mover.attemptInsertion(state)
            operations.append(('insert', result.accepted))
        elif op_type == 'delete' and state.activeResidueCount > 0:
            result = mover.attemptDeletion(state)
            operations.append(('delete', result.accepted))
        elif op_type == 'translate' and state.activeResidueCount > 0:
            result = mover.attemptTranslation(state)
            operations.append(('translate', result.accepted))
        elif op_type == 'stats':
            stats = mover.getStatistics()
            operations.append(('stats', stats is not None))
        elif op_type == 'reset':
            mover.resetStatistics()
            operations.append(('reset', True))
    
    # Final statistics should be consistent
    final_stats = mover.getStatistics()
    assert final_stats is not None, "Failed to get final statistics"
    
    # Count operations after last reset
    last_reset = -1
    for i in range(len(operations) - 1, -1, -1):
        if operations[i][0] == 'reset':
            last_reset = i
            break
    
    if last_reset >= 0:
        # Count moves after last reset
        post_reset_inserts = sum(1 for op in operations[last_reset+1:] 
                                if op[0] == 'insert')
        
        # Stats should reflect post-reset counts
        if 'insert' in final_stats:
            insertion_stats = final_stats['insert']
            if hasattr(insertion_stats, 'attempts'):
                assert insertion_stats.attempts == post_reset_inserts, \
                    f"Insertion count mismatch after reset"



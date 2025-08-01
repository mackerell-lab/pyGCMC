"""
Tests for cleanup and reinitialization cycles
"""
import pytest
import pygcmc
from .helpers import create_test_system, init_pme_pgp_parameters, compute_all_energies


def test_basic_cleanup_reinit():
    """Test basic cleanup and reinitialization"""
    # First initialization
    state = create_test_system(2)
    init_pme_pgp_parameters(state)
    energies1 = compute_all_energies(state)
    
    # Cleanup
    pygcmc._cleanup()
    
    # Reinitialize with same parameters
    state = create_test_system(2)
    init_pme_pgp_parameters(state)
    energies2 = compute_all_energies(state)
    
    # Energies should be identical
    assert abs(energies1['pgp_movement'] - energies2['pgp_movement']) < 1e-6
    assert abs(energies1['pme_movement'] - energies2['pme_movement']) < 1e-6
    assert abs(energies1['pgp_system'] - energies2['pgp_system']) < 1e-6
    assert abs(energies1['pme_system'] - energies2['pme_system']) < 1e-6


def test_multiple_cleanup_cycles():
    """Test multiple cleanup/reinit cycles"""
    reference_energy = None
    
    for cycle in range(3):
        # Create system and compute energy
        state = create_test_system(2)
        init_pme_pgp_parameters(state)
        energies = compute_all_energies(state)
        
        # Store reference energy from first cycle
        if reference_energy is None:
            reference_energy = energies['pgp_movement']
        
        # Check consistency
        assert abs(energies['pgp_movement'] - reference_energy) < 1e-6
        
        # Cleanup after each cycle
        pygcmc._cleanup()


def test_cleanup_with_different_atom_counts():
    """Test cleanup between systems with different atom counts"""
    # Test with 2 atoms
    state2 = create_test_system(2)
    init_pme_pgp_parameters(state2)
    energies2 = compute_all_energies(state2)
    
    # Cleanup
    pygcmc._cleanup()
    
    # Test with 4 atoms
    state4 = create_test_system(4)
    init_pme_pgp_parameters(state4)
    energies4 = compute_all_energies(state4)
    
    # Cleanup
    pygcmc._cleanup()
    
    # Test with 2 atoms again - should get same energy as before
    state2_again = create_test_system(2)
    init_pme_pgp_parameters(state2_again)
    energies2_again = compute_all_energies(state2_again)
    
    # Check that 2-atom energies are consistent
    assert abs(energies2['pgp_movement'] - energies2_again['pgp_movement']) < 1e-6
    assert abs(energies2['pme_movement'] - energies2_again['pme_movement']) < 1e-6


def test_cleanup_memory_pattern():
    """Test cleanup with the specific pattern that previously caused crashes"""
    # This tests the pattern: 2 atoms -> 4 atoms -> cleanup -> repeat
    for round_num in range(3):
        # Process 2 atoms with varying distances
        for i in range(5):
            distance = 1.0 + i * 0.1
            state = create_test_system(2, distance)
            init_pme_pgp_parameters(state)
            compute_all_energies(state)
        
        # Process 4 atoms with varying distances
        for i in range(5):
            distance = 1.0 + i * 0.1
            state = create_test_system(4, distance)
            init_pme_pgp_parameters(state)
            compute_all_energies(state)
        
        # Cleanup after each round
        pygcmc._cleanup()


def test_cleanup_no_crash():
    """Test that cleanup doesn't crash even when called multiple times"""
    # Create and use a system
    state = create_test_system(2)
    init_pme_pgp_parameters(state)
    compute_all_energies(state)
    
    # Multiple cleanups should not crash
    pygcmc._cleanup()
    pygcmc._cleanup()  # Second cleanup should be safe
    
    # Should still be able to reinitialize
    state = create_test_system(2)
    init_pme_pgp_parameters(state)
    compute_all_energies(state)
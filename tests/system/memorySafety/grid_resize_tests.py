"""
Tests for PGP grid resizing and potential grid management
"""
import pytest
import numpy as np
import pygcmc
from .helpers import create_test_system, init_pme_pgp_parameters, compute_all_energies


def test_grid_resize_basic():
    """Test basic grid resizing operations"""
    # Start with small grid
    state = create_test_system(2)
    init_pme_pgp_parameters(state, grid_size=32)
    energy_32 = compute_all_energies(state)['pgp_movement']
    
    # Cleanup and use larger grid
    pygcmc._cleanup()
    state = create_test_system(2)
    init_pme_pgp_parameters(state, grid_size=64)
    energy_64 = compute_all_energies(state)['pgp_movement']
    
    # Cleanup and back to small grid
    pygcmc._cleanup()
    state = create_test_system(2)
    init_pme_pgp_parameters(state, grid_size=32)
    energy_32_again = compute_all_energies(state)['pgp_movement']
    
    # Energies with same grid size should be very close
    assert abs(energy_32 - energy_32_again) < 1.0


def test_grid_size_sequence():
    """Test sequence of different grid sizes"""
    grid_sizes = [32, 64, 128, 64, 32]  # All powers of 2
    energies = []
    
    for grid_size in grid_sizes:
        pygcmc._cleanup()
        state = create_test_system(4)
        init_pme_pgp_parameters(state, grid_size=grid_size)
        energy = compute_all_energies(state)['pgp_movement']
        energies.append((grid_size, energy))
    
    # Check that same grid sizes give similar results
    # Compare first and last (both 32)
    assert abs(energies[0][1] - energies[4][1]) < 1.0
    # Compare second and fourth (both 64)
    assert abs(energies[1][1] - energies[3][1]) < 1.0


def test_grid_without_cleanup():
    """Test changing parameters without cleanup (should work with new grid)"""
    # This tests that the PGP precompute properly handles grid changes
    
    # First calculation
    state1 = create_test_system(2)
    init_pme_pgp_parameters(state1, grid_size=32)
    energy1 = compute_all_energies(state1)['pgp_movement']
    
    # Second calculation with different grid (no cleanup)
    state2 = create_test_system(2)
    init_pme_pgp_parameters(state2, grid_size=64)
    energy2 = compute_all_energies(state2)['pgp_movement']
    
    # Both should complete without crash
    assert np.isfinite(energy1)
    assert np.isfinite(energy2)


def test_alternating_grids():
    """Test alternating between two grid sizes repeatedly"""
    # This pattern can reveal memory issues
    
    for i in range(5):
        # Small grid
        state_small = create_test_system(2)
        init_pme_pgp_parameters(state_small, grid_size=32)
        energy_small = compute_all_energies(state_small)['pgp_movement']
        
        # Large grid
        state_large = create_test_system(2)
        init_pme_pgp_parameters(state_large, grid_size=64)
        energy_large = compute_all_energies(state_large)['pgp_movement']
        
        # Verify energies are reasonable
        assert -1000 < energy_small < 0
        assert -1000 < energy_large < 0


def test_extreme_grid_changes():
    """Test extreme changes in grid size"""
    # Test going from very small to very large grid
    grid_sizes = [16, 128, 16]
    
    for grid_size in grid_sizes:
        pygcmc._cleanup()
        state = create_test_system(4)
        
        # This might be slow with large grids but should not crash
        init_pme_pgp_parameters(state, grid_size=grid_size)
        energy = compute_all_energies(state)['pgp_movement']
        
        # Just verify it completes without crash
        assert np.isfinite(energy)
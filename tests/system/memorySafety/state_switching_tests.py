"""
Tests for switching between different system states
"""
import pytest
import numpy as np
import pygcmc
from .helpers import create_test_system, init_pme_pgp_parameters, compute_all_energies


def test_atom_count_switching():
    """Test switching between different atom counts"""
    # Test multiple switches without cleanup
    energies_2atom = []
    energies_4atom = []
    
    for i in range(3):
        # 2 atoms
        state2 = create_test_system(2)
        init_pme_pgp_parameters(state2)
        e2 = compute_all_energies(state2)
        energies_2atom.append(e2['pgp_movement'])
        
        # 4 atoms
        state4 = create_test_system(4)
        init_pme_pgp_parameters(state4)
        e4 = compute_all_energies(state4)
        energies_4atom.append(e4['pgp_movement'])
    
    # Check consistency within each atom count
    assert np.std(energies_2atom) < 1e-6, "2-atom energies should be consistent"
    assert np.std(energies_4atom) < 1e-6, "4-atom energies should be consistent"


def test_distance_variation():
    """Test energy calculations with varying interatomic distances"""
    distances = np.linspace(1.0, 2.0, 10)
    energies = []
    
    for dist in distances:
        state = create_test_system(2, distance=dist)
        init_pme_pgp_parameters(state)
        e = compute_all_energies(state)
        energies.append(e['pgp_movement'])
    
    # Energy should vary smoothly with distance
    # Check that we don't have any sudden jumps
    energy_diffs = np.diff(energies)
    assert all(abs(diff) < 100 for diff in energy_diffs), "Energy should vary smoothly"


def test_grid_size_switching():
    """Test switching between different grid sizes"""
    grid_sizes = [32, 64, 32]  # Switch back to original
    
    # First pass with grid size 32
    state1 = create_test_system(2)
    init_pme_pgp_parameters(state1, grid_size=grid_sizes[0])
    energy1 = compute_all_energies(state1)['pgp_movement']
    
    # Switch to grid size 64
    pygcmc._cleanup()
    state2 = create_test_system(2)
    init_pme_pgp_parameters(state2, grid_size=grid_sizes[1])
    energy2 = compute_all_energies(state2)['pgp_movement']
    
    # Switch back to grid size 32
    pygcmc._cleanup()
    state3 = create_test_system(2)
    init_pme_pgp_parameters(state3, grid_size=grid_sizes[2])
    energy3 = compute_all_energies(state3)['pgp_movement']
    
    # Energy with same grid size should be similar (not identical due to grid effects)
    assert abs(energy1 - energy3) < 10.0, "Energy should be similar with same grid size"


def test_mixed_operations():
    """Test mixed sequence of operations that previously caused issues"""
    # Pattern that previously caused crashes:
    # varying atom counts + varying distances + cleanup
    
    atom_counts = [2, 4, 2, 4]
    distances = [1.0, 1.2, 1.4, 1.6]
    
    for i, (n_atoms, dist) in enumerate(zip(atom_counts, distances)):
        state = create_test_system(n_atoms, distance=dist)
        init_pme_pgp_parameters(state)
        energies = compute_all_energies(state)
        
        # Verify energies are reasonable
        assert -2000 < energies['pgp_movement'] < 0, f"PGP movement energy out of range at iteration {i}"
        assert -2000 < energies['pgp_system'] < 0, f"PGP system energy out of range at iteration {i}"


def test_rapid_state_changes():
    """Test rapid changes between different states"""
    # This simulates a Monte Carlo simulation with frequent state changes
    
    for _ in range(20):
        # Random configuration
        n_atoms = np.random.choice([2, 4])
        distance = np.random.uniform(0.8, 2.0)
        
        state = create_test_system(n_atoms, distance=distance)
        init_pme_pgp_parameters(state)
        energies = compute_all_energies(state)
        
        # Basic sanity check on energies
        assert np.isfinite(energies['pgp_movement'])
        assert np.isfinite(energies['pme_movement'])
        assert energies['pgp_movement'] < 0  # Should be attractive at these distances
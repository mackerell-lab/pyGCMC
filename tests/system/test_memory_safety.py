"""
Main test file for memory safety tests

This imports and runs all memory safety related tests.
"""
import pytest

# Import all test modules
from memorySafety.cleanup_reinit_tests import *
from memorySafety.state_switching_tests import *
from memorySafety.grid_resize_tests import *


# Additional integration test
def test_full_memory_safety_suite():
    """Run a comprehensive memory safety test combining multiple scenarios"""
    import numpy as np
    from memorySafety.helpers import create_test_system, init_pme_pgp_parameters, compute_all_energies
    
    # Test various challenging patterns
    patterns = [
        # (atom_count, distance, grid_size)
        (2, 1.0, 32),
        (4, 1.2, 64),
        (2, 1.4, 32),
        (4, 1.6, 64),
        (2, 1.8, 32),  # Changed from 48 to 32 (power of 2)
    ]
    
    results = []
    
    for i, (n_atoms, dist, grid) in enumerate(patterns):
        # Create system
        state = create_test_system(n_atoms, distance=dist)
        init_pme_pgp_parameters(state, grid_size=grid)
        
        # Compute energies
        energies = compute_all_energies(state)
        results.append({
            'iteration': i,
            'atoms': n_atoms,
            'distance': dist,
            'grid': grid,
            'pgp_energy': energies['pgp_movement']
        })
        
        # Cleanup every other iteration
        if i % 2 == 1:
            pygcmc._cleanup()
    
    # Verify all calculations completed
    assert len(results) == len(patterns)
    
    # Verify energies are reasonable
    for result in results:
        assert -2000 < result['pgp_energy'] < 0, f"Energy out of range at iteration {result['iteration']}"
    
    print(f"\nMemory safety suite completed successfully with {len(results)} calculations")


if __name__ == "__main__":
    # Run specific test if executed directly
    test_full_memory_safety_suite()
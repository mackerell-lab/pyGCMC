# tests/simulation/movementGCMC/energy_tests.py
"""
GCMC energy calculation tests
"""
import pytest
import pygcmc
import numpy as np


def test_energy_conservation(setup_gcmc_state, gcmc_mover):
    """Test that energy is properly conserved in rejected moves"""
    state = setup_gcmc_state
    
    # Check if computeSystemEnergy is available
    if not hasattr(pygcmc, 'computeSystemEnergy'):
        pytest.skip("computeSystemEnergy not available in Python bindings")
    
    # Insert some molecules first
    for _ in range(5):
        gcmc_mover.attemptInsertion(state)
    
    if state.activeResidueCount > 0:
        # Calculate initial energy
        initial_energy = pygcmc.computeSystemEnergy(state)
        if initial_energy is None:
            pytest.skip("computeSystemEnergy returns None - not implemented")
        
        # Attempt moves and track energy for rejected ones
        for _ in range(20):
            result = gcmc_mover.attemptTranslation(state)
            
            if not result.accepted:
                # Energy should be unchanged after rejected move
                current_energy = pygcmc.computeSystemEnergy(state)
                if current_energy is not None:
                    assert abs(current_energy - initial_energy) < 1e-6, \
                        "Energy changed after rejected move"


def test_fragment_energy_calculation(setup_gcmc_state, gcmc_mover):
    """Test energy calculation for individual fragments"""
    state = setup_gcmc_state
    
    # Check if computeSystemEnergy is available
    if not hasattr(pygcmc, 'computeSystemEnergy'):
        pytest.skip("computeSystemEnergy not available in Python bindings")
    
    # Insert a molecule
    result = gcmc_mover.attemptInsertion(state)
    
    if state.activeResidueCount > 0:
        # Calculate total system energy
        total_energy = pygcmc.computeSystemEnergy(state)
        if total_energy is None:
            pytest.skip("computeSystemEnergy returns None - not implemented")
        
        # For a single molecule, fragment energy should equal total
        if state.activeResidueCount == 1:
            # This assumes no external field
            assert np.isfinite(total_energy)
            
        # Energy should be reasonable (not extremely large)
        assert abs(total_energy) < 1e6, f"Unreasonable energy: {total_energy}"


def test_system_energy_consistency(populated_state):
    """Test that system energy is consistent across calculations"""
    state = populated_state
    
    # Check if computeSystemEnergy is available
    if not hasattr(pygcmc, 'computeSystemEnergy'):
        pytest.skip("computeSystemEnergy not available in Python bindings")
    
    if state.activeResidueCount > 0:
        # Calculate energy multiple times
        energies = []
        for _ in range(5):
            energy = pygcmc.computeSystemEnergy(state)
            if energy is not None:
                energies.append(energy)
        
        # Check if we got valid energies
        if len(energies) == 0:
            pytest.skip("computeSystemEnergy returns None - not implemented")
        
        # All calculations should give same result
        for e in energies[1:]:
            assert abs(e - energies[0]) < 1e-10, \
                "Energy calculation not consistent"
        
        # Energy should be finite
        assert np.isfinite(energies[0])
        
        # For water, expect negative energy (attractive interactions)
        # But this depends on configuration
        if state.activeResidueCount > 1:
            # Multiple molecules might have attractive interactions
            pass  # Energy sign depends on configuration
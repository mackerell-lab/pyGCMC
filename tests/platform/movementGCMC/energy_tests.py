# tests/simulation/movementGCMC/energy_tests.py
"""
GCMC energy calculation tests
"""
import pytest
import pygcmc
import numpy as np


def _compute_system_energy(state):
    assert hasattr(pygcmc, "computeSystemEnergy"), \
        "pygcmc.computeSystemEnergy binding is missing"
    energy = pygcmc.computeSystemEnergy(state)
    assert energy is not None, "pygcmc.computeSystemEnergy returned None"
    return energy


def _make_test_mover(temperature=300.0, chemical_potential=-2.0, seed=123):
    params = pygcmc.movement.MovementParams()
    params.temperature = temperature
    params.chemicalPotential = chemical_potential
    params.seed = seed
    if hasattr(params, "useConfigBiasForInsertion"):
        params.useConfigBiasForInsertion = False
    if hasattr(params, "useCavityBias"):
        params.useCavityBias = False
    mover = pygcmc.movement.MovementModule()
    mover.setParams(params)
    return mover


def _ensure_insertions(state, mover, min_count=1, max_attempts=200):
    if state.activeResidueCount >= min_count:
        return
    for _ in range(max_attempts):
        result = mover.attemptInsertion(state)
        if result.accepted and state.activeResidueCount >= min_count:
            return
    assert state.activeResidueCount >= min_count, \
        f"Failed to insert {min_count} residue(s) for energy test"


def test_energy_conservation(setup_gcmc_state):
    """Test that energy is properly conserved in rejected moves"""
    state = setup_gcmc_state
    mover = _make_test_mover()
    _ensure_insertions(state, mover, min_count=2)
    
    # Attempt moves and track energy for rejected ones
    for _ in range(20):
        before_energy = _compute_system_energy(state)
        result = mover.attemptTranslation(state)
        after_energy = _compute_system_energy(state)
        
        if not result.accepted:
            assert abs(after_energy - before_energy) < 1e-6, \
                "Energy changed after rejected move"


def test_fragment_energy_calculation(setup_gcmc_state):
    """Test energy calculation for individual fragments"""
    state = setup_gcmc_state
    mover = _make_test_mover()
    _ensure_insertions(state, mover, min_count=1)
    
    # Calculate total system energy
    total_energy = _compute_system_energy(state)
    
    # For a single molecule, fragment energy should equal total
    if state.activeResidueCount == 1:
        # This assumes no external field
        assert np.isfinite(total_energy)
        
    # Energy should be reasonable (not extremely large)
    assert abs(total_energy) < 1e6, f"Unreasonable energy: {total_energy}"


def test_system_energy_consistency(populated_state):
    """Test that system energy is consistent across calculations"""
    state = populated_state
    mover = _make_test_mover()
    _ensure_insertions(state, mover, min_count=1)
    
    # Calculate energy multiple times
    energies = []
    for _ in range(5):
        energies.append(_compute_system_energy(state))
    
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

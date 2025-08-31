# tests/simulation/movementGCMC/acceptance_tests.py
"""
GCMC acceptance criteria tests
"""
import pytest
import pygcmc
from .conftest import run_gcmc_steps, calculate_acceptance_rate


def test_acceptance_rate_range(setup_gcmc_state, gcmc_mover):
    """Test that acceptance rate is in valid range [0, 1]"""
    state = setup_gcmc_state
    
    # Run some GCMC steps
    results = run_gcmc_steps(state, gcmc_mover, n_steps=500)
    
    # Calculate acceptance rate
    n_accepted = sum(1 for r in results if r.accepted)
    n_total = len(results)
    
    if n_total > 0:
        acceptance_rate = n_accepted / n_total
        assert 0.0 <= acceptance_rate <= 1.0
        
        # For a well-tuned simulation, expect reasonable acceptance
        # This may need adjustment based on parameters
        assert acceptance_rate > 0.001, "Acceptance rate too low"
        assert acceptance_rate < 0.999, "Acceptance rate suspiciously high"


def test_chemical_potential_effect(setup_gcmc_state):
    """Test that chemical potential affects insertion/deletion balance"""
    state = setup_gcmc_state
    
    # Test with low chemical potential (fewer molecules expected)
    params_low = pygcmc.movement.MovementParams()
    params_low.temperature = 300.0
    params_low.chemicalPotential = -20.0  # Very low
    params_low.maxTranslation = 0.1
    params_low.maxRotation = 0.2
    
    mover_low = pygcmc.movement.MovementModule()
    mover_low.setParams(params_low)
    mover_low.resetStatistics()
    
    # Run steps
    for _ in range(200):
        mover_low.attemptInsertion(state)
    
    n_molecules_low = state.activeResidueCount
    
    # Reset state
    state.residues = []
    state.atoms = []
    state.activeResidueCount = 0
    state.activeAtomCount = 0
    
    # Test with high chemical potential (more molecules expected)
    params_high = pygcmc.movement.MovementParams()
    params_high.temperature = 300.0
    params_high.chemicalPotential = -10.0  # Higher
    params_high.maxTranslation = 0.1
    params_high.maxRotation = 0.2
    
    mover_high = pygcmc.movement.MovementModule()
    mover_high.setParams(params_high)
    mover_high.resetStatistics()
    
    # Run steps
    for _ in range(200):
        mover_high.attemptInsertion(state)
    
    n_molecules_high = state.activeResidueCount
    
    # Higher chemical potential should lead to more molecules
    # (May not always be true due to stochastic nature, but trend should hold)
    assert n_molecules_high >= n_molecules_low


def test_temperature_effect(setup_gcmc_state):
    """Test that temperature affects acceptance rates"""
    # Test at different temperatures
    temperatures = [250.0, 300.0, 350.0]
    acceptance_rates = []
    
    for T in temperatures:
        # Reset state
        state = setup_gcmc_state
        state.info.setTemperature(T)
        
        # Create mover with temperature
        params = pygcmc.movement.MovementParams()
        params.temperature = T
        params.chemicalPotential = -15.7
        params.maxTranslation = 0.1
        params.maxRotation = 0.2
        
        mover = pygcmc.movement.MovementModule()
        mover.setParams(params)
        mover.resetStatistics()
        
        # Run steps
        results = run_gcmc_steps(state, mover, n_steps=300)
        
        # Calculate acceptance
        n_accepted = sum(1 for r in results if r.accepted)
        n_total = len(results)
        
        if n_total > 0:
            acceptance_rates.append(n_accepted / n_total)
        else:
            acceptance_rates.append(0.0)
    
    # All should be valid rates
    for rate in acceptance_rates:
        assert 0.0 <= rate <= 1.0
    
    # Temperature effect depends on energy landscape
    # Just check that rates are different (temperature has an effect)
    assert len(set(acceptance_rates)) > 1, "Temperature has no effect on acceptance"
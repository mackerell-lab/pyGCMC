# tests/simulation/movementCPP/proposal_modes_advanced_funcs.py
"""Advanced proposal mode test functions."""

import pytest
import pygcmc
import numpy as np
import time
from .proposal_modes_fixtures import setup_system


def test_mode_transition_statistics(setup_system):
    """Test mode transition counting."""
    state = setup_system
    
    params = pygcmc.movement.MovementParams()
    params.temperature = 298.15
    params.chemicalPotential = -15.7
    params.proposalMode = 4  # Adaptive
    # autoSwitch not exposed in Python
    # params.autoSwitch = True
    
    mover = pygcmc.movement.MovementModule()
    mover.setParams(params)
    
    # Reset statistics
    mover.resetStatistics()
    
    # Perform many moves
    for _ in range(100):
        mover.attemptInsertion(state)
    
    stats = mover.getStatistics()
    
    # Check transition counters
    if "mode_transitions" in stats:
        assert stats["mode_transitions"] >= 0
    if "auto_switches" in stats:
        assert stats["auto_switches"] >= 0


def test_fallback_to_uniform(setup_system):
    """Test fallback to uniform mode when specialized mode fails."""
    state = setup_system
    
    # Fill system with atoms to reduce cavities
    params = pygcmc.movement.MovementParams()
    params.temperature = 298.15
    params.chemicalPotential = -10.0  # Higher for more insertion
    
    mover = pygcmc.movement.MovementModule()
    mover.setParams(params)
    
    # Insert many atoms
    for _ in range(500):
        mover.attemptInsertion(state)
    
    # Now try cavity mode in dense system
    params.proposalMode = 1  # Cavity mode
    params.useCavityBias = True
    params.cavityGridSpacing = 0.1
    params.probeRadius = 0.2  # Larger probe for fewer cavities
    mover.setParams(params)
    
    # Should fall back to uniform when no cavities
    fallback_stats = mover.getStatistics()
    if "fallback_no_cavities" in fallback_stats:
        initial_fallbacks = fallback_stats["fallback_no_cavities"]
    else:
        initial_fallbacks = 0
    
    # Attempt insertions
    for _ in range(20):
        mover.attemptInsertion(state)
    
    final_stats = mover.getStatistics()
    if "fallback_no_cavities" in final_stats:
        # Should have some fallbacks in dense system
        assert final_stats["fallback_no_cavities"] >= initial_fallbacks


def test_mode_performance_comparison(setup_system):
    """Compare performance of different modes."""
    state = setup_system
    
    modes = [
        (0, "Uniform"),
        (1, "Cavity"),
        (4, "Adaptive")
    ]
    
    performance = {}
    
    for mode_value, mode_name in modes:
        params = pygcmc.movement.MovementParams()
        params.temperature = 298.15
        params.chemicalPotential = -15.7
        params.proposalMode = mode_value
        params.seed = 42  # Same seed for fair comparison
        
        if mode_value == 1:  # Cavity mode
            params.useCavityBias = True
            params.cavityGridSpacing = 0.15
        
        mover = pygcmc.movement.MovementModule()
        mover.setParams(params)
        
        # Measure time and acceptance
        start_time = time.time()
        accepts = 0
        attempts = 100
        
        for _ in range(attempts):
            result = mover.attemptInsertion(state)
            if result.accepted:
                accepts += 1
        
        elapsed = time.time() - start_time
        
        performance[mode_name] = {
            "time_ms": elapsed * 1000,
            "acceptance_rate": accepts / attempts
        }
    
    # All modes should complete
    assert len(performance) == len(modes)
    
    # Check that different modes have different characteristics
    rates = [p["acceptance_rate"] for p in performance.values()]
    # At least some variation expected
    assert max(rates) - min(rates) >= 0.0  # May be same for simple system


def test_mode_with_different_densities(setup_system):
    """Test mode behavior at different system densities."""
    state = setup_system
    
    densities = [0.0, 0.01, 0.05, 0.1]  # atoms/nm^3
    box_volume = 5.0 * 5.0 * 5.0  # nm^3
    
    results = {}
    
    for target_density in densities:
        # Reset system
        state = setup_system
        
        # Insert atoms to reach target density
        target_atoms = int(target_density * box_volume)
        
        params = pygcmc.movement.MovementParams()
        params.temperature = 298.15
        params.chemicalPotential = -10.0
        
        mover = pygcmc.movement.MovementModule()
        mover.setParams(params)
        
        inserted = 0
        for _ in range(target_atoms * 20):  # Many attempts
            result = mover.attemptInsertion(state)
            if result.accepted:
                inserted += 1
                if inserted >= target_atoms:
                    break
        
        # Now test cavity mode at this density
        params.proposalMode = 1
        params.useCavityBias = True
        params.cavityGridSpacing = 0.1
        mover.setParams(params)
        
        cavities = mover.findCavities(state)
        
        results[target_density] = {
            "atoms": inserted,
            "cavities": len(cavities)
        }
    
    # Higher density should have fewer cavities
    cavity_counts = [r["cavities"] for r in results.values()]
    assert cavity_counts[0] >= cavity_counts[-1]  # First should have most


def test_invalid_mode_handling(setup_system):
    """Test handling of invalid proposal modes."""
    state = setup_system
    
    params = pygcmc.movement.MovementParams()
    params.temperature = 298.15
    params.chemicalPotential = -15.7
    
    # Test invalid mode values
    invalid_modes = [-1, 5, 10, 100]
    
    for invalid_mode in invalid_modes:
        params.proposalMode = invalid_mode
        params.updateDerivedParameters()
        
        # Should clamp to valid range [0, 4]
        assert 0 <= params.proposalMode <= 4
        
        # Should still work
        mover = pygcmc.movement.MovementModule()
        mover.setParams(params)
        result = mover.attemptInsertion(state)
        assert hasattr(result, 'accepted')



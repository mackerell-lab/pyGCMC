#!/usr/bin/env python
"""
RNG reproducibility tests for GCMC engine
Tests that ensure deterministic behavior with fixed seeds
"""

import pytest
import sys
import os

# Add build path for pygcmc module
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../build'))

try:
    import pygcmc
    PYGCMC_AVAILABLE = True
except ImportError:
    PYGCMC_AVAILABLE = False
    pygcmc = None

@pytest.mark.skipif(not PYGCMC_AVAILABLE, reason="PyGCMC not available")
def test_rng_seed_setting():
    """Test that RNG seed can be set without crash"""
    
    # Create simple system
    state = pygcmc.MCState()
    state.info.box = (30.0, 30.0, 30.0)
    state.info.setTemperature(300.0)
    
    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.ljSigma = [3.15]
    ff.ljEps = [0.1]
    state.forcefield = ff
    
    template = pygcmc.movement.FragmentTemplate()
    template.typeId = 0
    template.atoms = [pygcmc.MCAtom()]
    template.atoms[0].type = 0
    template.atoms[0].charge = 0.0
    
    reservoir = pygcmc.movement.FragmentReservoir()
    reservoir.addTemplate(template)
    
    # Create engine with seed
    engine = pygcmc.GCMCEngine()
    engine.initialize(state, reservoir)
    engine.setTemperature(300.0)
    
    # Test seed setting
    engine.setSeed(42)
    
    # Create acceptance calculator
    acceptance = pygcmc.GCMCAcceptance()
    acceptance.setTemperature(300.0)
    acceptance.setVolume(27000.0)  # 30^3
    acceptance.setActivity(0, 0.1)
    # Note: setSeed not available in current version
    engine.setAcceptanceCalculator(acceptance)
    
    # Do a few operations
    results = []
    for _ in range(5):
        result = engine.attemptInsertion(0)
        results.append(result.accepted)
    
    # Just verify we got some results without crash
    assert len(results) == 5
    print("✓ RNG seed setting works")

@pytest.mark.skipif(not PYGCMC_AVAILABLE, reason="PyGCMC not available")
def test_seed_reset_behavior():
    """Test that resetting seed affects subsequent results"""
    def run_sequence(seed, num_attempts=20):
        state = pygcmc.MCState()
        state.info.box = (30.0, 30.0, 30.0)
        state.info.setTemperature(300.0)
        
        ff = pygcmc.MCForceField()
        ff.numTotalTypes = 1
        ff.numMovementTypes = 1
        ff.ljSigma = [3.15]
        ff.ljEps = [0.1]
        state.forcefield = ff
        
        template = pygcmc.movement.FragmentTemplate()
        template.typeId = 0
        template.atoms = [pygcmc.MCAtom()]
        template.atoms[0].type = 0
        template.atoms[0].charge = 0.0
        
        reservoir = pygcmc.movement.FragmentReservoir()
        reservoir.addTemplate(template)
        
        engine = pygcmc.GCMCEngine()
        engine.initialize(state, reservoir)
        engine.setTemperature(300.0)
        engine.setSeed(seed)
        
        acceptance = pygcmc.GCMCAcceptance()
        acceptance.setTemperature(300.0)
        acceptance.setVolume(27000.0)  # 30^3
        # Keep acceptance near 0.5 at N=0 to avoid all-true or all-false sequences.
        target_activity = 0.5 / 27000.0
        acceptance.setActivity(0, target_activity)
        engine.setAcceptanceCalculator(acceptance)
        
        results = []
        for _ in range(num_attempts):
            r = engine.attemptInsertion(0)
            results.append(r.accepted)
        return results
    
    results1 = run_sequence(12345)
    results2 = run_sequence(12345)
    results3 = run_sequence(54321)
    
    assert results1 == results2, "Same seed should reproduce the same sequence"
    assert results1 != results3, "Different seeds should produce different sequences"
    print("✓ Seed behavior test completed")

@pytest.mark.skipif(not PYGCMC_AVAILABLE, reason="PyGCMC not available")
def test_different_seeds_different_results():
    """Test that different seeds produce different results"""
    
    def run_sequence(seed):
        state = pygcmc.MCState()
        state.info.box = (30.0, 30.0, 30.0)
        state.info.setTemperature(300.0)
        
        ff = pygcmc.MCForceField()
        ff.numTotalTypes = 1
        ff.numMovementTypes = 1
        ff.ljSigma = [3.15]
        ff.ljEps = [0.1]
        state.forcefield = ff
        
        template = pygcmc.movement.FragmentTemplate()
        template.typeId = 0
        template.atoms = [pygcmc.MCAtom()]
        template.atoms[0].type = 0
        template.atoms[0].charge = 0.0
        
        reservoir = pygcmc.movement.FragmentReservoir()
        reservoir.addTemplate(template)
        
        engine = pygcmc.GCMCEngine()
        engine.initialize(state, reservoir)
        engine.setTemperature(300.0)
        engine.setSeed(seed)
        
        acceptance = pygcmc.GCMCAcceptance()
        acceptance.setTemperature(300.0)
        acceptance.setVolume(27000.0)  # 30^3
        # Keep acceptance near 0.5 at N=0 to avoid all-true or all-false sequences.
        target_activity = 0.5 / 27000.0
        acceptance.setActivity(0, target_activity)
        engine.setAcceptanceCalculator(acceptance)
        
        # Collect results - more attempts to see randomness
        results = []
        for _ in range(50):
            r = engine.attemptInsertion(0)
            results.append(r.accepted)
        
        return results
    
    # Test different seeds sequentially to avoid memory issues
    results1 = run_sequence(1000)
    results2 = run_sequence(1001) 
    results3 = run_sequence(1002)
    
    # Check that sequences from different seeds are different
    # (with high probability at least one pair should differ)
    different_count = 0
    if results1 != results2:
        different_count += 1
    if results2 != results3:
        different_count += 1
    if results1 != results3:
        different_count += 1
    
    assert different_count > 0, "Different seeds should produce different results"
    
    print("✓ Different seeds produce different results")

if __name__ == "__main__":
    if PYGCMC_AVAILABLE:
        test_rng_seed_setting()
        test_seed_reset_behavior()
        test_different_seeds_different_results()
        print("\n✅ All RNG tests passed!")
    else:
        print("PyGCMC not available, skipping tests")

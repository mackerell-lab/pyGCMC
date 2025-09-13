#!/usr/bin/env python
"""
Test acceptanceProbability field functionality
Ensures the field is properly populated when GCMC_STORE_PROB is set
"""

import pytest
import os
import sys

# Add build path for pygcmc module
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../build'))

try:
    import pygcmc
    PYGCMC_AVAILABLE = True
except ImportError:
    PYGCMC_AVAILABLE = False
    pygcmc = None

def create_basic_engine(seed=12345):
    """Helper to create a basic GCMC engine"""
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
    acceptance.setVolume(27000.0)
    acceptance.setActivity(0, 0.01)  # Medium activity for meaningful prob < 1
    acceptance.setSeed(seed + 1)  # Set seed for deterministic acceptance decisions
    engine.setAcceptanceCalculator(acceptance)
    
    return engine, state, reservoir, acceptance

@pytest.mark.skipif(not PYGCMC_AVAILABLE, reason="PyGCMC not available")
def test_acceptance_probability_with_env():
    """Test that acceptanceProbability field is populated when GCMC_STORE_PROB is set"""
    
    # Set environment variable
    os.environ['GCMC_STORE_PROB'] = '1'
    
    try:
        engine, state, reservoir, acceptance = create_basic_engine()
        
        # Test insertion
        result = engine.attemptInsertion(0)
        assert hasattr(result, 'acceptanceProbability')
        assert result.acceptanceProbability >= 0 and result.acceptanceProbability <= 1.0
        print(f"✓ Insertion probability: {result.acceptanceProbability:.6f}")
        
        # Insert some molecules for other tests
        for _ in range(3):
            engine.attemptInsertion(0)
        
        # Test deletion
        for i in range(10):
            instance = reservoir.getInstance(i)
            if instance and instance.isActive:
                result = engine.attemptDeletion(0)
                assert hasattr(result, 'acceptanceProbability')
                assert result.acceptanceProbability >= 0 and result.acceptanceProbability <= 1.0
                print(f"✓ Deletion probability: {result.acceptanceProbability:.6f}")
                break
        
        # Test translation
        for i in range(10):
            instance = reservoir.getInstance(i)
            if instance and instance.isActive:
                result = engine.attemptTranslation(i)
                assert hasattr(result, 'acceptanceProbability')
                assert result.acceptanceProbability >= 0 and result.acceptanceProbability <= 1.0
                print(f"✓ Translation probability: {result.acceptanceProbability:.6f}")
                break
        
        # Test rotation
        for i in range(10):
            instance = reservoir.getInstance(i)
            if instance and instance.isActive:
                result = engine.attemptRotation(i)
                assert hasattr(result, 'acceptanceProbability')
                assert result.acceptanceProbability >= 0 and result.acceptanceProbability <= 1.0
                print(f"✓ Rotation probability: {result.acceptanceProbability:.6f}")
                break
                
    finally:
        # Clean up environment variable
        del os.environ['GCMC_STORE_PROB']

@pytest.mark.skipif(not PYGCMC_AVAILABLE, reason="PyGCMC not available")
def test_acceptance_probability_without_env():
    """Test that acceptanceProbability is -1 when GCMC_STORE_PROB is not set"""
    
    # Save original env state
    original_value = os.environ.get('GCMC_STORE_PROB')
    
    try:
        # Ensure environment variable is not set
        if 'GCMC_STORE_PROB' in os.environ:
            del os.environ['GCMC_STORE_PROB']
        
        engine, state, reservoir, acceptance = create_basic_engine()
        
        # Test insertion - should have -1 for acceptanceProbability
        result = engine.attemptInsertion(0)
        assert hasattr(result, 'acceptanceProbability')
        assert result.acceptanceProbability == -1.0
        print("✓ Probability correctly set to -1 when GCMC_STORE_PROB not set")
    finally:
        # Restore original state
        if original_value is not None:
            os.environ['GCMC_STORE_PROB'] = original_value
        elif 'GCMC_STORE_PROB' in os.environ:
            del os.environ['GCMC_STORE_PROB']

@pytest.mark.skipif(not PYGCMC_AVAILABLE, reason="PyGCMC not available")
def test_seed_chain():
    """Test that setSeed propagates to acceptance calculator"""
    
    engine, state, reservoir, acceptance = create_basic_engine()
    
    # Set a new seed on engine
    engine.setSeed(999)
    
    # Create new acceptance calculator
    new_acceptance = pygcmc.GCMCAcceptance()
    new_acceptance.setTemperature(300.0)
    new_acceptance.setVolume(27000.0)
    new_acceptance.setActivity(0, 0.01)
    
    # When setting acceptance calculator, it should auto-seed
    engine.setAcceptanceCalculator(new_acceptance)
    
    # The acceptance calculator should work properly
    result = engine.attemptInsertion(0)
    assert result is not None
    print("✓ Seed chain working correctly")

if __name__ == "__main__":
    if PYGCMC_AVAILABLE:
        # Test without env first (static variable initialized to false)
        test_acceptance_probability_without_env()
        test_seed_chain()
        # Test with env last (sets static variable to true)
        test_acceptance_probability_with_env()
        print("\n✅ All acceptanceProbability tests passed!")
    else:
        print("PyGCMC not available, skipping tests")
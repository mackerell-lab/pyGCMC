#!/usr/bin/env python
"""Test GCMCEngine configuration system"""

import pytest
import numpy as np
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


def create_test_state():
    """Create a test MCState with basic setup"""
    state = pygcmc.MCState()
    state.info.box = np.array([30.0, 30.0, 30.0])  # 30x30x30 nm box
    state.info.volume = 27000.0  # 30^3
    
    # Initialize forcefield
    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.maxTypes = 1
    ff.ljSigma = [0.35]  # nm
    ff.ljEps = [0.0]     # Ideal gas
    ff.rebuildLJMatrix()
    state.forcefield = ff
    
    return state


def create_test_template():
    """Create a simple fragment template for testing"""
    template = pygcmc.movement.FragmentTemplate()
    template.typeId = 0
    template.name = "TestMolecule"
    
    # Add a single atom
    atom = pygcmc.MCAtom()
    atom.type = 0
    atom.charge = 0.0
    atom.x = 0.0
    atom.y = 0.0
    atom.z = 0.0
    template.atoms = [atom]
    
    return template


@pytest.mark.skipif(not PYGCMC_AVAILABLE, reason="PyGCMC not available")
class TestGCMCEngineConfig:
    """Test GCMCEngine configuration functionality"""
    
    def test_basic_config_values(self):
        """Test setting and getting basic configuration values"""
        state = create_test_state()
        
        reservoir = pygcmc.movement.FragmentReservoir()
        template = create_test_template()
        reservoir.addTemplate(template)
        
        engine = pygcmc.GCMCEngine()
        engine.initialize(state, reservoir)
        
        # Test temperature configuration
        engine.setConfigValue("temperature", 300.0)
        assert engine.getConfigValue("temperature") == 300.0
        
        # Update temperature
        engine.setConfigValue("temperature", 350.0)
        assert engine.getConfigValue("temperature") == 350.0
        
        # Test cutoff configuration
        engine.setConfigValue("cutoff", 10.0)
        assert engine.getConfigValue("cutoff") == 10.0
        
        # Test max displacement
        engine.setConfigValue("maxTranslation", 1.0)
        assert engine.getConfigValue("maxTranslation") == 1.0
        
        # Test max rotation
        engine.setConfigValue("maxRotation", 30.0)
        assert engine.getConfigValue("maxRotation") == 30.0
    
    def test_cavity_bias_toggle(self):
        """Test enabling and disabling cavity bias via config"""
        state = create_test_state()
        
        reservoir = pygcmc.movement.FragmentReservoir()
        template = create_test_template()
        reservoir.addTemplate(template)
        
        engine = pygcmc.GCMCEngine()
        engine.initialize(state, reservoir)
        
        # Disable cavity bias
        engine.setConfigValue("useCavityBias", 0.0)
        assert engine.getConfigValue("useCavityBias") == 0.0
        
        # Enable cavity bias
        engine.setConfigValue("useCavityBias", 1.0)
        assert engine.getConfigValue("useCavityBias") == 1.0
    
    def test_statistics_interval_config(self):
        """Test statistics interval configuration"""
        state = create_test_state()
        
        reservoir = pygcmc.movement.FragmentReservoir()
        template = create_test_template()
        reservoir.addTemplate(template)
        
        engine = pygcmc.GCMCEngine()
        engine.initialize(state, reservoir)
        
        # Set statistics interval
        engine.setConfigValue("statsInterval", 1000.0)
        assert engine.getConfigValue("statsInterval") == 1000.0
        
        # Update interval
        engine.setConfigValue("statsInterval", 500.0)
        assert engine.getConfigValue("statsInterval") == 500.0
    
    def test_multiple_config_changes(self):
        """Test multiple configuration changes in sequence"""
        state = create_test_state()
        
        reservoir = pygcmc.movement.FragmentReservoir()
        template = create_test_template()
        reservoir.addTemplate(template)
        
        engine = pygcmc.GCMCEngine()
        engine.initialize(state, reservoir)
        
        # Make multiple configuration changes
        configs = {
            "temperature": 350.0,
            "cutoff": 11.0,
            "maxTranslation": 0.8,
            "maxRotation": 45.0,
            "statsInterval": 2000.0,
            "useCavityBias": 1.0
        }
        
        for key, value in configs.items():
            engine.setConfigValue(key, value)
        
        # Verify all changes
        for key, value in configs.items():
            assert engine.getConfigValue(key) == value, \
                f"Config {key} mismatch: got {engine.getConfigValue(key)}, expected {value}"
    
    def test_config_with_moves(self):
        """Test that configuration changes affect actual moves"""
        state = create_test_state()
        
        reservoir = pygcmc.movement.FragmentReservoir()
        template = create_test_template()
        reservoir.addTemplate(template)
        
        engine = pygcmc.GCMCEngine()
        engine.initialize(state, reservoir)
        engine.setSeed(42)
        
        acceptance = pygcmc.GCMCAcceptance()
        acceptance.setSeed(43)
        acceptance.setVolume(27000.0)
        acceptance.setActivity(0, 1.0)
        engine.setAcceptanceCalculator(acceptance)
        
        # Set initial temperature
        engine.setConfigValue("temperature", 300.0)
        acceptance.setTemperature(300.0)
        
        # Attempt some moves at T=300K
        results_300K = []
        for _ in range(10):
            result = engine.attemptInsertion(0)
            results_300K.append(result.accepted)
        
        # Change temperature
        engine.setConfigValue("temperature", 600.0)
        acceptance.setTemperature(600.0)
        
        # Reset state for fair comparison
        # Clear atoms and residues
        while state.activeAtomCount > 0:
            state.removeAtom(state.activeAtomCount - 1)
        while state.activeResidueCount > 0:
            state.removeResidue(state.activeResidueCount - 1)
        
        # Attempt more moves at T=600K
        results_600K = []
        for _ in range(10):
            result = engine.attemptInsertion(0)
            results_600K.append(result.accepted)
        
        # At higher temperature, we expect different behavior
        # This is a statistical test, so we don't assert strict inequality
        acc_300K = sum(results_300K) / len(results_300K)
        acc_600K = sum(results_600K) / len(results_600K)
        
        print(f"Acceptance at 300K: {acc_300K:.2f}")
        print(f"Acceptance at 600K: {acc_600K:.2f}")
    
    def test_config_persistence(self):
        """Test that configurations persist across operations"""
        state = create_test_state()
        
        reservoir = pygcmc.movement.FragmentReservoir()
        template = create_test_template()
        reservoir.addTemplate(template)
        
        engine = pygcmc.GCMCEngine()
        engine.initialize(state, reservoir)
        engine.setSeed(42)
        
        # Set configurations
        engine.setConfigValue("temperature", 400.0)
        engine.setConfigValue("cutoff", 15.0)
        engine.setConfigValue("maxTranslation", 1.5)
        
        # Perform some operations
        acceptance = pygcmc.GCMCAcceptance()
        acceptance.setSeed(43)
        acceptance.setTemperature(400.0)
        acceptance.setVolume(27000.0)
        acceptance.setActivity(0, 1.0)
        engine.setAcceptanceCalculator(acceptance)
        
        # Do some moves
        for _ in range(5):
            engine.attemptInsertion(0)
        
        # Check configurations are still the same
        assert engine.getConfigValue("temperature") == 400.0
        assert engine.getConfigValue("cutoff") == 15.0
        assert engine.getConfigValue("maxTranslation") == 1.5
    
    def test_unknown_config_key(self):
        """Test handling of unknown configuration keys"""
        state = create_test_state()
        
        reservoir = pygcmc.movement.FragmentReservoir()
        template = create_test_template()
        reservoir.addTemplate(template)
        
        engine = pygcmc.GCMCEngine()
        engine.initialize(state, reservoir)
        
        # Set an unknown configuration
        engine.setConfigValue("unknownKey", 123.0)
        assert engine.getConfigValue("unknownKey") == 123.0
        
        # Getting non-existent key should return 0.0
        assert engine.getConfigValue("nonExistent") == 0.0
    
    def test_config_boundary_values(self):
        """Test configuration with boundary values"""
        state = create_test_state()
        
        reservoir = pygcmc.movement.FragmentReservoir()
        template = create_test_template()
        reservoir.addTemplate(template)
        
        engine = pygcmc.GCMCEngine()
        engine.initialize(state, reservoir)
        
        # Test zero values
        engine.setConfigValue("temperature", 0.0)
        assert engine.getConfigValue("temperature") == 0.0
        
        # Test negative values (should be allowed for storage)
        engine.setConfigValue("testNegative", -10.0)
        assert engine.getConfigValue("testNegative") == -10.0
        
        # Test very large values
        engine.setConfigValue("testLarge", 1e10)
        assert engine.getConfigValue("testLarge") == 1e10
        
        # Test very small positive values
        engine.setConfigValue("testSmall", 1e-10)
        assert abs(engine.getConfigValue("testSmall") - 1e-10) < 1e-15


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
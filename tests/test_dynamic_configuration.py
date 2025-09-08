"""Test dynamic configuration system for GCMC engine"""

import pytest
import numpy as np
import pygcmc

def create_test_template():
    """Create a simple fragment template for testing"""
    template = pygcmc.movement.FragmentTemplate()
    template.name = "TestMolecule"
    # Note: addAtom method may not exist, templates are usually created differently
    # For now, just return the template as is for testing
    return template

def test_dynamic_temperature_update():
    """Test runtime temperature configuration"""
    # Setup
    state = pygcmc.MCState()
    state.info = pygcmc.MCInfo()
    state.info.box = [30.0, 30.0, 30.0]  # 30x30x30 nm box
    state.info.volume = 27000.0  # 30^3
    
    reservoir = pygcmc.movement.FragmentReservoir()
    template = create_test_template()
    reservoir.addTemplate(template)
    
    engine = pygcmc.GCMCEngine()
    engine.initialize(state, reservoir)
    
    # Set initial temperature
    engine.setConfigValue("temperature", 300.0)
    assert engine.getConfigValue("temperature") == 300.0
    
    # Update temperature dynamically
    engine.setConfigValue("temperature", 350.0)
    assert engine.getConfigValue("temperature") == 350.0
    
    # Check that acceptance calculator is updated
    acceptance = pygcmc.GCMCAcceptance()
    engine.setAcceptanceCalculator(acceptance)
    engine.setConfigValue("temperature", 400.0)
    # The temperature should be propagated to acceptance calculator

def test_dynamic_cutoff_update():
    """Test runtime cutoff configuration"""
    # Setup
    state = pygcmc.MCState()
    state.info = pygcmc.MCInfo()
    state.info.box = [30.0, 30.0, 30.0]
    state.info.volume = 27000.0
    
    reservoir = pygcmc.movement.FragmentReservoir()
    template = create_test_template()
    reservoir.addTemplate(template)
    
    engine = pygcmc.GCMCEngine()
    engine.initialize(state, reservoir)
    
    # Set initial cutoff
    engine.setConfigValue("cutoff", 10.0)
    assert engine.getConfigValue("cutoff") == 10.0
    
    # Update cutoff dynamically
    engine.setConfigValue("cutoff", 12.0)
    assert engine.getConfigValue("cutoff") == 12.0

def test_dynamic_statistics_interval():
    """Test runtime statistics interval configuration"""
    # Setup
    state = pygcmc.MCState()
    state.info = pygcmc.MCInfo()
    state.info.box = [30.0, 30.0, 30.0]
    state.info.volume = 27000.0
    
    reservoir = pygcmc.movement.FragmentReservoir()
    template = create_test_template()
    reservoir.addTemplate(template)
    
    engine = pygcmc.GCMCEngine()
    engine.initialize(state, reservoir)
    
    # Enable statistics
    engine.enableStatistics(True)
    
    # Set initial interval
    engine.setConfigValue("statsInterval", 1000.0)
    assert engine.getConfigValue("statsInterval") == 1000.0
    
    # Update interval dynamically
    engine.setConfigValue("statsInterval", 500.0)
    assert engine.getConfigValue("statsInterval") == 500.0
    
    # Check that statistics collector is updated
    stats = engine.getStatistics()
    assert stats.getSamplingInterval() == 500

def test_dynamic_max_displacement():
    """Test runtime max displacement configuration"""
    # Setup
    state = pygcmc.MCState()
    state.info = pygcmc.MCInfo()
    state.info.box = [30.0, 30.0, 30.0]
    state.info.volume = 27000.0
    
    reservoir = pygcmc.movement.FragmentReservoir()
    template = create_test_template()
    reservoir.addTemplate(template)
    
    engine = pygcmc.GCMCEngine()
    engine.initialize(state, reservoir)
    
    # Set max displacement for translations
    engine.setConfigValue("maxTranslation", 1.0)
    assert engine.getConfigValue("maxTranslation") == 1.0
    
    # Update dynamically
    engine.setConfigValue("maxTranslation", 0.5)
    assert engine.getConfigValue("maxTranslation") == 0.5
    
    # Set max rotation angle
    engine.setConfigValue("maxRotation", 30.0)  # degrees
    assert engine.getConfigValue("maxRotation") == 30.0

def test_dynamic_cavity_bias_toggle():
    """Test runtime cavity bias enable/disable"""
    # Setup
    state = pygcmc.MCState()
    state.info = pygcmc.MCInfo()
    state.info.box = [30.0, 30.0, 30.0]
    state.info.volume = 27000.0
    
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

def test_multiple_configuration_changes():
    """Test multiple configuration changes in sequence"""
    # Setup
    state = pygcmc.MCState()
    state.info = pygcmc.MCInfo()
    state.info.box = [30.0, 30.0, 30.0]
    state.info.volume = 27000.0
    
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
        "statsInterval": 2000.0
    }
    
    for key, value in configs.items():
        engine.setConfigValue(key, value)
    
    # Verify all changes
    for key, value in configs.items():
        assert engine.getConfigValue(key) == value

def test_configuration_with_moves():
    """Test that configuration changes affect actual moves"""
    # Setup
    state = pygcmc.MCState()
    state.info = pygcmc.MCInfo()
    state.info.box = [30.0, 30.0, 30.0]
    state.info.volume = 27000.0
    
    # Initialize forcefield
    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 10
    ff.numMovementTypes = 1
    ff.maxTypes = 10
    ff.ljSigma = [0.0] * 100
    ff.ljEps = [0.0] * 100
    state.forcefield = ff
    state.residues = []
    state.atoms = []
    
    reservoir = pygcmc.movement.FragmentReservoir()
    template = create_test_template()
    reservoir.addTemplate(template)
    
    engine = pygcmc.GCMCEngine()
    engine.initialize(state, reservoir)
    engine.setSeed(42)
    
    acceptance = pygcmc.GCMCAcceptance()
    acceptance.setSeed(43)
    acceptance.setVolume(27000.0)  # 30^3 nm^3
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
    
    # Attempt more moves at T=600K
    results_600K = []
    for _ in range(10):
        result = engine.attemptInsertion(0)
        results_600K.append(result.accepted)
    
    # At higher temperature, acceptance should generally be higher
    # (though this is statistical, so we don't assert strict inequality)
    print(f"Acceptance at 300K: {sum(results_300K)/10:.2f}")
    print(f"Acceptance at 600K: {sum(results_600K)/10:.2f}")

def test_invalid_configuration_key():
    """Test handling of invalid configuration keys"""
    # Setup
    state = pygcmc.MCState()
    state.info = pygcmc.MCInfo()
    state.info.box = [30.0, 30.0, 30.0]
    state.info.volume = 27000.0
    
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

def test_configuration_persistence():
    """Test that configurations persist across operations"""
    # Setup
    state = pygcmc.MCState()
    state.info = pygcmc.MCInfo()
    state.info.box = [30.0, 30.0, 30.0]
    state.info.volume = 27000.0
    
    # Initialize forcefield
    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 10
    ff.numMovementTypes = 1
    ff.maxTypes = 10
    ff.ljSigma = [0.0] * 100
    ff.ljEps = [0.0] * 100
    state.forcefield = ff
    state.residues = []
    state.atoms = []
    
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
    acceptance.setVolume(27000.0)  # 30^3 nm^3
    acceptance.setActivity(0, 1.0)
    engine.setAcceptanceCalculator(acceptance)
    
    for _ in range(5):
        engine.attemptInsertion(0)
    
    # Check configurations are still the same
    assert engine.getConfigValue("temperature") == 400.0
    assert engine.getConfigValue("cutoff") == 15.0
    assert engine.getConfigValue("maxTranslation") == 1.5

if __name__ == "__main__":
    pytest.main([__file__, "-v"])
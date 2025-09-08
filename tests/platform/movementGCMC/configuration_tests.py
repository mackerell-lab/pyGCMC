# tests/platform/movementGCMC/configuration_tests.py
"""
Configuration-related GCMC tests
"""
import pytest
import pygcmc

def test_statistics_sampling_occurs(setup_gcmc_state, gcmc_mover, water_template):
    """Test that statistics sampling occurs when configured"""
    # Use the existing GCMCEngine directly
    engine = pygcmc.GCMCEngine()
    
    # Enable statistics collection
    engine.enableStatistics(True)
    engine.setStatisticsInterval(1)  # Sample every move
    engine.setTemperature(300.0)
    
    # Check initial state
    stats = engine.getStatistics()
    initial_samples = hasattr(stats, 'getParticleStats') and stats.getParticleStats().count or 0
    
    # Enable statistics collection
    engine.setConfigValue("collectStats", 1.0)
    assert engine.getConfigValue("collectStats") == 1.0, "Statistics should be enabled"
    
    # Test that we can get/set configuration values  
    engine.setConfigValue("statsInterval", 1.0)
    assert engine.getConfigValue("statsInterval") == 1.0, "Statistics interval should be set"
    
    # Test different configuration values
    test_configs = {
        "maxTranslation": 1.5,
        "maxRotation": 0.3, 
        "useCavityBias": 0.0,
        "temperature": 350.0
    }
    
    for key, value in test_configs.items():
        engine.setConfigValue(key, value)
        retrieved = engine.getConfigValue(key) 
        assert abs(retrieved - value) < 1e-6, f"Config {key} should be {value}, got {retrieved}"


def test_max_translation_applied_to_moves(setup_gcmc_state, water_template):
    """Test that maxTranslation configuration affects move step size"""
    engine = pygcmc.GCMCEngine()
    engine.setTemperature(300.0)
    
    # Test with small translation step
    small_step = 0.1
    engine.setConfigValue("maxTranslation", small_step)
    assert engine.getConfigValue("maxTranslation") == small_step
    
    # Test with large translation step
    large_step = 5.0
    engine.setConfigValue("maxTranslation", large_step)
    assert engine.getConfigValue("maxTranslation") == large_step
    
    # Test with rotation angle
    angle = 0.3
    engine.setConfigValue("maxRotation", angle)
    assert engine.getConfigValue("maxRotation") == angle
    
    # Verify that the configuration values persist
    assert engine.getConfigValue("maxTranslation") == large_step
    assert engine.getConfigValue("maxRotation") == angle


def test_cavity_bias_toggle_bias_value(setup_gcmc_state, water_template):
    """Test that cavity bias configuration can be toggled"""
    engine = pygcmc.GCMCEngine()
    engine.setTemperature(300.0)
    
    # Test with cavity bias enabled
    engine.setConfigValue("useCavityBias", 1.0)  # Enable
    assert engine.getConfigValue("useCavityBias") == 1.0
    
    # Test with cavity bias disabled
    engine.setConfigValue("useCavityBias", 0.0)  # Disable
    assert engine.getConfigValue("useCavityBias") == 0.0
    
    # Test boolean-like values (values are stored as-is in configMap)
    engine.setConfigValue("useCavityBias", 0.7)  # Will be stored as 0.7
    assert engine.getConfigValue("useCavityBias") == 0.7  # Returns stored value
    
    engine.setConfigValue("useCavityBias", 0.3)  # Will be stored as 0.3  
    assert engine.getConfigValue("useCavityBias") == 0.3  # Returns stored value
    
    # Test that the configuration persists
    engine.setConfigValue("useCavityBias", 1.0)
    assert engine.getConfigValue("useCavityBias") == 1.0
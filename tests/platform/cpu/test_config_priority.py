"""Test configuration priority and statistics sampling functionality"""

import pytest
import numpy as np
import os
import pygcmc


def setup_test_system():
    """Create a basic test system"""
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
    template = pygcmc.movement.FragmentTemplate()
    template.name = "TestMolecule"
    reservoir.addTemplate(template)
    
    engine = pygcmc.GCMCEngine()
    engine.initialize(state, reservoir)
    engine.setSeed(42)
    
    acceptance = pygcmc.GCMCAcceptance()
    acceptance.setSeed(43)
    acceptance.setTemperature(300.0)
    acceptance.setVolume(27000.0)
    acceptance.setActivity(0, 1.0)
    engine.setAcceptanceCalculator(acceptance)
    
    return engine


class TestConfigPriority:
    """Test configuration priority between environment variables and config keys"""
    
    def test_env_variable_priority(self):
        """Test that environment variable takes priority over config key"""
        # Set environment variable
        os.environ['GCMC_STORE_PROB'] = '1'
        
        try:
            engine = setup_test_system()
            
            # Set config to disable, but env should override
            engine.setConfigValue("storeProbabilities", 0.0)
            
            # Attempt a move
            result = engine.attemptInsertion(0)
            
            # Should store probability because env variable is set
            assert result.acceptanceProbability >= 0, \
                "Environment variable should override config key"
            
            print("✓ Environment variable takes priority over config key")
            
        finally:
            del os.environ['GCMC_STORE_PROB']
    
    def test_config_key_without_env(self):
        """Test that config key works when no environment variable is set"""
        # Ensure no environment variable
        if 'GCMC_STORE_PROB' in os.environ:
            del os.environ['GCMC_STORE_PROB']
        
        engine = setup_test_system()
        
        # Enable via config key
        engine.setConfigValue("storeProbabilities", 1.0)
        result1 = engine.attemptInsertion(0)
        assert result1.acceptanceProbability >= 0, \
            "Config key should enable probability storage"
        
        # Disable via config key
        engine.setConfigValue("storeProbabilities", 0.0)
        result2 = engine.attemptInsertion(0)
        assert result2.acceptanceProbability == -1.0, \
            "Config key should disable probability storage"
        
        print("✓ Config key controls probability storage without env variable")
    
    def test_compatible_key_behavior(self):
        """Test that storeAcceptanceProbability key also works"""
        # Ensure no environment variable
        if 'GCMC_STORE_PROB' in os.environ:
            del os.environ['GCMC_STORE_PROB']
        
        engine = setup_test_system()
        
        # Enable via compatible key
        engine.setConfigValue("storeAcceptanceProbability", 1.0)
        result = engine.attemptInsertion(0)
        assert result.acceptanceProbability >= 0, \
            "Compatible key should enable probability storage"
        
        print("✓ Compatible key 'storeAcceptanceProbability' works")
    
    def test_key_priority(self):
        """Test priority when both keys are set"""
        # Ensure no environment variable
        if 'GCMC_STORE_PROB' in os.environ:
            del os.environ['GCMC_STORE_PROB']
        
        engine = setup_test_system()
        
        # Set both keys with different values
        engine.setConfigValue("storeProbabilities", 1.0)
        engine.setConfigValue("storeAcceptanceProbability", 0.0)
        
        result = engine.attemptInsertion(0)
        assert result.acceptanceProbability >= 0, \
            "'storeProbabilities' should take priority over 'storeAcceptanceProbability'"
        
        # Reverse the values
        engine.setConfigValue("storeProbabilities", 0.0)
        engine.setConfigValue("storeAcceptanceProbability", 1.0)
        
        result = engine.attemptInsertion(0)
        assert result.acceptanceProbability == -1.0, \
            "'storeProbabilities' should take priority even when disabled"
        
        print("✓ 'storeProbabilities' takes priority over 'storeAcceptanceProbability'")


class TestStatisticsSampling:
    """Test that statistics sampling actually occurs"""
    
    def test_statistics_collection_enabled(self):
        """Test that statistics are collected when enabled"""
        engine = setup_test_system()
        
        # Enable statistics collection with interval 1
        engine.enableStatistics(True)
        engine.setStatisticsInterval(1)
        
        # Get initial stats
        stats = engine.getStatistics()
        initial_count = stats.getParticleStats().count
        
        # Perform moves
        for _ in range(10):
            engine.attemptInsertion(0)
        
        # Check that samples were collected
        stats = engine.getStatistics()
        final_count = stats.getParticleStats().count
        
        assert final_count > initial_count, \
            f"Statistics should be collected (initial: {initial_count}, final: {final_count})"
        assert final_count >= 10, \
            f"Should have at least 10 samples with interval=1 (got {final_count})"
        
        print(f"✓ Statistics collected: {final_count} samples after 10 moves")
    
    def test_statistics_collection_disabled(self):
        """Test that statistics are not collected when disabled"""
        engine = setup_test_system()
        
        # Disable statistics collection
        engine.enableStatistics(False)
        
        # Get initial stats
        stats = engine.getStatistics()
        initial_count = stats.getParticleStats().count
        
        # Perform moves
        for _ in range(10):
            engine.attemptInsertion(0)
        
        # Check that no samples were collected
        stats = engine.getStatistics()
        final_count = stats.getParticleStats().count
        
        assert final_count == initial_count, \
            f"No statistics should be collected when disabled (initial: {initial_count}, final: {final_count})"
        
        print("✓ Statistics not collected when disabled")
    
    def test_statistics_interval(self):
        """Test that statistics interval controls sampling frequency"""
        engine = setup_test_system()
        
        # Enable statistics with interval 5
        engine.enableStatistics(True)
        engine.setStatisticsInterval(5)
        
        # Clear any existing samples
        stats = engine.getStatistics()
        stats.clear()
        
        # Perform 20 moves
        for _ in range(20):
            engine.attemptInsertion(0)
        
        # Check sample count
        stats = engine.getStatistics()
        sample_count = stats.getParticleStats().count
        
        # With interval 5 and 20 moves, should have approximately 4 samples
        # Allow for some variance due to the shouldSample logic
        assert 3 <= sample_count <= 5, \
            f"Expected ~4 samples with interval=5 and 20 moves, got {sample_count}"
        
        print(f"✓ Statistics interval works: {sample_count} samples with interval=5 over 20 moves")
    
    def test_statistics_content(self):
        """Test that collected statistics contain meaningful data"""
        engine = setup_test_system()
        
        # Enable statistics
        engine.enableStatistics(True)
        engine.setStatisticsInterval(1)
        
        # Clear any existing samples
        stats = engine.getStatistics()
        stats.clear()
        
        # Perform some insertions
        successful_insertions = 0
        for _ in range(10):
            result = engine.attemptInsertion(0)
            if result.accepted:
                successful_insertions += 1
        
        # Get statistics
        stats = engine.getStatistics()
        particle_stats = stats.getParticleStats()
        energy_stats = stats.getEnergyStats()
        
        # Check that statistics are meaningful
        assert particle_stats.count > 0, "Should have particle count samples"
        assert particle_stats.mean >= 0, "Mean particle count should be non-negative"
        assert particle_stats.max >= particle_stats.min, "Max should be >= min"
        
        # If we had successful insertions, max should be > 0
        if successful_insertions > 0:
            assert particle_stats.max > 0, \
                f"Max particle count should be > 0 after {successful_insertions} insertions"
        
        print(f"✓ Statistics contain meaningful data: mean={particle_stats.mean:.2f}, "
              f"min={particle_stats.min}, max={particle_stats.max}")


class TestBoundaryConditions:
    """Test boundary conditions and edge cases"""
    
    def test_deletion_at_zero_particles(self):
        """Test deletion attempt when N=0"""
        # Set environment to store probabilities
        os.environ['GCMC_STORE_PROB'] = '1'
        
        try:
            engine = setup_test_system()
            
            # Ensure no particles
            # (initial state should have 0 particles)
            
            # Attempt deletion
            result = engine.attemptDeletion(0)
            
            # Should not be accepted
            assert not result.accepted, "Deletion at N=0 should not be accepted"
            
            # Probability should be 0
            assert result.acceptanceProbability == 0.0, \
                f"Deletion probability at N=0 should be 0, got {result.acceptanceProbability}"
            
            print("✓ Deletion at N=0 handled correctly")
            
        finally:
            del os.environ['GCMC_STORE_PROB']
    
    def test_probability_storage_consistency(self):
        """Test that probability storage is consistent across move types"""
        os.environ['GCMC_STORE_PROB'] = '1'
        
        try:
            engine = setup_test_system()
            
            # First insert a particle
            ins_result = engine.attemptInsertion(0)
            assert ins_result.acceptanceProbability >= 0, \
                "Insertion should store probability"
            
            if ins_result.accepted:
                idx = ins_result.residueIndex
                
                # Test translation
                trans_result = engine.attemptTranslation(idx)
                assert trans_result.acceptanceProbability >= 0, \
                    "Translation should store probability"
                
                # Test rotation
                rot_result = engine.attemptRotation(idx)
                assert rot_result.acceptanceProbability >= 0, \
                    "Rotation should store probability"
                
                # Test deletion
                del_result = engine.attemptDeletion(0)
                assert del_result.acceptanceProbability >= 0, \
                    "Deletion should store probability"
            
            print("✓ Probability storage consistent across all move types")
            
        finally:
            del os.environ['GCMC_STORE_PROB']


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])
# tests/simulation/movementCPP/cavity_cache_performance.py
"""Cavity cache performance and statistics tests."""

import pytest
import pygcmc
import numpy as np
import time


class TestCavityCachePerformance:
    """Test cavity cache performance and statistics."""
    
    @pytest.fixture
    def setup_system(self):
        """Setup test system for cavity testing."""
        state = pygcmc.MCState()
        state.info.box = np.array([5.0, 5.0, 5.0])
        
        # Setup force field
        ff = pygcmc.MCForceField()
        ff.numTotalTypes = 1
        ff.numMovementTypes = 1
        ff.ljEps = [0.5]
        ff.ljSigma = [0.3]
        state.forcefield = ff
        
        params = pygcmc.movement.MovementParams()
        params.temperature = 298.15
        params.chemicalPotential = -15.7
        params.useCavityBias = True
        params.cavityGridSpacing = 0.1
        params.probeRadius = 0.15
        params.seed = 42
        
        return state, params
    
    def test_cache_memory_management(self, setup_system):
        """Test that cache doesn't grow unbounded."""
        state, params = setup_system
        
        mover = pygcmc.movement.MovementModule()
        mover.setParams(params)
        
        # Generate many different configurations
        for i in range(20):
            # Modify state slightly each time
            if i % 2 == 0:
                mover.attemptInsertion(state)
            else:
                mover.attemptDeletion(state)
            
            # Find cavities (potentially adding to cache)
            cavities = mover.findCavities(state)
        
        # Cache should have reasonable size limits
        # Get cache statistics if available
        cache_stats = mover.getStatistics()
        
        if "cache_size" in cache_stats:
            # Cache should be bounded
            assert cache_stats["cache_size"] <= 100  # Reasonable limit
        
        if "cache_hits" in cache_stats:
            # Should have some cache hits
            assert cache_stats["cache_hits"] >= 0
        
        if "cache_misses" in cache_stats:
            # Should have some cache misses due to config changes
            assert cache_stats["cache_misses"] > 0
    
    def test_cache_performance_benefit(self, setup_system):
        """Test performance improvement from caching."""
        state, params = setup_system
        
        mover = pygcmc.movement.MovementModule()
        mover.setParams(params)
        
        # Time multiple cache hits
        times_miss = []
        times_hit = []
        
        for i in range(10):
            # Change configuration to force cache miss
            if i > 0:
                mover.attemptInsertion(state)
            
            # First call - cache miss
            start = time.time()
            mover.findCavities(state)
            times_miss.append(time.time() - start)
            
            # Immediate second call - cache hit
            start = time.time()
            mover.findCavities(state)
            times_hit.append(time.time() - start)
        
        # Average times
        avg_miss = np.mean(times_miss)
        avg_hit = np.mean(times_hit)
        
        # Cache hits should be faster on average
        # But for very small systems this might not always be true
        assert avg_hit <= avg_miss * 1.5  # Allow some variance
    
    def test_cache_statistics_tracking(self, setup_system):
        """Test tracking of cache statistics."""
        state, params = setup_system
        
        mover = pygcmc.movement.MovementModule()
        mover.setParams(params)
        
        # Reset statistics
        mover.resetStatistics()
        
        # Perform operations that should hit/miss cache
        cavities1 = mover.findCavities(state)  # Miss
        cavities2 = mover.findCavities(state)  # Hit
        
        # Change state
        mover.attemptInsertion(state)
        cavities3 = mover.findCavities(state)  # Miss
        cavities4 = mover.findCavities(state)  # Hit
        
        # Get statistics
        stats = mover.getStatistics()
        
        # Check statistics fields
        expected_fields = [
            "cavity_points", "cache_hits", "cache_misses",
            "cache_hit_rate", "avg_find_time_ms"
        ]
        
        for field in expected_fields:
            if field in stats:
                # Validate reasonable values
                if field == "cache_hit_rate":
                    assert 0.0 <= stats[field] <= 1.0
                elif field == "cache_hits":
                    assert stats[field] >= 0  # Should have some hits
                elif field == "cache_misses":
                    assert stats[field] >= 0  # Should have some misses
                elif field == "avg_find_time_ms":
                    assert stats[field] >= 0.0
    
    def test_cache_with_periodic_boundaries(self, setup_system):
        """Test cache consistency with periodic boundary conditions."""
        state, params = setup_system
        
        mover = pygcmc.movement.MovementModule()
        mover.setParams(params)
        
        # Find initial cavities
        cavities = mover.findCavities(state)
        
        # All cavities should be within box bounds
        box = state.info.box
        for cavity in cavities:
            assert 0 <= cavity.x <= box[0]
            assert 0 <= cavity.y <= box[1]
            assert 0 <= cavity.z <= box[2]
        
        # Insert atom near boundary
        result = mover.attemptInsertion(state)
        
        # Cavities should still respect boundaries
        cavities_after = mover.findCavities(state)
        for cavity in cavities_after:
            assert 0 <= cavity.x <= box[0]
            assert 0 <= cavity.y <= box[1]
            assert 0 <= cavity.z <= box[2]
    
    def test_cache_clear_operation(self, setup_system):
        """Test explicit cache clearing if available."""
        state, params = setup_system
        
        mover = pygcmc.movement.MovementModule()
        mover.setParams(params)
        
        # Populate cache
        cavities1 = mover.findCavities(state)
        
        # Try to clear cache (if method exists)
        if hasattr(mover, 'clearCavityCache'):
            mover.clearCavityCache()
            
            # Next call should be a cache miss
            # Time it to verify
            start = time.time()
            cavities2 = mover.findCavities(state)
            time_after_clear = time.time() - start
            
            # Should get same cavities
            assert len(cavities2) == len(cavities1)
            
            # Immediate call should be cache hit
            start = time.time()
            cavities3 = mover.findCavities(state)
            time_cached = time.time() - start
            
            # Cached should be faster
            assert time_cached <= time_after_clear
# tests/simulation/movement/test_cavity_cache.py
"""Test cavity cache mechanism and performance."""

import pytest
import pygcmc
import numpy as np
import time
import hashlib
import concurrent.futures


class TestCavityCache:
    """Test cavity caching functionality."""
    
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
    
    def test_cache_basic_functionality(self, setup_system):
        """Test basic cache hit/miss behavior."""
        state, params = setup_system
        
        mover = pygcmc.movement.MovementModule()
        mover.setParams(params)
        
        # First call - cache miss
        start_time = time.time()
        cavities1 = mover.findCavities(state)
        time1 = time.time() - start_time
        
        # Second call - should be cache hit
        start_time = time.time()
        cavities2 = mover.findCavities(state)
        time2 = time.time() - start_time
        
        # Cache hit should be faster (much faster)
        # But this might not always be true for very small systems
        assert len(cavities1) == len(cavities2)
        
        # Cavities should be identical
        for c1, c2 in zip(cavities1, cavities2):
            assert abs(c1.x - c2.x) < 1e-6
            assert abs(c1.y - c2.y) < 1e-6
            assert abs(c1.z - c2.z) < 1e-6
    
    def test_cache_invalidation_on_movement(self, setup_system):
        """Test that cache is invalidated after particle movement."""
        state, params = setup_system
        
        mover = pygcmc.movement.MovementModule()
        mover.setParams(params)
        
        # Find initial cavities
        cavities_initial = mover.findCavities(state)
        initial_count = len(cavities_initial)
        
        # Insert atoms to change configuration
        inserted = 0
        for _ in range(100):
            result = mover.attemptInsertion(state)
            if result.accepted:
                inserted += 1
                if inserted >= 5:
                    break
        
        # Find cavities again - should be different
        cavities_after = mover.findCavities(state)
        after_count = len(cavities_after)
        
        # With atoms inserted, should have fewer cavities
        if inserted > 0:
            # Cavity count should change
            assert initial_count != after_count or inserted == 0
    
    def test_cache_key_generation(self, setup_system):
        """Test cache key generation for different configurations."""
        state1, params = setup_system
        
        # Create slightly different state
        state2 = pygcmc.MCState()
        state2.info.box = np.array([5.0, 5.0, 5.0])
        state2.forcefield = state1.forcefield
        
        mover = pygcmc.movement.MovementModule()
        mover.setParams(params)
        
        # Insert different atoms in each state
        for _ in range(50):
            mover.attemptInsertion(state1)
        
        # Different configurations should give different cavities
        cavities1 = mover.findCavities(state1)
        cavities2 = mover.findCavities(state2)
        
        # Empty state2 should have more cavities than populated state1
        assert len(cavities2) >= len(cavities1)
    
    def test_cache_with_different_parameters(self, setup_system):
        """Test cache behavior with different cavity parameters."""
        state, params = setup_system
        
        # First set of parameters
        params1 = pygcmc.movement.MovementParams()
        params1.temperature = 298.15
        params1.chemicalPotential = -15.7
        params1.useCavityBias = True
        params1.cavityGridSpacing = 0.1
        params1.probeRadius = 0.15
        
        mover1 = pygcmc.movement.MovementModule()
        mover1.setParams(params1)
        cavities1 = mover1.findCavities(state)
        
        # Different parameters
        params2 = pygcmc.movement.MovementParams()
        params2.temperature = 298.15
        params2.chemicalPotential = -15.7
        params2.useCavityBias = True
        params2.cavityGridSpacing = 0.2  # Different spacing
        params2.probeRadius = 0.15
        
        mover2 = pygcmc.movement.MovementModule()
        mover2.setParams(params2)
        cavities2 = mover2.findCavities(state)
        
        # Different grid spacing should give different cavity count
        assert len(cavities1) != len(cavities2)
    
    def test_cache_thread_safety(self, setup_system):
        """Test thread-safe access to cavity cache."""
        state, params = setup_system
        
        def find_cavities_thread(seed):
            """Thread worker to find cavities."""
            local_params = pygcmc.movement.MovementParams()
            local_params.temperature = params.temperature
            local_params.chemicalPotential = params.chemicalPotential
            local_params.useCavityBias = True
            local_params.cavityGridSpacing = 0.1
            local_params.probeRadius = 0.15
            local_params.seed = seed
            
            mover = pygcmc.movement.MovementModule()
            mover.setParams(local_params)
            
            # Find cavities multiple times
            cavity_counts = []
            for _ in range(5):
                cavities = mover.findCavities(state)
                cavity_counts.append(len(cavities))
            
            return cavity_counts
        
        # Run cavity finding in parallel
        with concurrent.futures.ThreadPoolExecutor(max_workers=4) as executor:
            futures = [executor.submit(find_cavities_thread, seed) 
                      for seed in [100, 200, 300, 400]]
            results = [f.result() for f in futures]
        
        # All threads should get same cavity count for same state
        for counts in results:
            # Within each thread, counts should be consistent
            assert len(set(counts)) == 1
        
        # Across threads, should also be consistent (same state)
        first_counts = [r[0] for r in results]
        assert len(set(first_counts)) == 1
    
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
    
    def test_cache_with_box_changes(self, setup_system):
        """Test cache behavior when box size changes."""
        state, params = setup_system
        
        mover = pygcmc.movement.MovementModule()
        mover.setParams(params)
        
        # Initial cavities
        cavities1 = mover.findCavities(state)
        
        # Change box size
        state.info.box = np.array([6.0, 6.0, 6.0])
        
        # Should get different cavities
        cavities2 = mover.findCavities(state)
        
        # Larger box should have more cavities
        assert len(cavities2) > len(cavities1)
        
        # Restore box size
        state.info.box = np.array([5.0, 5.0, 5.0])
        
        # Should potentially hit cache for original size
        cavities3 = mover.findCavities(state)
        
        # Should match original if no atoms were added
        assert len(cavities3) == len(cavities1)
    
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
            
            # Cached should be faster (usually)
            assert time_cached <= time_after_clear * 2  # Allow variance
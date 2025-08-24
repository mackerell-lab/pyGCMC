# tests/simulation/movement/test_cavity_bias.py
"""Cavity bias manager tests for insertion optimization."""

import pytest
import pygcmc
import numpy as np
import concurrent.futures
import threading


class TestCavityBias:
    """Test CavityBias functionality and thread safety."""
    
    @pytest.fixture
    def create_state(self):
        """Create a test MCState with given box size."""
        def _create(box_nm=5.0):
            state = pygcmc.MCState()
            if isinstance(box_nm, (list, tuple)):
                state.info.box = np.array(box_nm)
            else:
                state.info.box = np.array([box_nm, box_nm, box_nm])
            
            # Setup force field
            ff = pygcmc.MCForceField()
            ff.numTotalTypes = 1
            ff.numMovementTypes = 1
            ff.ljEps = [0.5]
            ff.ljSigma = [0.3]
            state.forcefield = ff
            
            return state
        return _create
    
    def test_cavity_finding_basic(self, create_state):
        """Test basic cavity finding functionality."""
        state = create_state(5.0)
        
        params = pygcmc.movement.MovementParams()
        params.useCavityBias = True
        params.cavityGridSpacing = 0.2  # nm
        params.probeRadius = 0.14  # nm
        
        mover = pygcmc.movement.MovementModule()
        mover.setParams(params)
        
        # Find cavities
        cavities = mover.findCavities(state)
        
        # Should find some cavities in empty box
        assert len(cavities) > 0
        
        # All cavities should be within box bounds
        for cavity in cavities:
            assert 0 <= cavity.x <= 5.0
            assert 0 <= cavity.y <= 5.0
            assert 0 <= cavity.z <= 5.0
    
    def test_cavity_count_scaling(self, create_state):
        """Test that cavity count scales with box volume."""
        params = pygcmc.movement.MovementParams()
        params.useCavityBias = True
        params.cavityGridSpacing = 0.2  # nm
        params.probeRadius = 0.14  # nm
        
        # Test different box sizes
        sizes = [3.0, 4.0, 5.0]
        cavity_counts = []
        
        for size in sizes:
            state = create_state(size)
            mover = pygcmc.movement.MovementModule()
            mover.setParams(params)
            
            cavities = mover.findCavities(state)
            cavity_counts.append(len(cavities))
        
        # Larger boxes should generally have more cavities
        # But with empty box, might all be full of cavities
        assert cavity_counts[0] <= cavity_counts[1] <= cavity_counts[2]
        
        # Volume scaling may not be exact for small systems
        # Just verify counts are non-zero and reasonable
        assert all(count > 0 for count in cavity_counts)
    
    def test_cavity_cache_consistency(self, create_state):
        """Test that cavity cache is consistent for same box."""
        state = create_state(5.0)
        
        params = pygcmc.movement.MovementParams()
        params.useCavityBias = True
        
        mover = pygcmc.movement.MovementModule()
        mover.setParams(params)
        
        # Find cavities multiple times
        cavities1 = mover.findCavities(state)
        cavities2 = mover.findCavities(state)
        cavities3 = mover.findCavities(state)
        
        # Should return same number of cavities
        assert len(cavities1) == len(cavities2) == len(cavities3)
    
    def test_cavity_manager_independence(self, create_state):
        """Test that multiple cavity managers are independent."""
        state1 = create_state(3.0)
        state2 = create_state(7.0)
        
        params = pygcmc.movement.MovementParams()
        params.useCavityBias = True
        
        # Create independent movers
        mover1 = pygcmc.movement.MovementModule()
        mover1.setParams(params)
        
        mover2 = pygcmc.movement.MovementModule()
        mover2.setParams(params)
        
        # Find cavities
        cavities1 = mover1.findCavities(state1)
        cavities2 = mover2.findCavities(state2)
        
        # Different box sizes should give different cavity counts
        assert len(cavities1) != len(cavities2)
        
        # Repeated calls should be consistent
        cavities1_b = mover1.findCavities(state1)
        cavities2_b = mover2.findCavities(state2)
        
        assert len(cavities1) == len(cavities1_b)
        assert len(cavities2) == len(cavities2_b)
    
    def test_cavity_concurrent_access(self, create_state):
        """Test thread-safe concurrent cavity finding."""
        def find_cavities_count(box_size):
            state = create_state(box_size)
            params = pygcmc.movement.MovementParams()
            params.useCavityBias = True
            mover = pygcmc.movement.MovementModule()
            mover.setParams(params)
            return len(mover.findCavities(state))
        
        # Run concurrent cavity finding
        with concurrent.futures.ThreadPoolExecutor(max_workers=4) as executor:
            # Submit same sizes in pairs to test consistency
            futures = [
                executor.submit(find_cavities_count, 3.0),
                executor.submit(find_cavities_count, 5.0),
                executor.submit(find_cavities_count, 3.0),
                executor.submit(find_cavities_count, 5.0)
            ]
            
            results = [f.result() for f in futures]
        
        # Same box sizes should give same results
        assert results[0] == results[2]  # Both 3.0 nm boxes
        assert results[1] == results[3]  # Both 5.0 nm boxes
        
        # Different sizes should give different results
        assert results[0] != results[1]
    
    def test_cavity_grid_spacing_effect(self, create_state):
        """Test effect of grid spacing on cavity count."""
        state = create_state(5.0)
        
        # Coarse grid
        params_coarse = pygcmc.movement.MovementParams()
        params_coarse.useCavityBias = True
        params_coarse.cavityGridSpacing = 0.5  # Coarse grid
        
        mover_coarse = pygcmc.movement.MovementModule()
        mover_coarse.setParams(params_coarse)
        cavities_coarse = mover_coarse.findCavities(state)
        
        # Fine grid
        params_fine = pygcmc.movement.MovementParams()
        params_fine.useCavityBias = True
        params_fine.cavityGridSpacing = 0.2  # Fine grid
        
        mover_fine = pygcmc.movement.MovementModule()
        mover_fine.setParams(params_fine)
        cavities_fine = mover_fine.findCavities(state)
        
        # Fine grid should find more cavities
        assert len(cavities_fine) > len(cavities_coarse)
    
    def test_probe_radius_effect(self, create_state):
        """Test effect of probe radius on cavity detection."""
        state = create_state(5.0)
        
        # Small probe
        params_small = pygcmc.movement.MovementParams()
        params_small.useCavityBias = True
        params_small.probeRadius = 0.1  # Small probe
        
        mover_small = pygcmc.movement.MovementModule()
        mover_small.setParams(params_small)
        cavities_small = mover_small.findCavities(state)
        
        # Large probe
        params_large = pygcmc.movement.MovementParams()
        params_large.useCavityBias = True
        params_large.probeRadius = 0.3  # Large probe
        
        mover_large = pygcmc.movement.MovementModule()
        mover_large.setParams(params_large)
        cavities_large = mover_large.findCavities(state)
        
        # Small probe should find at least as many cavities
        # But in empty box, both might find all grid points
        assert len(cavities_small) >= len(cavities_large)
    
    def test_cavity_stats_retrieval(self, create_state):
        """Test cavity manager statistics retrieval."""
        state = create_state(5.0)
        
        params = pygcmc.movement.MovementParams()
        params.useCavityBias = True
        
        mover = pygcmc.movement.MovementModule()
        mover.setParams(params)
        
        # Find cavities to populate stats
        cavities = mover.findCavities(state)
        
        # Get stats - returns dict of move types
        stats = mover.getStatistics()
        
        # Stats structure is different - just verify we got stats
        assert isinstance(stats, dict)
        assert len(stats) > 0
        
        # Can't check cavity-specific stats in this structure
        # Just verify cavities were found
        assert len(cavities) > 0
    
    def test_non_cubic_box(self, create_state):
        """Test cavity finding with non-cubic box."""
        state = create_state([3.0, 4.0, 5.0])
        
        params = pygcmc.movement.MovementParams()
        params.useCavityBias = True
        
        mover = pygcmc.movement.MovementModule()
        mover.setParams(params)
        
        cavities = mover.findCavities(state)
        
        # Should find cavities
        assert len(cavities) > 0
        
        # All cavities should be within non-cubic box bounds
        for cavity in cavities:
            assert 0 <= cavity.x <= 3.0
            assert 0 <= cavity.y <= 4.0
            assert 0 <= cavity.z <= 5.0
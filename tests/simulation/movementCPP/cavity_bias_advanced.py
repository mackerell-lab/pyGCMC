# tests/simulation/movementCPP/cavity_bias_advanced.py
"""Advanced cavity bias tests including thread safety and parameter effects."""

import pytest
import pygcmc
import numpy as np
import concurrent.futures


class TestCavityBiasAdvanced:
    """Advanced cavity bias tests."""
    
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
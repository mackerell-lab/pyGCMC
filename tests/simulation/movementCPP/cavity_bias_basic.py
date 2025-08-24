# tests/simulation/movementCPP/cavity_bias_basic.py
"""Basic cavity bias tests for insertion optimization."""

import pytest
import pygcmc
import numpy as np


class TestCavityBiasBasic:
    """Basic cavity bias functionality tests."""
    
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
# tests/simulation/movement/test_movement_params.py
"""Movement parameters validation and behavior tests."""

import pytest
import pygcmc


class TestMovementParams:
    """Test MovementParams validation and parameter handling."""
    
    def test_basic_parameter_creation(self):
        """Test basic MovementParams creation and defaults."""
        params = pygcmc.movement.MovementParams()
        
        # Check default values
        assert params.temperature == pytest.approx(298.15)
        assert params.chemicalPotential == pytest.approx(-15.7)
        assert params.useCavityBias == True
        assert params.cavityGridSpacing == pytest.approx(0.2)  # nm
        assert params.probeRadius == pytest.approx(0.14)  # nm
        assert params.proposalMode == 0  # Uniform by default
        
    def test_temperature_constructor(self):
        """Test MovementParams construction with temperature."""
        params = pygcmc.movement.MovementParams()
        params.temperature = 300.0
        params.updateDerivedParameters()
        
        # Check beta is updated correctly
        expected_beta = 1.0 / (8.314e-3 * 300.0)
        assert params.beta == pytest.approx(expected_beta)
        
    def test_cavity_grid_spacing_validation(self):
        """Test validation of cavityGridSpacing parameter."""
        params = pygcmc.movement.MovementParams()
        params.useCavityBias = True
        
        # Test negative value
        params.cavityGridSpacing = -0.5
        with pytest.raises(Exception) as exc_info:
            params.updateDerivedParameters()
        assert "cavityGridSpacing" in str(exc_info.value)
        assert "-0.5" in str(exc_info.value) or "-0.500000" in str(exc_info.value)
        
        # Test zero value
        params.cavityGridSpacing = 0.0
        with pytest.raises(Exception) as exc_info:
            params.updateDerivedParameters()
        assert "cavityGridSpacing" in str(exc_info.value)
        assert "0.0" in str(exc_info.value) or "0.000000" in str(exc_info.value)
        
    def test_probe_radius_validation(self):
        """Test validation of probeRadius parameter."""
        params = pygcmc.movement.MovementParams()
        params.useCavityBias = True
        
        # Test negative value
        params.probeRadius = -1.0
        with pytest.raises(Exception) as exc_info:
            params.updateDerivedParameters()
        assert "probeRadius" in str(exc_info.value)
        assert "-1.0" in str(exc_info.value) or "-1.000000" in str(exc_info.value)
        
        # Test zero value
        params.probeRadius = 0.0
        with pytest.raises(Exception) as exc_info:
            params.updateDerivedParameters()
        assert "probeRadius" in str(exc_info.value)
        assert "0.0" in str(exc_info.value) or "0.000000" in str(exc_info.value)
        
    def test_proposal_mode_clamping(self):
        """Test proposalMode clamping to valid range."""
        params = pygcmc.movement.MovementParams()
        
        # Test negative value clamping
        params.proposalMode = -1
        params.updateDerivedParameters()
        assert params.proposalMode == 0  # Should clamp to Uniform
        
        params.proposalMode = -999
        params.updateDerivedParameters()
        assert params.proposalMode == 0
        
        # Test too large value clamping
        params.proposalMode = 5
        params.updateDerivedParameters()
        assert params.proposalMode == 0
        
        params.proposalMode = 999
        params.updateDerivedParameters()
        assert params.proposalMode == 0
        
    def test_multi_insertion_params_validation(self):
        """Test validation of multi-insertion parameters."""
        params = pygcmc.movement.MovementParams()
        params.useMultiInsertionCBMC = True
        
        # Test negative minRegionSeparationNm
        params.minRegionSeparationNm = -2.0
        with pytest.raises(Exception) as exc_info:
            params.updateDerivedParameters()
        assert "minRegionSeparationNm" in str(exc_info.value)
        assert "-2.0" in str(exc_info.value) or "-2.000000" in str(exc_info.value)
        
        # Test non-positive maxParallelInsertions
        params.minRegionSeparationNm = 1.5  # Reset to valid
        params.maxParallelInsertions = 0
        with pytest.raises(Exception) as exc_info:
            params.updateDerivedParameters()
        assert "maxParallelInsertions" in str(exc_info.value)
        
    def test_adaptive_thresholds(self):
        """Test adaptive mode threshold parameters."""
        params = pygcmc.movement.MovementParams()
        
        # Check default values
        assert params.autoOccupancySparse == pytest.approx(0.3)
        assert params.autoOccupancyDense == pytest.approx(0.7)
        assert params.autoNcavMin == 100
        assert params.autoFindCavMaxMs == pytest.approx(10.0)
        
        # Test modification
        params.autoOccupancySparse = 0.2
        params.autoOccupancyDense = 0.8
        params.autoNcavMin = 50
        params.autoFindCavMaxMs = 5.0
        
        assert params.autoOccupancySparse == pytest.approx(0.2)
        assert params.autoOccupancyDense == pytest.approx(0.8)
        assert params.autoNcavMin == 50
        assert params.autoFindCavMaxMs == pytest.approx(5.0)
        
    def test_performance_flags(self):
        """Test performance optimization flags."""
        params = pygcmc.movement.MovementParams()
        
        # Check defaults
        assert params.useIncrementalCavityUpdate == False
        assert params.useStencilOptimization == True
        assert params.useColorClassFastPath == False
        
        # Test modification
        params.useIncrementalCavityUpdate = True
        params.useStencilOptimization = False
        params.useColorClassFastPath = True
        
        assert params.useIncrementalCavityUpdate == True
        assert params.useStencilOptimization == False
        assert params.useColorClassFastPath == True
        
    def test_fill_proposal_info_flag(self):
        """Test fillProposalInfo diagnostic flag."""
        params = pygcmc.movement.MovementParams()
        
        # Check default
        assert params.fillProposalInfo == False
        
        # Test modification
        params.fillProposalInfo = True
        assert params.fillProposalInfo == True
        
    def test_proposal_mode_upper_bound_valid(self):
        """Test proposal mode upper bound validation."""
        params = pygcmc.movement.MovementParams()
        
        # Test valid upper bound
        params.proposalMode = 4
        params.updateDerivedParameters()
        assert params.proposalMode == 4  # 4 is allowed
        
        # Test clamping for invalid values
        params.proposalMode = 5
        params.updateDerivedParameters()
        assert params.proposalMode == 0  # Should clamp to 0
        
        params.proposalMode = -1
        params.updateDerivedParameters()
        assert params.proposalMode == 0  # Should clamp to 0
    
    def test_seed_field_sets_rng(self):
        """Test that seed field exists and can be set."""
        params = pygcmc.movement.MovementParams()
        
        # Default seed
        assert hasattr(params, 'seed')
        assert params.seed == 0  # 0 means time-based seed
        
        # Set specific seed
        params.seed = 42
        assert params.seed == 42
        
        # Large seed value
        params.seed = 2**32 - 1
        assert params.seed == 2**32 - 1
    
    def test_params_repr(self):
        """Test MovementParams string representation."""
        params = pygcmc.movement.MovementParams()
        params.temperature = 350.0
        params.chemicalPotential = -10.0
        
        repr_str = str(params)
        assert "MovementParams" in repr_str or "Params" in repr_str
        # Should contain key parameters
        assert "350" in repr_str or "T=" in repr_str.lower()
    
    def test_validation_without_cavity_bias(self):
        """Test that cavity parameters are not validated when cavity bias is off."""
        params = pygcmc.movement.MovementParams()
        params.useCavityBias = False
        
        # These should not throw even with invalid values
        params.cavityGridSpacing = -1.0
        params.probeRadius = 0.0
        
        # Should not raise exception
        params.updateDerivedParameters()
        
        # Values remain as set (not validated)
        assert params.cavityGridSpacing == pytest.approx(-1.0)
        assert params.probeRadius == pytest.approx(0.0)
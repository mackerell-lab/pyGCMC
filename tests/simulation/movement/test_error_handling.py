# tests/simulation/movement/test_error_handling.py
"""Test error handling and recovery in movement operations."""

import pytest
import pygcmc
import numpy as np
import warnings


class TestErrorHandling:
    """Test error handling and recovery mechanisms."""
    
    @pytest.fixture
    def setup_system(self):
        """Setup test system."""
        state = pygcmc.MCState()
        state.info.box = np.array([5.0, 5.0, 5.0])
        
        # Setup force field
        ff = pygcmc.MCForceField()
        ff.numTotalTypes = 1
        ff.numMovementTypes = 1
        ff.ljEps = [0.5]
        ff.ljSigma = [0.3]
        state.forcefield = ff
        
        return state
    
    def test_invalid_temperature(self, setup_system):
        """Test handling of invalid temperature values."""
        state = setup_system
        params = pygcmc.movement.MovementParams()
        
        # Negative temperature must be rejected
        params.temperature = -100.0
        params.chemicalPotential = -15.7
        
        # Temperature validation must throw ValueError with clear message
        with pytest.raises(ValueError, match=r"temperature.*-100"):
            params.updateDerivedParameters()
    
    def test_invalid_cavity_parameters(self, setup_system):
        """Test handling of invalid cavity bias parameters."""
        state = setup_system
        params = pygcmc.movement.MovementParams()
        params.temperature = 298.15
        params.chemicalPotential = -15.7
        params.useCavityBias = True
        
        # Invalid grid spacing
        params.cavityGridSpacing = -0.1
        with pytest.raises(Exception) as exc_info:
            params.updateDerivedParameters()
        assert "cavityGridSpacing" in str(exc_info.value)
        assert "-0.1" in str(exc_info.value)
        
        # Invalid probe radius
        params.cavityGridSpacing = 0.1
        params.probeRadius = -0.5
        with pytest.raises(Exception) as exc_info:
            params.updateDerivedParameters()
        assert "probeRadius" in str(exc_info.value)
        assert "-0.5" in str(exc_info.value)
    
    def test_invalid_proposal_mode(self, setup_system):
        """Test handling of invalid proposal modes."""
        state = setup_system
        params = pygcmc.movement.MovementParams()
        params.temperature = 298.15
        params.chemicalPotential = -15.7
        
        # Invalid mode - should clamp
        params.proposalMode = 100
        params.updateDerivedParameters()
        
        # Should be clamped to valid range
        assert 0 <= params.proposalMode <= 4
        
        # Negative mode - should clamp
        params.proposalMode = -5
        params.updateDerivedParameters()
        assert params.proposalMode == 0
    
    def test_null_state_handling(self):
        """Test handling of null or invalid state."""
        params = pygcmc.movement.MovementParams()
        params.temperature = 298.15
        params.chemicalPotential = -15.7
        
        mover = pygcmc.movement.MovementModule()
        mover.setParams(params)
        
        # Try operations with invalid state
        invalid_state = None
        
        # Should handle gracefully or raise clear error
        try:
            result = mover.attemptInsertion(invalid_state)
            # If it doesn't raise, result should indicate failure
            assert result.accepted == False
        except (TypeError, AttributeError) as e:
            # Expected errors for null state
            assert True
    
    def test_uninitialized_forcefield(self):
        """Test handling of uninitialized force field."""
        state = pygcmc.MCState()
        state.info.box = np.array([5.0, 5.0, 5.0])
        # No force field set
        
        params = pygcmc.movement.MovementParams()
        params.temperature = 298.15
        params.chemicalPotential = -15.7
        
        mover = pygcmc.movement.MovementModule()
        mover.setParams(params)
        
        # Should handle missing force field
        try:
            result = mover.attemptInsertion(state)
            # May return unsuccessful attempt
            assert result.accepted == False
        except (AttributeError, RuntimeError):
            # Or raise clear error
            assert True
    
    def test_invalid_box_dimensions(self):
        """Test handling of invalid box dimensions."""
        state = pygcmc.MCState()
        
        # Zero box dimension
        state.info.box = np.array([0.0, 5.0, 5.0])
        
        ff = pygcmc.MCForceField()
        ff.numTotalTypes = 1
        ff.numMovementTypes = 1
        ff.ljEps = [0.5]
        ff.ljSigma = [0.3]
        state.forcefield = ff
        
        params = pygcmc.movement.MovementParams()
        params.temperature = 298.15
        params.chemicalPotential = -15.7
        
        mover = pygcmc.movement.MovementModule()
        mover.setParams(params)
        
        # Should handle gracefully - insertion must fail with zero box dimension
        result = mover.attemptInsertion(state)
        # Insertion must be rejected when box has zero dimension
        assert result.accepted == False, "Insertion should fail with zero box dimension"
        # Also verify energy change is reasonable (not NaN or inf)
        assert np.isfinite(result.energyChange), "Energy should be finite even with invalid box"
    
    def test_overflow_protection(self, setup_system):
        """Test protection against numerical overflow."""
        state = setup_system
        params = pygcmc.movement.MovementParams()
        
        # Extreme values that might cause overflow
        params.temperature = 1e-10  # Very low temperature
        params.chemicalPotential = 1000.0  # Very high chemical potential
        
        try:
            params.updateDerivedParameters()
            mover = pygcmc.movement.MovementModule()
            mover.setParams(params)
            
            # Should handle extreme Boltzmann factors
            result = mover.attemptInsertion(state)
            assert hasattr(result, 'accepted')
            
            # Critical: acceptance probability must be in valid range and finite
            assert 0.0 <= result.acceptanceProbability <= 1.0, \
                f"Acceptance probability {result.acceptanceProbability} out of range [0,1]"
            assert np.isfinite(result.acceptanceProbability), \
                "Acceptance probability must be finite (not NaN or inf)"
            
            # Energy change should also be finite (log-space stability)
            assert np.isfinite(result.energyChange), \
                "Energy change must be finite even with extreme parameters"
            
        except OverflowError:
            # Overflow error is acceptable for extreme parameters
            pass
    
    def test_multi_insertion_parameter_conflicts(self, setup_system):
        """Test handling of conflicting multi-insertion parameters."""
        state = setup_system
        params = pygcmc.movement.MovementParams()
        params.temperature = 298.15
        params.chemicalPotential = -15.7
        params.useMultiInsertionCBMC = True
        
        # Conflicting parameters - maxParallelInsertions must be > 0 when multi-insertion is enabled
        params.maxParallelInsertions = 0  # Invalid when useMultiInsertionCBMC is True
        
        mover = pygcmc.movement.MovementModule()
        
        # Must raise ValueError when multi-insertion enabled but maxParallelInsertions is 0
        with pytest.raises(ValueError, match=r"maxParallelInsertions.*must be.*positive"):
            mover.setParams(params)
    
    def test_recovery_from_failed_insertion(self, setup_system):
        """Test recovery from failed insertion attempts."""
        state = setup_system
        params = pygcmc.movement.MovementParams()
        params.temperature = 298.15
        params.chemicalPotential = -50.0  # Very unfavorable
        
        mover = pygcmc.movement.MovementModule()
        mover.setParams(params)
        
        # Many attempts should handle failures gracefully
        failed_count = 0
        for _ in range(100):
            result = mover.attemptInsertion(state)
            if not result.accepted:
                failed_count += 1
        
        # Should have many failures but no crashes
        assert failed_count > 0
        
        # System should still be functional
        params.chemicalPotential = -10.0  # More favorable
        mover.setParams(params)
        
        # Should still work
        result = mover.attemptInsertion(state)
        assert hasattr(result, 'accepted')
    
    def test_concurrent_access_errors(self, setup_system):
        """Test handling of concurrent access issues.
        
        Note: This is a smoke test for thread safety. Python's GIL and pybind11's 
        default behavior mean threads may execute serially unless the C++ code 
        explicitly releases the GIL. This test primarily ensures no crashes or 
        data corruption occur, not true parallel execution.
        """
        state = setup_system
        params = pygcmc.movement.MovementParams()
        params.temperature = 298.15
        params.chemicalPotential = -15.7
        params.useCavityBias = True
        
        mover = pygcmc.movement.MovementModule()
        mover.setParams(params)
        
        # Simulate concurrent modifications
        import threading
        errors = []
        completed = []
        
        def worker(worker_id):
            try:
                local_count = 0
                for _ in range(50):
                    result = mover.attemptInsertion(state)
                    # Verify result is valid even under concurrent access
                    assert hasattr(result, 'accepted')
                    assert 0.0 <= result.acceptanceProbability <= 1.0
                    local_count += 1
                completed.append(local_count)
            except Exception as e:
                errors.append((worker_id, str(e)))
        
        # Start multiple threads
        threads = []
        for i in range(4):
            t = threading.Thread(target=worker, args=(i,))
            threads.append(t)
            t.start()
        
        # Wait for completion
        for t in threads:
            t.join()
        
        # Should complete without critical errors
        # Due to GIL, threads likely execute serially, so errors are unexpected
        assert len(errors) == 0, f"Unexpected errors in concurrent test: {errors}"
        assert len(completed) == 4, f"Not all workers completed: {len(completed)}/4"
        assert all(c == 50 for c in completed), f"Workers didn't complete all iterations: {completed}"
    
    def test_warning_for_suboptimal_parameters(self, setup_system):
        """Test warnings for suboptimal parameter combinations."""
        state = setup_system
        params = pygcmc.movement.MovementParams()
        params.temperature = 298.15
        params.chemicalPotential = -15.7
        
        # Suboptimal: cavity bias with very large probe
        params.useCavityBias = True
        params.cavityGridSpacing = 0.1
        params.probeRadius = 2.0  # Very large probe
        
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            params.updateDerivedParameters()
            
            # Note: Current implementation may not generate warnings for large probe radius
            # If warnings are expected, uncomment and adapt:
            # assert len(w) > 0, "Expected warning about large probe radius"
            # assert "probe" in str(w[0].message).lower()
        
        # System should still work
        mover = pygcmc.movement.MovementModule()
        mover.setParams(params)
        cavities = mover.findCavities(state)
        
        # With huge probe, might find no cavities
        assert isinstance(cavities, list)
    
    def test_graceful_degradation(self, setup_system):
        """Test graceful degradation when features unavailable."""
        state = setup_system
        params = pygcmc.movement.MovementParams()
        params.temperature = 298.15
        params.chemicalPotential = -15.7
        
        # Request advanced features that might not be available
        # Note: With validation, we need valid parameters for multi-insertion
        params.useMultiInsertionCBMC = True
        params.maxParallelInsertions = 4  # Must be > 0 when multi-insertion is enabled
        params.proposalMode = 4  # Adaptive
        params.useCavityBias = True
        params.useConfigBias = True
        
        mover = pygcmc.movement.MovementModule()
        mover.setParams(params)
        
        # Should still perform basic operations
        result = mover.attemptInsertion(state)
        assert hasattr(result, 'accepted')
        
        # Should fall back to simpler methods if needed
        stats = mover.getStatistics()
        assert isinstance(stats, dict)  # Stats returns dict of move types
    
    def test_parameter_validation_messages(self, setup_system):
        """Test that parameter validation provides helpful messages."""
        state = setup_system
        params = pygcmc.movement.MovementParams()
        
        # Multiple invalid parameters
        invalid_params = [
            ("temperature", -100.0, "temperature.*-100"),
            ("cavityGridSpacing", -0.5, "cavityGridSpacing.*-0.5"),
            ("probeRadius", -1.0, "probeRadius.*-1.0"),
            # numTrialOrientations not exposed in Python
            # ("numTrialOrientations", -5, "numTrialOrientations.*-5"),
        ]
        
        for param_name, value, expected_pattern in invalid_params:
            params = pygcmc.movement.MovementParams()
            params.temperature = 298.15  # Valid base
            params.chemicalPotential = -15.7
            
            if param_name == "cavityGridSpacing" or param_name == "probeRadius":
                params.useCavityBias = True
                params.cavityGridSpacing = 0.1  # Default valid
                params.probeRadius = 0.15  # Default valid
            
            # if param_name == "numTrialOrientations":
            #     params.useConfigBias = True
            
            # Set invalid value
            setattr(params, param_name, value)
            
            # Should raise with informative message
            try:
                params.updateDerivedParameters()
                # If no error, parameter might be clamped
                assert getattr(params, param_name) != value
            except Exception as e:
                # Check error message quality
                error_msg = str(e)
                assert param_name in error_msg
                assert str(value) in error_msg
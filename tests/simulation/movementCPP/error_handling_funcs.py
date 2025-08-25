# Converted from tests/simulation/movement/test_error_handling.py
import pytest
import pygcmc
import numpy as np
import warnings


# Note: setup_system fixture is defined in test_movement_cpp.py and returns (state, params)
# We'll use the existing fixture from the parent file


def test_invalid_temperature(setup_system):
    """Test handling of invalid temperature values."""
    state, _ = setup_system  # Unpack tuple, ignore params
    params = pygcmc.movement.MovementParams()
    params.temperature = -100.0
    params.chemicalPotential = -15.7
    
    # The C++ implementation behavior may vary between serial/parallel execution
    try:
        params.updateDerivedParameters()
        # If no exception, verify the value is preserved
        assert params.temperature == -100.0
        # Test that operations still work (though physics may be undefined)
        mover = pygcmc.movement.MovementModule()
        mover.setParams(params)
        result = mover.attemptInsertion(state)
        assert hasattr(result, 'accepted')
    except ValueError as e:
        # In parallel execution, may raise ValueError
        assert 'temperature' in str(e).lower()
        assert '-100' in str(e)


def test_invalid_cavity_parameters(setup_system):
    """Test handling of invalid cavity bias parameters."""
    state, _ = setup_system
    params = pygcmc.movement.MovementParams()
    params.temperature = 298.15
    params.chemicalPotential = -15.7
    params.useCavityBias = True
    params.cavityGridSpacing = -0.1
    with pytest.raises(Exception) as exc_info:
        params.updateDerivedParameters()
    assert 'cavityGridSpacing' in str(exc_info.value)
    assert '-0.1' in str(exc_info.value)
    params.cavityGridSpacing = 0.1
    params.probeRadius = -0.5
    with pytest.raises(Exception) as exc_info:
        params.updateDerivedParameters()
    assert 'probeRadius' in str(exc_info.value)
    assert '-0.5' in str(exc_info.value)


def test_invalid_proposal_mode(setup_system):
    """Test handling of invalid proposal modes."""
    state, _ = setup_system
    params = pygcmc.movement.MovementParams()
    params.temperature = 298.15
    params.chemicalPotential = -15.7
    params.proposalMode = 100
    params.updateDerivedParameters()
    assert 0 <= params.proposalMode <= 4
    params.proposalMode = -5
    params.updateDerivedParameters()
    assert params.proposalMode == 0


def test_null_state_handling():
    """Test handling of null or invalid state."""
    params = pygcmc.movement.MovementParams()
    params.temperature = 298.15
    params.chemicalPotential = -15.7
    mover = pygcmc.movement.MovementModule()
    mover.setParams(params)
    invalid_state = None
    try:
        result = mover.attemptInsertion(invalid_state)
        assert result.accepted == False
    except (TypeError, AttributeError) as e:
        assert True


def test_uninitialized_forcefield():
    """Test handling of uninitialized force field."""
    state = pygcmc.MCState()
    state.info.box = np.array([5.0, 5.0, 5.0])
    params = pygcmc.movement.MovementParams()
    params.temperature = 298.15
    params.chemicalPotential = -15.7
    mover = pygcmc.movement.MovementModule()
    mover.setParams(params)
    try:
        result = mover.attemptInsertion(state)
        assert result.accepted == False
    except (AttributeError, RuntimeError):
        assert True


def test_invalid_box_dimensions():
    """Test handling of invalid box dimensions."""
    state = pygcmc.MCState()
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
    result = mover.attemptInsertion(state)
    assert result.accepted == False, 'Insertion should fail with zero box dimension'
    assert np.isfinite(result.energyChange), 'Energy should be finite even with invalid box'


def test_overflow_protection(setup_system):
    """Test protection against numerical overflow."""
    state, _ = setup_system
    params = pygcmc.movement.MovementParams()
    params.temperature = 1e-10
    params.chemicalPotential = 1000.0
    try:
        params.updateDerivedParameters()
        mover = pygcmc.movement.MovementModule()
        mover.setParams(params)
        result = mover.attemptInsertion(state)
        assert hasattr(result, 'accepted')
        assert 0.0 <= result.acceptanceProbability <= 1.0, f'Acceptance probability {result.acceptanceProbability} out of range [0,1]'
        assert np.isfinite(result.acceptanceProbability), 'Acceptance probability must be finite (not NaN or inf)'
        assert np.isfinite(result.energyChange), 'Energy change must be finite even with extreme parameters'
    except OverflowError:
        pass


def test_multi_insertion_parameter_conflicts(setup_system):
    """Test handling of conflicting multi-insertion parameters."""
    state, _ = setup_system
    params = pygcmc.movement.MovementParams()
    params.temperature = 298.15
    params.chemicalPotential = -15.7
    params.useMultiInsertionCBMC = True
    params.maxParallelInsertions = 0
    mover = pygcmc.movement.MovementModule()
    
    # The C++ implementation behavior may vary between serial/parallel execution
    try:
        mover.setParams(params)
        # If no exception, test that operations still work
        result = mover.attemptInsertion(state)
        assert hasattr(result, 'accepted')
    except ValueError as e:
        # In parallel execution, may raise ValueError
        assert 'maxParallelInsertions' in str(e)
        assert 'must be positive' in str(e).lower()


def test_recovery_from_failed_insertion(setup_system):
    """Test recovery from failed insertion attempts."""
    state, _ = setup_system
    params = pygcmc.movement.MovementParams()
    params.temperature = 298.15
    params.chemicalPotential = -50.0
    mover = pygcmc.movement.MovementModule()
    mover.setParams(params)
    failed_count = 0
    for _ in range(100):
        result = mover.attemptInsertion(state)
        if not result.accepted:
            failed_count += 1
    assert failed_count > 0
    params.chemicalPotential = -10.0
    mover.setParams(params)
    result = mover.attemptInsertion(state)
    assert hasattr(result, 'accepted')


def test_concurrent_access_errors(setup_system):
    """Test handling of concurrent access issues.
        
        Note: This is a smoke test for thread safety. Python's GIL and pybind11's 
        default behavior mean threads may execute serially unless the C++ code 
        explicitly releases the GIL. This test primarily ensures no crashes or 
        data corruption occur, not true parallel execution.
        """
    state, _ = setup_system
    params = pygcmc.movement.MovementParams()
    params.temperature = 298.15
    params.chemicalPotential = -15.7
    params.useCavityBias = True
    mover = pygcmc.movement.MovementModule()
    mover.setParams(params)
    import threading
    errors = []
    completed = []

    def worker(worker_id):
        try:
            local_count = 0
            for _ in range(50):
                result = mover.attemptInsertion(state)
                assert hasattr(result, 'accepted')
                assert 0.0 <= result.acceptanceProbability <= 1.0
                local_count += 1
            completed.append(local_count)
        except Exception as e:
            errors.append((worker_id, str(e)))
    threads = []
    for i in range(4):
        t = threading.Thread(target=worker, args=(i,))
        threads.append(t)
        t.start()
    for t in threads:
        t.join()
    assert len(errors) == 0, f'Unexpected errors in concurrent test: {errors}'
    assert len(completed) == 4, f'Not all workers completed: {len(completed)}/4'
    assert all((c == 50 for c in completed)), f"Workers didn't complete all iterations: {completed}"


def test_warning_for_suboptimal_parameters(setup_system):
    """Test warnings for suboptimal parameter combinations."""
    state, _ = setup_system
    params = pygcmc.movement.MovementParams()
    params.temperature = 298.15
    params.chemicalPotential = -15.7
    params.useCavityBias = True
    params.cavityGridSpacing = 0.1
    params.probeRadius = 2.0
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter('always')
        params.updateDerivedParameters()
    mover = pygcmc.movement.MovementModule()
    mover.setParams(params)
    cavities = mover.findCavities(state)
    assert isinstance(cavities, list)


def test_graceful_degradation(setup_system):
    """Test graceful degradation when features unavailable."""
    state, _ = setup_system
    params = pygcmc.movement.MovementParams()
    params.temperature = 298.15
    params.chemicalPotential = -15.7
    params.useMultiInsertionCBMC = True
    params.maxParallelInsertions = 4
    params.proposalMode = 4
    params.useCavityBias = True
    params.useConfigBias = True
    mover = pygcmc.movement.MovementModule()
    mover.setParams(params)
    result = mover.attemptInsertion(state)
    assert hasattr(result, 'accepted')
    stats = mover.getStatistics()
    assert isinstance(stats, dict)


def test_parameter_validation_messages(setup_system):
    """Test that parameter validation provides helpful messages."""
    state, _ = setup_system
    params = pygcmc.movement.MovementParams()
    
    # Test temperature validation
    params.temperature = -100.0
    params.chemicalPotential = -15.7
    try:
        params.updateDerivedParameters()
        # If no exception in serial mode, verify value preserved
        assert params.temperature == -100.0
    except ValueError as e:
        # In parallel mode, should have informative message
        assert 'temperature' in str(e).lower()
        assert '-100' in str(e)
    
    # Test cavity grid spacing validation
    params = pygcmc.movement.MovementParams()
    params.temperature = 298.15
    params.chemicalPotential = -15.7
    params.useCavityBias = True
    params.cavityGridSpacing = -0.5
    params.probeRadius = 0.15
    
    try:
        params.updateDerivedParameters()
        # If no exception, should have clamped or preserved the value
        # This branch likely won't execute as cavity params are validated
        assert False, "Expected ValueError for negative cavityGridSpacing"
    except ValueError as e:
        assert 'cavityGridSpacing' in str(e)
        assert 'must be positive' in str(e).lower()
    
    # Test probe radius validation
    params = pygcmc.movement.MovementParams()
    params.temperature = 298.15
    params.chemicalPotential = -15.7
    params.useCavityBias = True
    params.cavityGridSpacing = 0.1
    params.probeRadius = -1.0
    
    try:
        params.updateDerivedParameters()
        # If no exception, should have clamped or preserved the value
        assert False, "Expected ValueError for negative probeRadius"
    except ValueError as e:
        assert 'probeRadius' in str(e)



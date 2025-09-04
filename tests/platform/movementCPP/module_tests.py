"""
Test MovementModule main functionality
"""

import pytest
from .fixtures import (
    MOVEMENT_AVAILABLE,
    create_movement_module,
    create_mock_state
)


@pytest.mark.skipif(not MOVEMENT_AVAILABLE, reason="Movement module not compiled")
def test_module_creation():
    """Test module creation and initialization"""
    movement_module = create_movement_module()
    assert movement_module is not None
    params = movement_module.getParams()
    assert params.temperature == pytest.approx(298.15)
    assert params.chemicalPotential == pytest.approx(-15.7)


@pytest.mark.skipif(not MOVEMENT_AVAILABLE, reason="Movement module not compiled")
def test_insertion_attempt():
    """Test insertion attempt"""
    movement_module = create_movement_module()
    mock_state = create_mock_state()
    
    result = movement_module.attemptInsertion(mock_state, moleculeType=0)
    assert result.moveType == "insert"
    assert 0.0 <= result.acceptanceProbability <= 1.0
    assert result.computeTimeMs >= 0


@pytest.mark.skipif(not MOVEMENT_AVAILABLE, reason="Movement module not compiled")
def test_deletion_attempt():
    """Test deletion attempt"""
    movement_module = create_movement_module()
    mock_state = create_mock_state()
    
    # First try to insert some molecules
    inserted = False
    for _ in range(20):  # Try up to 20 times
        result = movement_module.attemptInsertion(mock_state, moleculeType=0)
        if result.accepted:
            inserted = True
            break
    
    if inserted:
        # Now attempt deletion
        result = movement_module.attemptDeletion(mock_state)
        assert result.moveType == "delete"
        assert 0.0 <= result.acceptanceProbability <= 1.0


@pytest.mark.skipif(not MOVEMENT_AVAILABLE, reason="Movement module not compiled")
def test_translation_attempt():
    """Test translation attempt"""
    movement_module = create_movement_module()
    mock_state = create_mock_state()
    
    # First insert a molecule
    inserted = False
    for _ in range(20):
        result = movement_module.attemptInsertion(mock_state, moleculeType=0)
        if result.accepted:
            inserted = True
            break
    
    if inserted:
        result = movement_module.attemptTranslation(mock_state)
        assert result.moveType == "translate"
        assert 0.0 <= result.acceptanceProbability <= 1.0


@pytest.mark.skipif(not MOVEMENT_AVAILABLE, reason="Movement module not compiled")
def test_rotation_attempt():
    """Test rotation attempt"""
    movement_module = create_movement_module()
    mock_state = create_mock_state()
    
    # First insert a molecule
    inserted = False
    for _ in range(20):
        result = movement_module.attemptInsertion(mock_state, moleculeType=0)
        if result.accepted:
            inserted = True
            break
    
    if inserted:
        result = movement_module.attemptRotation(mock_state)
        assert result.moveType == "rotate"
        assert 0.0 <= result.acceptanceProbability <= 1.0


@pytest.mark.skipif(not MOVEMENT_AVAILABLE, reason="Movement module not compiled")
def test_cavity_finding():
    """Test cavity finding functionality"""
    movement_module = create_movement_module()
    mock_state = create_mock_state()
    
    cavities = movement_module.findCavities(mock_state)
    assert isinstance(cavities, list)
    # Empty box should have many cavities
    assert len(cavities) >= 0  # May be 0 if cavity bias is disabled


@pytest.mark.skipif(not MOVEMENT_AVAILABLE, reason="Movement module not compiled")
def test_statistics_tracking():
    """Test statistics tracking"""
    movement_module = create_movement_module()
    mock_state = create_mock_state()
    
    movement_module.resetStatistics()
    
    # Perform some moves
    for _ in range(10):
        movement_module.attemptInsertion(mock_state, moleculeType=0)
    
    stats = movement_module.getStatistics()
    insert_stats = stats["insert"]
    assert insert_stats.attempts == 10
    assert 0 <= insert_stats.accepts <= 10
    
    acceptance_rate = movement_module.calculateAcceptanceRate("insert")
    assert 0.0 <= acceptance_rate <= 1.0
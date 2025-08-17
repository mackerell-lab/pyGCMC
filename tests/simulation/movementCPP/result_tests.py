"""
Test MovementResult structure
"""

import pytest

try:
    import pygcmc
    MOVEMENT_AVAILABLE = hasattr(pygcmc, 'movement')
except ImportError:
    pygcmc = None
    MOVEMENT_AVAILABLE = False


@pytest.mark.skipif(not MOVEMENT_AVAILABLE, reason="Movement module not compiled")
def test_result_creation():
    """Test result creation and fields"""
    result = pygcmc.movement.MovementResult(
        accepted=True,
        energyChange=-5.0,
        acceptanceProbability=0.8,
        moveType="insert"
    )
    
    assert result.accepted == True
    assert result.energyChange == pytest.approx(-5.0)
    assert result.acceptanceProbability == pytest.approx(0.8)
    assert result.moveType == "insert"
    assert result.isSuccessful()


@pytest.mark.skipif(not MOVEMENT_AVAILABLE, reason="Movement module not compiled")
def test_result_summary():
    """Test result summary string"""
    result = pygcmc.movement.MovementResult(
        accepted=False,
        energyChange=10.0,
        acceptanceProbability=0.1,
        moveType="delete"
    )
    
    summary = result.summary()
    assert "delete" in summary
    assert "rejected" in summary or "accepted" in summary
    assert "10.0" in summary or "10.000" in summary
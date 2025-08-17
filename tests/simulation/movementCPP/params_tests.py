"""
Test MovementParams configuration
"""

import pytest

try:
    import pygcmc
    MOVEMENT_AVAILABLE = hasattr(pygcmc, 'movement')
except ImportError:
    pygcmc = None
    MOVEMENT_AVAILABLE = False


@pytest.mark.skipif(not MOVEMENT_AVAILABLE, reason="Movement module not compiled")
def test_default_params():
    """Test default parameter initialization"""
    params = pygcmc.movement.MovementParams()
    assert params.temperature == pytest.approx(298.15)
    assert params.chemicalPotential == pytest.approx(-15.7)
    assert params.useCavityBias == True
    assert params.useConfigBias == True


@pytest.mark.skipif(not MOVEMENT_AVAILABLE, reason="Movement module not compiled")
def test_custom_temperature():
    """Test parameter initialization with custom temperature"""
    params = pygcmc.movement.MovementParams(300.0)
    assert params.temperature == pytest.approx(300.0)
    assert params.beta == pytest.approx(1.0 / (8.314e-3 * 300.0))


@pytest.mark.skipif(not MOVEMENT_AVAILABLE, reason="Movement module not compiled")
def test_update_derived_parameters():
    """Test updating derived parameters after temperature change"""
    params = pygcmc.movement.MovementParams()
    params.temperature = 350.0
    params.updateDerivedParameters()
    assert params.beta == pytest.approx(1.0 / (8.314e-3 * 350.0))
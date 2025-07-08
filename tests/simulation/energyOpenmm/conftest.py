"""
Shared pytest fixtures for OpenMM energy tests.

This file contains fixtures that are automatically available to all tests
in the energyOpenmm directory.
"""

import warnings
import pytest
import pygcmc

# Suppress SWIG-related deprecation warnings
# These warnings come from the SWIG-generated Python bindings
warnings.filterwarnings("ignore", message="builtin type SwigPyPacked has no __module__ attribute", category=DeprecationWarning)
warnings.filterwarnings("ignore", message="builtin type SwigPyObject has no __module__ attribute", category=DeprecationWarning)
warnings.filterwarnings("ignore", message="builtin type swigvarlink has no __module__ attribute", category=DeprecationWarning)


def pytest_configure(config):
    """Configure pytest to ignore SWIG warnings."""
    config.addinivalue_line(
        "filterwarnings", "ignore:builtin type SwigPyPacked has no __module__ attribute:DeprecationWarning"
    )
    config.addinivalue_line(
        "filterwarnings", "ignore:builtin type SwigPyObject has no __module__ attribute:DeprecationWarning"
    )
    config.addinivalue_line(
        "filterwarnings", "ignore:builtin type swigvarlink has no __module__ attribute:DeprecationWarning"
    )


@pytest.fixture(autouse=True)
def clear_pme_state():
    """Clear PME engine state before and after each test to prevent cross-test contamination.
    
    This fixture is automatically applied to all tests in the energyOpenmm directory.
    It ensures that PME global state (parameters, grids, FFT weights) is cleared
    between tests to prevent failures when tests are run in parallel.
    
    This is necessary because PME uses global state in C++ that can cause
    segmentation faults when tests modify mesh sizes or other parameters.
    """
    # Clear before each test
    try:
        pygcmc.clearPMEEngine()
    except:
        pass  # No engine to clear on first run
    
    yield  # Run the test
    
    # Clear after each test
    try:
        pygcmc.clearPMEEngine()
    except:
        pass  # No engine to clear if test failed
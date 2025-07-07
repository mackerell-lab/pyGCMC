"""
Shared pytest fixtures for OpenMM energy tests.

This file contains fixtures that are automatically available to all tests
in the energyOpenmm directory.
"""

import pytest
import pygcmc


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
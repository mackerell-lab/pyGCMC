"""
Pytest configuration for PGP tests.

This file configures pytest behavior for PGP tests, particularly those that
need special handling due to global state issues.
"""

import pytest

def pytest_configure(config):
    """Register custom markers."""
    config.addinivalue_line(
        "markers", 
        "pgp_complete: marks tests as PGP Complete tests that modify global state"
    )

@pytest.fixture(autouse=False)  # Changed to False - tests will call reset themselves
def reset_pgp_state_before_test():
    """Reset PGP state before test if requested."""
    import pygcmc
    # Reset before test
    pygcmc.resetPGPState()
    yield
    # Could also reset after test if needed
    # pygcmc.resetPGPState()
"""
Configuration for OpenMM tests to avoid parallel execution issues.
"""

import pytest

def pytest_configure(config):
    """Configure pytest to mark OpenMM tests for serial execution."""
    config.addinivalue_line(
        "markers", "serial: mark test to run in serial mode (not parallel)"
    )

@pytest.fixture(scope="session", autouse=True)
def setup_openmm_tests():
    """Setup for OpenMM tests to avoid GC-related crashes."""
    import gc
    
    # Disable garbage collection during OpenMM tests to avoid segfaults
    # This is a workaround for the interaction between Python GC and OpenMM's C++ objects
    gc_was_enabled = gc.isenabled()
    if gc_was_enabled:
        gc.disable()
    
    yield
    
    # Re-enable GC after tests
    if gc_was_enabled:
        gc.enable()
        gc.collect()
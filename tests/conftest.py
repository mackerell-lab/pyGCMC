"""
Root conftest.py for all tests.

This file contains configurations that apply to all tests in the project.
"""

import warnings

# Suppress SWIG-related deprecation warnings globally
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
# tests/simulation/energyOpenmm/__init__.py
"""
Energy tests comparing PyGCMC with OpenMM
"""

from __future__ import annotations

import warnings


# OpenMM (SWIG) on Python 3.12 emits noisy DeprecationWarnings like:
# "builtin type SwigPyPacked has no __module__ attribute".
# Filter here so all OpenMM-dependent tests stay quiet without relying on pytest.ini.
warnings.filterwarnings(
    "ignore",
    category=DeprecationWarning,
    message=r"builtin type SwigPyPacked has no __module__ attribute",
)
warnings.filterwarnings(
    "ignore",
    category=DeprecationWarning,
    message=r"builtin type SwigPyObject has no __module__ attribute",
)
warnings.filterwarnings(
    "ignore",
    category=DeprecationWarning,
    message=r"builtin type swigvarlink has no __module__ attribute",
)

"""
GCMC Move Physics Tests - Main Entry Point

This file aggregates all move-level suites (acceptance rates, geometry checks,
region constraints). Run: pytest tests/platform/test_cpu_moves.py
"""

from cpu_moves.acceptance import TestMoveAcceptance
from cpu_moves.geometry import TestMoveGeometry
from cpu_moves.regions import (
    test_gcmc_region_box,
    test_gcmc_region_cylinder,
    test_gcmc_region_sphere,
)

__all__ = [
    "TestMoveAcceptance",
    "TestMoveGeometry",
    "test_gcmc_region_box",
    "test_gcmc_region_cylinder",
    "test_gcmc_region_sphere",
]

# tests/platform/test_movement_cpu.py
"""
CPU Movement Module Tests - Main Entry Point

This file imports all CPU-specific movement tests from movementCPU/ directory.
These tests cover CPU platform-specific functionality including:
- Quaternion rotation operations
- Periodic boundary conditions (PBC)
- Energy calculations and consistency
- NBFIX parameters and mixing rules
- State consistency and RNG reproducibility
- Performance benchmarks

Run: pytest tests/platform/test_movement_cpu.py
"""

import pytest
import pygcmc

# =====================================================
# Test classes from movementCPU/
# =====================================================

# Quaternion and rotation tests
from movementCPU.quaternion_rotation import TestQuaternionRotation
from movementCPU.quaternion_extended import TestQuaternionExtended
from movementCPU.rotation_uniformity import TestRotationUniformity

# PBC tests
from movementCPU.pbc_equivalence import TestPBCEquivalence
from movementCPU.noncubic_pbc import TestNonCubicPBC

# Energy and NBFIX tests
from movementCPU.energy_consistency import TestEnergyConsistency
from movementCPU.nbfix_mixing_rules import TestNBFIXMixingRules
from movementCPU.nbfix_extended import TestNBFIXExtended

# System state and configuration tests
from movementCPU.state_consistency import TestStateConsistency
from movementCPU.config_priority import (
    TestConfigPriority,
    TestStatisticsSampling,
    TestBoundaryConditions
)

# RNG reproducibility tests (functions)
from movementCPU.rng_reproducibility import (
    test_rng_seed_setting,
    test_seed_reset_behavior,
    test_different_seeds_different_results
)

# Performance and comprehensive tests
from movementCPU.performance_benchmark import TestPerformance
from movementCPU.comprehensive_improvements import TestComprehensiveImprovements
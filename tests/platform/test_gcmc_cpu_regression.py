"""
GCMC CPU Regression Tests - Main Entry Point

Collects all higher-level regression/feature/compatibility suites.
Run: pytest tests/platform/test_gcmc_cpu_regression.py
"""

from cpu_regression.legacy_inputs import TestGCMCCompatibility, TestGCMCStatistics
from cpu_regression.examples import TestGCMCCPUExamples
from cpu_regression.features import (
    TestWaterInsertion,
    TestCavityBias,
    TestRegionConstraints,
    TestTargetControl,
    TestParameterParsing,
    TestOutputValidation,
)
from cpu_regression.integration import TestGCMCIntegration
from cpu_regression.outputs_dat import TestOutputDAT

__all__ = [
    "TestGCMCCompatibility",
    "TestGCMCStatistics",
    "TestGCMCCPUExamples",
    "TestWaterInsertion",
    "TestCavityBias",
    "TestRegionConstraints",
    "TestTargetControl",
    "TestParameterParsing",
    "TestOutputValidation",
    "TestGCMCIntegration",
    "TestOutputDAT",
]

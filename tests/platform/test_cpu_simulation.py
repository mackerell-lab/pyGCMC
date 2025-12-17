"""
GCMC CPU Simulation Tests - Main Entry Point

Collects all end-to-end gcmc_cpu simulation/theory/regression suites that used
to live under tests/simulation/.

Run: pytest tests/platform/test_cpu_simulation.py
"""

from cpu_simulation.gcmc_benchmark import TestMemoryAndStability, TestPerformanceBenchmark
from cpu_simulation.gcmc_examples import TestGCMCExamples
from cpu_simulation.gcmc_functionality import (
    TestAdvancedFeatures,
    TestBasicGCMC,
    TestINPParameters,
    TestMcMoveProb,
    TestNumericalCorrectness,
    TestStatisticsOutput,
)
from cpu_simulation.gcmc_regression import TestContinuousIntegration, TestRegression
from cpu_simulation.gcmc_repro import TestSeedReproducibility
from cpu_simulation.gcmc_theory import TestGCMCTheory
from cpu_simulation.multicomponent_activity import TestMulticomponentActivity

__all__ = [
    "TestPerformanceBenchmark",
    "TestMemoryAndStability",
    "TestGCMCExamples",
    "TestMcMoveProb",
    "TestStatisticsOutput",
    "TestBasicGCMC",
    "TestINPParameters",
    "TestAdvancedFeatures",
    "TestNumericalCorrectness",
    "TestRegression",
    "TestContinuousIntegration",
    "TestSeedReproducibility",
    "TestGCMCTheory",
    "TestMulticomponentActivity",
]


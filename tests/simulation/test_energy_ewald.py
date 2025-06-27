# tests/simulation/test_energy_ewald_new.py
"""
Energy Ewald Tests - Main Entry Point

This file imports all Energy Ewald tests from modular sub-files.
Run: pytest tests/simulation/test_energy_ewald_new.py

Modular structure:
- basic_comparison.py: Basic energy comparison tests (3 functions)
- parameter_studies.py: Parameter studies tests (3 functions)
- advanced_analysis.py: Advanced analysis tests (2 functions)
- accuracy_validation.py: Accuracy validation tests (2 functions)
"""

# Basic energy comparison tests
from energyEwald.basic_comparison import (
    test_energy_methods_comparison,
    test_pbc_cutoff_ewald_comparison,
    test_ewald_symmetry
)

# Parameter studies tests
from energyEwald.parameter_studies import (
    test_ewald_parameter_sensitivity,
    test_ewald_error_convergence,
    test_ewald_error_tolerance
)

# Advanced analysis tests
from energyEwald.advanced_analysis import (
    test_madelung_constant,
    test_ewald_charge_neutrality
)

# Accuracy validation tests
from energyEwald.accuracy_validation import (
    test_ewald_exact,
    test_erfc_approx
)

# Support direct execution for testing
if __name__ == "__main__":
    import pytest
    pytest.main([__file__])

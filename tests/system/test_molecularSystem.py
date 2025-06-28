# tests/system/test_molecularSystem.py
"""
MolecularSystem Tests - Main Entry Point

This file imports all MolecularSystem tests from modular sub-files.
Run: pytest tests/system/test_molecularSystem.py

Modular structure (all files under 300 lines):
- Original test_molecularSystem.py (571 lines, 11 functions) split into:
  * molecularSystem/basic_combination.py: Basic combine tests (4 functions)
  * molecularSystem/error_handling.py: Error handling tests (3 functions)
  * molecularSystem/incompatible_tests.py: Incompatible file tests (2 functions)
  * molecularSystem/multiple_vs_single.py: Multiple vs single topology test (1 function)
  * molecularSystem/top_psf_comparison.py: TOP vs PSF comparison test (1 function)
- Helper fixtures and common code: molecularSystem/helpers.py

Total: 11 test functions across modular files (all modules under 280 lines each).
Original total: 571 lines → New structure: 595 lines (+4% for better organization)
"""

# Import fixtures from helpers
from molecularSystem.helpers import test_structure, test_topology

# Basic combination tests (4 functions)
from molecularSystem.basic_combination import (
    test_combine_structure_data,
    test_combine_topology_data,
    test_combine_bonds_and_angles,
    test_combine_mapping_data
)

# Error handling tests (3 functions)
from molecularSystem.error_handling import (
    test_combine_null_inputs,
    test_combine_mismatched_data,
    test_combine_empty_data
)

# Incompatible file tests (2 functions)
from molecularSystem.incompatible_tests import (
    test_combine_incompatible_files,
    test_combine_multiple_topologies
)

# Multiple vs single topology test (1 function)
from molecularSystem.multiple_vs_single import test_combine_multiple_vs_single

# TOP vs PSF comparison test (1 function)
from molecularSystem.top_psf_comparison import test_compare_top_psf_systems

# Support direct execution for testing
if __name__ == "__main__":
    import pytest
    pytest.main([__file__])
# tests/simulation/test_openmm.py
"""
OpenMM Tests - Main Entry Point

This file imports all OpenMM tests from modular sub-files.
Run: pytest tests/simulation/test_openmm.py

Modular structure:
- Original test_openmm_naive_nonbonded.py (554 lines, 5 functions) → energyOpenmm/naive_comparison.py
- Original test_openmm_nonbonded.py (2399 lines, 16 functions) split into:
  * energyOpenmm/simple_interactions.py: Basic interaction tests (4 functions)
  * energyOpenmm/intermediate_tests.py: Intermediate energy tests (2 functions)
  * energyOpenmm/method_comparisons.py: Method comparison tests (2 functions)
  * energyOpenmm/analysis_tests.py: Energy analysis tests (3 functions)
  * energyOpenmm/periodic_tests.py: Periodic boundary tests (5 functions)
- Original test_openmm_nonbonded_file.py (207 lines, 1 function) → energyOpenmm/file_based.py

Total: 22 test functions across modular files (all modules under 600 lines each).
Original total: 3160 lines → New structure: 4554 lines (+44% for better organization)
"""

# Naive nonbonded comparison tests (5 functions)
from energyOpenmm.naive_comparison import (
    test_openmm_energy_components,
    test_naive_energy_components,
    test_compare_openmm_naive_nonbonded,
    test_compare_pbc_energies,
    test_compare_cutoff_effects
)

# Simple interaction tests (4 functions)
from energyOpenmm.simple_interactions import (
    test_attractive_interaction,
    test_repulsive_interaction,
    test_pbc_interaction,
    test_cutoff_effect
)

# Intermediate energy tests (2 functions)
from energyOpenmm.intermediate_tests import (
    test_energy_symmetry,
    test_compare_custom_vs_standard_nonbonded
)

# Method comparison tests (2 functions)
from energyOpenmm.method_comparisons import (
    test_compare_nonbonded_methods,
    test_compare_switching_functions
)

# Energy analysis tests (3 functions)
from energyOpenmm.analysis_tests import (
    test_compare_separate_terms,
    test_compare_naive_vs_cutoff_energy,
    test_detailed_energy_comparison_simple_vs_openmm_cutoff
)

# Periodic boundary tests (5 functions)
from energyOpenmm.periodic_tests import (
    test_compare_all_methods_with_self_energy,
    test_analyze_openmm_energy_terms,
    test_cutoff_periodic_comparison,
    test_separate_lj_coulomb_periodic,
    test_compare_force_parameters_and_energies
)

# File-based test (1 function)
from energyOpenmm.file_based import (
    test_verify_openmm_expressions
)

# Support direct execution for testing
if __name__ == "__main__":
    import pytest
    pytest.main([__file__])
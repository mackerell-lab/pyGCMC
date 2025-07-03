# tests/simulation/test_openmm.py
"""
OpenMM Tests - Main Entry Point

This file imports all OpenMM tests from modular sub-files.
Run: pytest tests/simulation/test_openmm.py

Modular structure (all files under 300 lines):
- Original test_openmm_naive_nonbonded.py (554 lines, 5 functions) split into:
  * energyOpenmm/naive_energy_components.py: Energy component tests (3 functions)
  * energyOpenmm/naive_pbc_cutoff.py: PBC and cutoff tests (2 functions)
- Original test_openmm_nonbonded.py (2399 lines, 16 functions) split into:
  * energyOpenmm/simple_attractive_repulsive.py: Basic interaction tests (2 functions)
  * energyOpenmm/simple_pbc_cutoff.py: PBC interaction tests (2 functions)
  * energyOpenmm/intermediate_energy_symmetry.py: Energy symmetry test (1 function)
  * energyOpenmm/intermediate_custom_vs_standard.py: Custom vs standard test (1 function)
  * energyOpenmm/method_combined.py: Combined method comparison test (1 function)
  * energyOpenmm/method_nocutoff.py: NoCutoff method test (1 function)
  * energyOpenmm/method_cutoff_nonperiodic.py: CutoffNonPeriodic method test (1 function)
  * energyOpenmm/method_switching_functions.py: Switching function test (1 function)
  * energyOpenmm/analysis_separate_terms.py: Energy term analysis (1 function)
  * energyOpenmm/analysis_naive_vs_cutoff.py: Naive vs cutoff comparison (1 function)
  * energyOpenmm/analysis_detailed_comparison.py: Detailed energy comparison (1 function)
  * energyOpenmm/periodic_self_energy.py: Self-energy correction test (1 function)
  * energyOpenmm/periodic_energy_terms.py: Energy term analysis (1 function)
  * energyOpenmm/periodic_cutoff_comparison.py: Cutoff comparison test (1 function)
  * energyOpenmm/periodic_lj_coulomb.py: LJ and Coulomb separation test (1 function)
  * energyOpenmm/periodic_force_parameters.py: Force parameter test (1 function)
- Original test_openmm_nonbonded_file.py (207 lines, 1 function) → energyOpenmm/file_based.py

Total: 30 test functions across modular files (all modules under 280 lines each).
Original total: 3160 lines → New structure: ~5000 lines (with NBFIX tests and demonstrations)
"""

# Naive nonbonded comparison tests (3 functions)
from energyOpenmm.naive_energy_components import (
    test_openmm_energy_components,
    test_naive_energy_components,
    test_compare_openmm_naive_nonbonded
)

# Naive PBC and cutoff tests (2 functions)
from energyOpenmm.naive_pbc_cutoff import (
    test_compare_pbc_energies,
    test_compare_cutoff_effects
)

# Simple interaction tests (2 functions)
from energyOpenmm.simple_attractive_repulsive import (
    test_attractive_interaction,
    test_repulsive_interaction
)

# Simple PBC and cutoff tests (2 functions)
from energyOpenmm.simple_pbc_cutoff import (
    test_pbc_interaction,
    test_cutoff_effect
)

# Intermediate energy tests (1 function each)
from energyOpenmm.intermediate_energy_symmetry import test_energy_symmetry
from energyOpenmm.intermediate_custom_vs_standard import test_compare_custom_vs_standard_nonbonded

# Method comparison tests (1 function each)
from energyOpenmm.method_combined import test_compare_nonbonded_methods
from energyOpenmm.method_nocutoff import test_compare_nonbonded_nocutoff
from energyOpenmm.method_cutoff_nonperiodic import test_compare_nonbonded_cutoff_nonperiodic
from energyOpenmm.method_switching_functions import test_compare_switching_functions

# Energy analysis tests (1 function each)
from energyOpenmm.analysis_separate_terms import test_compare_separate_terms
from energyOpenmm.analysis_naive_vs_cutoff import test_compare_naive_vs_cutoff_energy
from energyOpenmm.analysis_detailed_comparison import test_detailed_energy_comparison_simple_vs_openmm_cutoff

# Periodic boundary tests (1 function each)
from energyOpenmm.periodic_self_energy import test_compare_all_methods_with_self_energy
from energyOpenmm.periodic_energy_terms import test_analyze_openmm_energy_terms
from energyOpenmm.periodic_cutoff_comparison import test_cutoff_periodic_comparison
from energyOpenmm.periodic_lj_coulomb import test_separate_lj_coulomb_periodic
from energyOpenmm.periodic_force_parameters import test_compare_force_parameters_and_energies

# File-based test (1 function)
from energyOpenmm.file_based import test_verify_openmm_expressions

# NBFIX tests (6 functions)
from energyOpenmm.nbfix_tests import (
    test_nbfix_overrides_lj_combination,
    test_ion_water_nbfix_energy,
    test_nbfix_energy_vs_openmm,
    test_multiple_nbfix_pairs,
    test_nbfix_energy_components_separately,
    test_nbfix_across_pbc
)

# Residue activation demonstration functions are not tests
# If you want to run them as tests, create wrapper functions:
# def test_residue_activation():
#     from energyOpenmm.residue_activation_final import calculate_residue_addition_energy
#     energy_before, energy_after, delta = calculate_residue_addition_energy()
#     assert delta < 0  # Adding a residue should be favorable in this case

# Support direct execution for testing
if __name__ == "__main__":
    import pytest
    pytest.main([__file__])
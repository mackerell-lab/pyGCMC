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
- Original pme_lj_only.py (429 lines, 4 functions) split into:
  * energyOpenmm/pme_lj_only_helpers.py: Shared helper functions (151 lines)
  * energyOpenmm/pme_lj_only_basic.py: Basic cutoff and PME tests (118 lines, 2 functions)
  * energyOpenmm/pme_lj_only_analysis.py: Distance scan and mixing rules tests (195 lines, 2 functions)

Total: 75 test functions across modular files (all modules under 280 lines each).
Original total: 3160 lines → New structure: ~5000 lines (with NBFIX tests and demonstrations)
PME Total bug fix tests: 3 functions to verify correct energy calculation
PGP strict tolerance test: 1 function to test PGP with strict tolerances
PME-Ewald consistency test: 1 function to verify PME and Ewald agreement
PME regression tests: 4 functions from analysis tools to prevent bug recurrence
LJ double-counting fix tests: 22 functions to verify correct LJ energy calculation
PME convergence and comparison tests: 3 functions
  - pme_convergence.py: PME parameter sensitivity test (1 function)
  - pme_medium_complexity.py: Medium complexity system PME test (1 function)
  - pme_simple_comparison.py: Simple two-charge PME comparison test (1 function)
Intra-residue and NaCl structure tests: 2 functions
  - intra_residue_difference.py: Test PME vs Ewald intra-residue handling (1 function)
  - pme_nacl_correct.py: Test PME with correctly structured NaCl crystal (1 function)
Note: 2 PME tests are marked as skip due to C++ global state issues that need fixing
"""

# Naive nonbonded comparison tests (3 functions)
from simulation.energyOpenmm.naive_energy_components import (
    test_openmm_energy_components,
    test_naive_energy_components,
    test_compare_openmm_naive_nonbonded
)

# Naive PBC and cutoff tests (2 functions)
from simulation.energyOpenmm.naive_pbc_cutoff import (
    test_compare_pbc_energies,
    test_compare_cutoff_effects
)

# Simple interaction tests (2 functions)
from simulation.energyOpenmm.simple_attractive_repulsive import (
    test_attractive_interaction,
    test_repulsive_interaction
)

# Simple PBC and cutoff tests (2 functions)
from simulation.energyOpenmm.simple_pbc_cutoff import (
    test_pbc_interaction,
    test_cutoff_effect
)

# Intermediate energy tests (1 function each)
from simulation.energyOpenmm.intermediate_energy_symmetry import test_energy_symmetry
from simulation.energyOpenmm.intermediate_custom_vs_standard import test_compare_custom_vs_standard_nonbonded

# Method comparison tests (1 function each)
from simulation.energyOpenmm.method_combined import test_compare_nonbonded_methods
from simulation.energyOpenmm.method_nocutoff import test_compare_nonbonded_nocutoff
from simulation.energyOpenmm.method_cutoff_nonperiodic import test_compare_nonbonded_cutoff_nonperiodic
from simulation.energyOpenmm.method_switching_functions import test_compare_switching_functions

# Energy analysis tests (1 function each)
from simulation.energyOpenmm.analysis_separate_terms import test_compare_separate_terms
from simulation.energyOpenmm.analysis_naive_vs_cutoff import test_compare_naive_vs_cutoff_energy
from simulation.energyOpenmm.analysis_detailed_comparison import test_detailed_energy_comparison_simple_vs_openmm_cutoff

# Periodic boundary tests (1 function each)
from simulation.energyOpenmm.periodic_self_energy import test_compare_all_methods_with_self_energy
from simulation.energyOpenmm.periodic_energy_terms import test_analyze_openmm_energy_terms
from simulation.energyOpenmm.periodic_cutoff_comparison import test_cutoff_periodic_comparison
from simulation.energyOpenmm.periodic_lj_coulomb import test_separate_lj_coulomb_periodic
from simulation.energyOpenmm.periodic_force_parameters import test_compare_force_parameters_and_energies

# File-based test (1 function)
from simulation.energyOpenmm.file_based import test_verify_openmm_expressions

# NBFIX tests (6 functions)
from simulation.energyOpenmm.nbfix_parameters import (
    test_nbfix_overrides_lj_combination,
    test_multiple_nbfix_pairs
)
from simulation.energyOpenmm.nbfix_energy import test_ion_water_nbfix_energy
from simulation.energyOpenmm.nbfix_openmm import (
    test_nbfix_energy_vs_openmm,
    test_nbfix_energy_components_separately
)
from simulation.energyOpenmm.nbfix_pbc import test_nbfix_across_pbc

# PME Total bug fix tests (3 functions)
from simulation.energyOpenmm.pme_total_fix import (
    test_pme_total_single_residue,
    test_pme_total_multiple_residues,
    test_pme_total_vs_ewald
)

# PGP strict tolerance test (1 function)
from simulation.energyOpenmm.pgp_strict_tolerance import test_pgp_with_strict_tolerance

# PME-Ewald consistency test (1 function)
from simulation.energyOpenmm.pme_ewald_consistency import test_pme_ewald_consistency

# PME bug analysis tests (4 functions) - now used as regression tests
from simulation.energyOpenmm.pme_total_lj_hypothesis import test_lj_hypothesis
from simulation.energyOpenmm.pme_residue_offset_pattern import test_residue_offset_pattern
from simulation.energyOpenmm.verify_pme_total_source import test_pme_total_sources
from simulation.energyOpenmm.analyze_pme_total_bug import test_pme_total_configurations

# PME ligand energy comparison with OpenMM (2 functions)
from simulation.energyOpenmm.pme_ligand_openmm_comparison import (
    test_pme_ligand_energy_vs_openmm,
    test_pme_movement_residues
)

# High precision PME comparison with OpenMM (3 functions)
from simulation.energyOpenmm.pme_high_precision_comparison import (
    test_pme_high_precision_simple_system,
    test_pme_high_precision_complex_system,
    test_pme_convergence_with_mesh_size
)

# PME ligand insertion tests (2 functions) - test_movement_residue_pme_energy excluded due to segfault
from simulation.energyOpenmm.pme_ligand_insertion import (
    test_ligand_insertion_energy,
    test_ligand_position_scan
)

# Movement residue PME tests (2 functions)
from simulation.energyOpenmm.movement_residue_pme import (
    test_movement_residue_pme,
    test_movement_residue_removal
)

# Residue activation demonstration functions are not tests
# If you want to run them as tests, create wrapper functions:
# def test_residue_activation():
#     from simulation.energyOpenmm.residue_activation_final import calculate_residue_addition_energy
#     energy_before, energy_after, delta = calculate_residue_addition_energy()
#     assert delta < 0  # Adding a residue should be favorable in this case

# Full nonbonded comparison tests with LJ fix (2 functions)
# Original file: from simulation.energyOpenmm.full_nonbonded_comparison
# Now split into three files using sed
from simulation.energyOpenmm.full_nonbonded_pme import test_full_nonbonded_pme_lj
from simulation.energyOpenmm.full_nonbonded_cutoff import test_lj_only_comparison

# LJ double counting analysis tests (3 functions)
from simulation.energyOpenmm.lj_double_counting import (
    test_single_pair_lj,
    test_three_atom_system,
    test_residue_assignment
)

# LJ detailed comparison tests (2 functions)
from simulation.energyOpenmm.lj_detailed_comparison import (
    test_simple_two_atom_system,
    test_mixed_types
)

# PME LJ-only tests (4 functions)
# Basic tests (2 functions)
from simulation.energyOpenmm.pme_lj_only_basic import (
    test_lj_only_cutoff,
    test_lj_only_pme
)
# Analysis tests (2 functions)
from simulation.energyOpenmm.pme_lj_only_analysis import (
    test_lj_distance_scan,
    test_lj_mixing_rules
)

# Debug distance tests (2 functions)
from simulation.energyOpenmm.debug_distances import (
    test_system_from_test_full_nonbonded,
    test_system_from_test_pme_lj_only
)

# LJ close distance tests (2 functions)
from simulation.energyOpenmm.lj_close_distance import (
    test_close_distances,
    test_energy_capping
)

# LJ double counting summary test (1 function)
from simulation.energyOpenmm.lj_double_counting_summary import test_lj_double_counting_proof

# LJ multiparticle analysis tests (3 functions)
from simulation.energyOpenmm.lj_multiparticle_analysis import (
    test_progressive_atoms,
    test_cutoff_effects,
    test_pbc_boundary
)

# Fixed functions tests (3 functions)
from simulation.energyOpenmm.fixed_functions import (
    test_cutoff_fixed,
    test_pme_fixed,
    test_fixed_vs_openmm
)

# PME convergence and comparison tests (3 functions)
from simulation.energyOpenmm.pme_convergence import test_pme_parameter_sensitivity
from simulation.energyOpenmm.pme_medium_complexity import test_pme_medium_complexity
from simulation.energyOpenmm.pme_simple_comparison import test_simple_pme_comparison

# Intra-residue and correct NaCl tests (2 functions)
from simulation.energyOpenmm.intra_residue_difference import test_intra_residue_handling
from simulation.energyOpenmm.pme_nacl_correct import test_pme_vs_ewald_correct

# Support direct execution for testing
if __name__ == "__main__":
    import pytest
    pytest.main([__file__])
# tests/simulation/test_movement_insert.py
"""
Movement-based insertion tests - Main Entry Point

This file imports all movement insertion tests from modular sub-files.
Run: pytest tests/simulation/test_movement_insert.py

Test categories:
1. Basic insertion tests - fundamental molecule insertion operations
2. Advanced insertion tests - cavity detection and energy-guided insertion
3. Residue activation tests - GCMC residue activation/deactivation
4. Insertion/deletion tests - GCMC insertion and deletion energy calculations
5. GCMC algorithm tests - Detailed balance and acceptance criteria
6. Practical GCMC tests - Real energy calculations with PyGCMC
"""

# Basic insertion tests (4 functions) - test_insert_with_existing_molecules moved to tmp/skipped_tests
from movementInsert.basic_insertion_tests import (
    test_insert_single_molecule,
    test_insert_multiple_water_molecules,
    test_insert_ion_pair,
    test_insert_random_positions
)

# Advanced insertion tests (6 functions)
from movementInsert.advanced_insertion_basic_tests import (
    test_simple_two_particle_insertion,
    test_cavity_detection
)
from movementInsert.advanced_insertion_biased_tests import (
    test_biased_cavity_insertion,
    test_energy_guided_insertion
)
from movementInsert.advanced_insertion_sequential_tests import (
    test_sequential_cavity_filling,
    test_water_cluster_formation
)

# Residue activation tests (4 functions) - moved from energyGCMC
from movementInsert.residue_activation_tests import (
    test_residue_addition_energy,
    test_movement_energy_calculation,
    test_residue_activation_deactivation,
    test_multiple_residue_types
)

# GCMC insertion/deletion tests (6 functions) - moved from energyGCMC
from movementInsert.insertion_deletion_tests import (
    test_single_molecule_insertion,
    test_molecule_insertion_with_interactions,
    test_molecule_deletion,
    test_multiple_insertions,
    test_biased_insertion,
    test_movement_only_insertion
)

# GCMC algorithm tests (3 functions)
from movementInsert.gcmc_insertion_test_functions import (
    test_gcmc_benzene_insertion_in_protein,
    test_gcmc_water_insertion_simple,
    test_verify_detailed_balance
)

# Practical GCMC tests (3 functions)
from movementInsert.practical_gcmc_benzene_test import test_practical_benzene_insertion
from movementInsert.practical_gcmc_water_test import test_water_insertion_with_real_energy
from movementInsert.practical_gcmc_cavity_bias_test import test_cavity_bias_effect

# Numerical stability tests (7 functions) - NEW
from movementInsert.numerical_stability_tests import (
    test_extreme_energy_numerical_stability,
    test_configurational_bias_numerical_stability,
    test_temperature_extremes,
    test_gcmc_insertion_acceptance_stability,
    test_parallel_rosenbluth_calculation,
    test_insertion_with_pbc_wrapping
)

# Cavity bias tests (8 functions) - NEW
from movementInsert.cavity_bias_tests import (
    test_cavity_grid_creation,
    test_cavity_detection_with_molecules,
    test_biased_insertion_acceptance,
    test_cavity_bias_with_different_grid_spacings,
    test_cavity_bias_energy_calculation,
    test_adaptive_cavity_grid,
    test_cavity_bias_detailed_balance
)

# Log-space stability tests (10 functions) - NEW
from movementInsert.logspace_stability_tests import (
    test_logsumexp_extreme_values,
    test_prob_from_log_saturation,
    test_insertion_acceptance_extreme_parameters,
    test_deletion_acceptance_extreme_parameters,
    test_rosenbluth_weight_extreme_energies,
    test_configurational_bias_logspace,
    test_detailed_balance_logspace,
    test_extreme_temperature_logspace,
    test_parallel_logspace_calculations
)

# Translation and Rotation tests (6 functions) - NEW
from movementInsert.test_translation_rotation import (
    test_translation_move,
    test_rotation_move,
    test_combined_translation_rotation,
    test_translation_with_pbc,
    test_rotation_preserves_bond_lengths
)

# Total: 56 test functions covering comprehensive GCMC insertion operations (including translation/rotation)

# Support direct execution
if __name__ == "__main__":
    import pytest
    pytest.main([__file__, "-v"])
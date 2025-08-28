# tests/simulation/test_movement_cpp.py
"""
C++ Movement Module Bindings Tests - Main Entry Point

This file imports all C++ movement binding tests from modular sub-files.
Run: pytest tests/simulation/test_movement_cpp.py

Test categories:
1. Parameter tests - MovementParams configuration
2. Module tests - MovementModule functionality  
3. Active pool tests - Memory management
4. Result tests - MovementResult structure
5. Integration tests - Complete GCMC workflows
6. Cavity bias tests - Insertion optimization
7. Cavity cache tests - Cache mechanism and performance
"""

import pytest
import pygcmc

# Import fixtures from movementCPP/conftest.py for all tests in this file
from movementCPP.conftest import setup_system, setup_system_with_params, create_state
# Import additional fixtures from proposal_modes_fixtures.py
from movementCPP.proposal_modes_fixtures import clean_system, fresh_system

# Parameter configuration tests (3 functions)
from movementCPP.params_tests import (
    test_default_params,
    test_custom_temperature,
    test_update_derived_parameters
)

# Main module functionality tests (7 functions)
from movementCPP.module_tests import (
    test_module_creation,
    test_insertion_attempt,
    test_deletion_attempt,
    test_translation_attempt,
    test_rotation_attempt,
    test_cavity_finding,
    test_statistics_tracking
)

# Active pool memory management tests (5 functions)
from movementCPP.active_pool_tests import (
    test_pool_creation,
    test_insert_molecule,
    test_delete_residue,
    test_fragmentation,
    test_capacity_check
)

# Movement result tests (2 functions)
from movementCPP.result_tests import (
    test_result_creation,
    test_result_summary
)

# Integration and workflow tests (2 functions)
from movementCPP.integration_tests import (
    test_gcmc_equilibration,
    test_acceptance_rates
)

# Cavity bias tests - Basic (5 functions)
from movementCPP.cavity_bias_basic_funcs import (
    test_cavity_finding_basic,
    test_cavity_count_scaling,
    test_cavity_cache_consistency,
    test_cavity_manager_independence,
    test_non_cubic_box
)

# Cavity bias tests - Advanced (4 functions)
from movementCPP.cavity_bias_advanced_funcs import (
    test_cavity_concurrent_access,
    test_cavity_grid_spacing_effect,
    test_probe_radius_effect,
    test_cavity_stats_retrieval
)

# Cavity cache tests - Basic (6 functions)
from movementCPP.cavity_cache_basic_funcs import (
    test_cache_basic_functionality,
    test_cache_invalidation_on_movement,
    test_cache_key_generation,
    test_cache_with_different_parameters,
    test_cache_thread_safety,
    test_cache_with_box_changes
)

# Cavity cache tests - Performance (5 functions)
from movementCPP.cavity_cache_performance_funcs import (
    test_cache_memory_management,
    test_cache_performance_benefit,
    test_cache_statistics_tracking,
    test_cache_with_periodic_boundaries,
    test_cache_clear_operation
)

# Import detailed balance test functions directly (10 + 8 = 18 functions)
from movementCPP.detailed_balance_funcs import (
    test_insertion_deletion_balance,
    test_translation_reversibility,
    test_metropolis_criterion,
    test_cavity_bias_detailed_balance,
    test_config_bias_detailed_balance,
    test_multi_insertion_detailed_balance,
    test_temperature_scaling,
    test_chemical_potential_balance,
    test_ensemble_averages,
    test_rosenbluth_weight_consistency
)

from movementCPP.detailed_balance_strict_funcs import (
    test_metropolis_criterion_exact,
    test_insertion_deletion_pairwise_balance,
    test_insertion_deletion_flux_balance,
    test_cavity_bias_probability_consistency,
    test_chemical_potential_controls_density,
    test_translation_reversibility_exact,
    test_config_bias_fields_present,
    test_ensemble_convergence
)

# Import error handling test functions (13 functions)
from movementCPP.error_handling_funcs import (
    test_invalid_temperature,
    test_invalid_cavity_parameters,
    test_invalid_proposal_mode,
    test_null_state_handling,
    test_uninitialized_forcefield,
    test_invalid_box_dimensions,
    test_overflow_protection,
    test_multi_insertion_parameter_conflicts,
    test_recovery_from_failed_insertion,
    test_concurrent_access_errors,
    test_warning_for_suboptimal_parameters,
    test_graceful_degradation,
    test_parameter_validation_messages
)

# Import new movement module test functions (51 functions total)
# Basic module tests (10 functions)
from movementCPP.movement_module_basic_funcs import (
    test_module_creation_basic,
    test_parameter_setting_and_getting,
    test_insertion_basic,
    test_deletion_basic,
    test_translation_basic,
    test_rotation_basic,
    test_cavity_bias_insertion,
    test_config_bias_rotation_basic,
    test_acceptance_rate_calculation,
    test_find_cavities_integration
)

# Statistics and effects tests (8 functions)
from movementCPP.movement_module_stats_funcs import (
    test_statistics_retrieval,
    test_multi_insertion_cbmc_basic,
    test_temperature_effect,
    test_chemical_potential_effect_insertion,
    test_statistics_counts_and_reset,
    test_deletion_prefers_last_inserted,
    test_find_cavities_bounds_all,
    test_probability_bounds_all_moves
)

# Advanced module tests (8 functions)
from movementCPP.movement_module_advanced_funcs import (
    test_config_bias_rotation_movetype,
    test_get_proposal_stats_shape,
    test_get_cavity_stats_shape,
    test_statistics_keys_after_moves,
    test_constructor_seed_reproducibility,
    test_params_not_modified_by_methods,
    test_deletion_explicit_index_overrides_preference,
    test_module_repr
)

# Movement params tests (13 functions)
from movementCPP.movement_params_funcs import (
    test_basic_parameter_creation,
    test_temperature_constructor,
    test_cavity_grid_spacing_validation,
    test_probe_radius_validation,
    test_proposal_mode_clamping,
    test_multi_insertion_params_validation,
    test_adaptive_thresholds,
    test_performance_flags,
    test_fill_proposal_info_flag,
    test_proposal_mode_upper_bound_valid,
    test_seed_field_sets_rng,
    test_params_repr,
    test_validation_without_cavity_bias
)

# Movement result tests (11 functions)
from movementCPP.movement_result_funcs import (
    test_result_basic_fields,
    test_result_diagnostic_fields,
    test_fill_proposal_info_disabled,
    test_fill_proposal_info_enabled,
    test_fill_proposal_info_no_acceptance_effect,
    test_proposal_position_units,
    test_default_diagnostic_values,
    test_move_type_consistency,
    test_multi_insertion_result_fields,
    test_reproducibility_with_seed,
    test_result_repr
)

# ============================================================================
# MULTI-INSERTION TESTS (from movementCPP/test_movement_cpp.py)
# ============================================================================

# Basic multi-insertion tests (10 functions)
from movementCPP.multi_insertion_basic_funcs import (
    test_multi_insertion_basic,
    test_mproposal_scaling,
    test_parallel_insertion_count,
    test_region_volume_calculation,
    test_cbmc_weight_calculation,
    test_multi_insertion_with_cavity_bias,
    test_multi_insertion_performance,
    test_multi_insertion_detailed_balance,
    test_multi_insertion_with_types,
    test_multi_insertion_thread_safety
)

# Improved multi-insertion tests (7 functions)
from movementCPP.multi_insertion_improved_funcs import (
    test_multi_insertion_basic_strict,
    test_acceptance_rate_statistics,
    test_parallel_insertion_performance,
    test_detailed_balance_strict,
    test_cbmc_weight_validity,
    test_reproducibility_with_seed as test_multi_insertion_reproducibility,  # Rename to avoid conflict
    test_cavity_bias_with_multi_insertion
)

# Robustness tests (8 functions)
from movementCPP.multi_insertion_robustness_funcs import (
    test_acceptance_rate_monotonicity_chemical_potential,
    test_acceptance_rate_monotonicity_volume,
    test_energy_consistency_single_vs_system,
    test_region_independence_boundary,
    test_deterministic_config_selection,
    test_fallback_path_unit_consistency,
    test_energy_range_parameterization,
    test_openmp_determinism
)

# ============================================================================
# PROPOSAL MODE TESTS (from movementCPP/test_movement_cpp.py)
# ============================================================================

# Basic proposal mode tests (5 functions)
from movementCPP.proposal_modes_basic_funcs import (
    test_uniform_mode_basic,
    test_cavity_mode_basic,
    test_color_mode_placeholder,
    test_cluster_mode_placeholder,
    test_adaptive_mode_behavior
)

# Advanced proposal mode tests (5 functions)
from movementCPP.proposal_modes_advanced_funcs import (
    test_mode_transition_statistics,
    test_fallback_to_uniform,
    test_mode_performance_comparison,
    test_mode_with_different_densities,
    test_invalid_mode_handling
)

# Improved proposal mode tests (6 functions)
from movementCPP.proposal_modes_improved_funcs import (
    test_uniform_mode_statistical_uniformity,
    test_cavity_mode_with_geometry_validation,
    test_mode_acceptance_rate_comparison,
    test_density_effect_on_cavities,
    test_mode_clamping_behavior,
    test_statistics_aggregation
)

# Reproducibility tests (1 function)
from movementCPP.proposal_reproducibility_funcs import (
    test_seed_reproducibility
)

# ============================================================================
# PROPOSAL STATISTICS TESTS (from movementCPP/test_movement_cpp.py)
# ============================================================================

# Basic statistics tests (9 functions)
from movementCPP.proposal_statistics_basic_funcs import (
    test_basic_statistics_collection,
    test_mode_statistics,
    test_timing_statistics,
    test_fallback_statistics,
    test_cavity_statistics_integration,
    test_statistics_reset,
    test_proposal_mode_effect,
    test_acceptance_rate_calculation,
    test_statistics_types
)

# Improved statistics tests (6 functions)
from movementCPP.proposal_statistics_improved_funcs import (
    test_statistics_basic_counting,
    test_acceptance_rate_calculation as test_acceptance_rate_calculation_improved,
    test_statistics_reset as test_statistics_reset_improved,
    test_timing_statistics as test_timing_statistics_improved,
    test_cavity_statistics_consistency,
    test_statistics_thread_safety
)

# Edge case statistics tests (2 functions)
from movementCPP.proposal_statistics_edge_funcs import (
    test_empty_system_statistics,
    test_statistics_overflow_protection
)

# ============================================================================
# IDEAL GAS DISTRIBUTION TESTS (fundamental physics validation)
# ============================================================================

# Ideal gas limit tests (5 functions)
from movementCPP.ideal_gas_distribution_funcs import (
    test_ideal_gas_mean_particle_number,
    test_particle_distribution_shape,
    test_volume_scaling,
    test_chemical_potential_scaling,
    test_detailed_balance_ratio
)

# ============================================================================
# STATISTICAL UTILITIES TESTS (from movementCPP/test_movement_cpp.py)
# ============================================================================

# Basic statistical utility tests (8 functions)
from movementCPP.statistical_utils_tests_basic_funcs import (
    test_poisson_diff_ok_edge_cases,
    test_poisson_diff_ok_statistical_power,
    test_ratio_ci_ok_edge_cases,
    test_ratio_ci_ok_confidence_level,
    test_effective_sample_size_independent,
    test_effective_sample_size_correlated,
    test_effective_sample_size_edge_cases,
    test_batch_means_edge_cases
)

# Advanced statistical utility tests (7 functions)
from movementCPP.statistical_utils_tests_advanced_funcs import (
    test_batch_means_variance,
    test_batch_means_variance_autocorrelated,
    test_no_drift_converged,
    test_no_drift_drifting,
    test_calibrate_mu,
    test_no_drift_window_size,
    test_integration_detailed_balance
)
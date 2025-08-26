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
import numpy as np

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

# Common fixtures for all tests
@pytest.fixture
def setup_system():
    """Setup test system with known parameters for detailed balance tests."""
    state = pygcmc.MCState()
    state.info.box = np.array([4.0, 4.0, 4.0])
    
    # Setup force field with known interactions
    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.ljEps = [1.0]  # kJ/mol
    ff.ljSigma = [0.3]  # nm
    state.forcefield = ff
    
    params = pygcmc.movement.MovementParams()
    params.temperature = 300.0  # K
    params.chemicalPotential = -15.0  # kJ/mol - higher for better insertion rate
    params.seed = 42
    
    return state, params

@pytest.fixture
def create_state():
    """Create a test MCState with given box size."""
    def _create(box_nm=5.0):
        state = pygcmc.MCState()
        if isinstance(box_nm, (list, tuple)):
            state.info.box = np.array(box_nm)
        else:
            state.info.box = np.array([box_nm, box_nm, box_nm])
        
        # Setup force field
        ff = pygcmc.MCForceField()
        ff.numTotalTypes = 1
        ff.numMovementTypes = 1
        ff.ljEps = [0.5]
        ff.ljSigma = [0.3]
        state.forcefield = ff
        
        return state
    return _create


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
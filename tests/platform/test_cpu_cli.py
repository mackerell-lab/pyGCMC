"""
GCMC CLI Tests - Main Entry Point

This file imports all CLI tests from modular sub-files.
Run: pytest tests/platform/test_cpu_cli.py

Test categories:
1. Basic tests - Help, missing input, basic run
2. Output tests - Verbose output, file generation
3. Validation tests - Deterministic seed, invalid input, parameter validation
4. Compatibility tests - gcmc_gpu INP/units behaviors
5. Regression tests - builder/geometry edge cases
"""

# Import fixtures from cpu_cli/conftest.py for all tests in this file
from cpu_cli.conftest import (
    gcmc_cpu,
    test_data_dir,
    temp_dir
)

# Basic CLI tests (3 functions)
from cpu_cli.basic_tests import (
    test_gcmc_cpu_help,
    test_gcmc_cpu_missing_inp,
    test_gcmc_cpu_basic_run
)

# Output and file generation tests (3 functions)
from cpu_cli.output_tests import (
    test_gcmc_cpu_verbose_output,
    test_gcmc_cpu_output_files,
    test_gcmc_cpu_logging_regression
)

# Validation and error handling tests (3 functions)
from cpu_cli.validation_tests import (
    test_gcmc_cpu_deterministic_seed,
    test_gcmc_cpu_invalid_inp,
    test_gcmc_cpu_parameter_validation
)

# gcmc_gpu-style INP compatibility and unit semantics tests (10 functions + param cases)
from cpu_cli.inp_units_compat import (
    test_inp_units_gcmc_gpu_box_roundtrip_cryst1,
    test_inp_units_gcmc_gpu_fragmuex_kcal_to_kj_log,
    test_fragmuex_scales_activity_and_acceptance_in_ideal_gas_limit,
    test_dump_accept_log_consistent_with_final_pdb_counts,
    test_inp_random_seed_used_when_cli_missing,
    test_inp_units_gcmc_gpu_converts_grid_dx_and_cutoffs_and_target_volume,
    test_inp_units_gcmc_gpu_gcmc_region_numeric_conversion_affects_volume,
    test_nbar_volume_based_mode_sets_activity_from_concentration,
    test_nbar_const_water_nbar_scales_all_fragments_by_fragconc,
    test_nbar_number_water_nbar_updates_activity_after_first_insertion,
    test_gcmc_gpu_active_and_muex_files_are_written,
)

# gcmc_gpu/opencl example deck smoke tests (4 functions)
from cpu_cli.gcmc_gpu_inp_examples import (
    test_gcmc_gpu_water_tip3p_runs_and_honors_op_outputs,
    test_gcmc_gpu_version_gcmc_2_0_implies_angstrom_units_and_op_outputs,
    test_gcmc_gpu_multi_salt_runs_and_converts_muex_list,
    test_gcmc_opencl_style_alias_keys_are_honored,
)

# Initial system builder regressions and bias-term consistency (4 functions)
from cpu_cli.initial_system_loading import (
    test_builder_initial_system_loaded_with_itp_par,
    test_insertion_does_not_overwrite_initial_residues,
    test_cavity_bias_fraction_matches_simple_grid_occupancy,
    test_cbmc_and_cavity_terms_reconstruct_insertion_pacc,
)

# Fragment template geometry regression (1 function)
from cpu_cli.fragment_template_geometry import (
    test_fragment_template_geometry_from_monomerdir,
)

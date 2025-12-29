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
    test_gcmc_cpu_parameter_validation,
    test_dump_params_reports_unknown_inp_keys,
    test_dump_params_reports_ignored_inp_keys,
    test_strict_inp_keys_fails_on_unknown_or_ignored_keys,
    test_strict_inp_keys_passes_when_no_unknown_or_ignored_keys,
)

# Unit/geometry heuristic warnings and strict mode (9 functions)
from cpu_cli.inp_units_heuristic_warnings import (
    test_nm_mode_suspicious_lengths_emit_unit_warnings,
    test_strict_inp_fails_on_unit_warnings,
    test_pdb_cryst1_vs_inp_box_size_factor10_mismatch_is_reported,
    test_nm_mode_probe_radius_large_emits_heuristic_warning_and_converts_value,
    test_nm_mode_grid_spacing_large_emits_warning,
    test_nm_mode_gcmc_region_sphere_large_emits_warning,
    test_target_volume_box_mismatch_emits_warning,
    test_pairlist_cutoff_relation_emits_warning,
    test_reasonable_nm_deck_has_no_warnings,
)

# gcmc_gpu-style INP compatibility and unit semantics tests (11 functions + param cases)
from cpu_cli.inp_units_compat import (
    test_inp_units_gcmc_gpu_box_roundtrip_cryst1,
    test_inp_units_gcmc_gpu_fragmuex_kcal_to_kj_log,
    test_fragmuex_scales_activity_and_acceptance_in_ideal_gas_limit,
    test_dump_accept_log_consistent_with_final_pdb_counts,
    test_inp_random_seed_used_when_cli_missing,
    test_inp_units_auto_defaults_to_gcmc_gpu_angstrom_for_inp_files,
    test_inp_units_gcmc_gpu_converts_grid_dx_and_cutoffs_and_target_volume,
    test_inp_units_gcmc_gpu_gcmc_region_numeric_conversion_affects_volume,
    test_nbar_volume_based_mode_sets_activity_from_concentration,
    test_nbar_const_water_nbar_scales_all_fragments_by_fragconc,
    test_nbar_const_water_nbar_scales_multiple_fragments_in_one_run,
    test_nbar_number_water_nbar_updates_activity_after_first_insertion,
    test_gcmc_gpu_active_and_muex_files_are_written,
)

# gcmc_gpu/opencl example deck smoke tests (4 functions)
from cpu_cli.gcmc_gpu_inp_examples import (
    test_gcmc_gpu_inp_key_inventory_matches_templates,
    test_gcmc_gpu_water_tip3p_runs_and_honors_op_outputs,
    test_gcmc_gpu_version_gcmc_2_0_implies_angstrom_units_and_op_outputs,
    test_gcmc_gpu_multi_salt_runs_and_converts_muex_list,
    test_gcmc_opencl_style_alias_keys_are_honored,
)

# MH proposal-bias regressions (1 function)
from cpu_cli.asymmetric_proposal_bias import (
    test_asymmetric_attempt_prob_enters_mh_proposal_ratio_and_pacc,
)

# gcmc_opencl example deck smoke tests (3 functions)
from cpu_cli.gcmc_opencl_examples import (
    test_opencl_inp_key_inventory_matches_examples,
    test_opencl_twowater_example_runs_and_preserves_pdb_cryst1,
    test_opencl_waterbox_hollow_example_smoke_runs,
    test_opencl_test_mg_example_smoke_runs_and_converts_cutoff_units,
    test_opencl_waterbox_example_smoke_runs_and_enables_cavity_bias,
    test_opencl_benz_example_smoke_runs_and_parses_multi_fragment_lists,
    test_opencl_protein_example_dump_params_parses_complex_deck,
    test_opencl_protein_example_active_muex_outputs,
    test_opencl_cdk2_example_smoke_runs_and_reports_ignored_gcmc_cutoff,
    test_opencl_lysozyme_example_smoke_runs_and_reports_ignored_map_keys,
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

# Movement acceptance + energy closure for TRANSLATE/ROTATE (2 functions)
from cpu_cli.movement_energy_translation_rotation import (
    test_translation_move_deltaU_matches_analytic_lj_and_is_logged_in_dump_accept,
    test_rotation_move_deltaU_matches_analytic_lj_for_multi_atom_fragment,
    test_translation_move_deltaU_matches_energy_backend_components_with_coulomb,
    test_rotation_move_deltaU_matches_energy_backend_components_with_coulomb,
)

# Physical-contract GCMC checks (3 functions)
from cpu_cli.physical_contracts import (
    test_acceptance_formula_insertion_and_deletion_ideal_gas,
    test_acceptance_log_pacc_matches_formula_all_moves,
    test_cavity_bias_pacc_matches_acceptance_formula,
    test_cbmc_cavity_bias_acceptance_closure,
    test_mu_targets_stable_concentration_with_insertion_and_deletion,
    test_poisson_number_distribution_ideal_gas,
    test_cutoff_excludes_far_lj_interactions_in_insertion_energy,
    test_insertion_deltaU_matches_analytic_coulomb_energy,
    test_insertion_deltaU_matches_analytic_lj_plus_coulomb_energy,
    test_deletion_deltaU_is_negative_of_insertion_for_charged_system,
)

# Multi-fragment physical contracts (2 functions)
from cpu_cli.physical_contracts_multifragment import (
    test_multifragment_cross_species_deltaU_matches_energy_components,
    test_multifragment_poisson_distribution_ideal_gas,
)

# CBMC + region contract regression (1 function)
from cpu_cli.cbmc_region_contract import (
    test_cbmc_insertion_respects_gcmc_region,
)

# CBMC trial-energy observability and Rosenbluth closure (1 function)
from cpu_cli.cbmc_trial_energies import (
    test_dump_accept_exposes_cbmc_trial_energies_and_closes_rosenbluth_terms,
)

# CBMC trial-count key compatibility (1 function)
from cpu_cli.cbmc_trial_count_alias import (
    test_num_conf_bias_trial_sets_cbmc_trials_when_fragconf_missing,
)

# CBMC trial-count visibility + strictness (3 functions)
from cpu_cli.cbmc_trials_visibility import (
    test_dump_params_exposes_conf_bias_trials_per_fragment,
    test_dump_params_uses_num_conf_bias_trial_when_fragconf_missing,
    test_strict_inp_fails_when_use_conf_bias_has_no_trial_count_source,
)

# Cavity-bias knob regressions (2 functions)
from cpu_cli.cavity_bias_knobs import (
    test_cavity_bias_probe_radius_changes_cavity_fraction,
    test_cavity_bias_exclude_hydrogens_changes_cavity_fraction,
    test_cavity_bias_use_vdw_radius_changes_cavity_fraction,
)

# MaxCount hard-cap policy regression (1 function)
from cpu_cli.maxcount_policy import (
    test_maxcount_is_initial_plus_mcsteps_plus_buffer,
)

# Movement step-size parameter compatibility (1 function)
from cpu_cli.movement_step_params_compat import (
    test_max_translation_and_rotation_keys_are_parsed_converted_and_applied,
)

# ITP nonbonded -> forcefield/LJ energy wiring (2 functions)
from cpu_cli.itp_nonbonded_forcefield import (
    test_itp_atomtypes_builds_lj_matrix_and_deltaU_matches_analytic_lj,
    test_itp_nonbond_params_override_applies_nbfix_in_lj_matrix,
    test_itp_pairtypes_override_applies_nbfix_in_lj_matrix,
    test_itp_pairtypes_strict_mode_does_not_override_inter_residue,
    test_itp_nonbond_params_take_precedence_over_pairtypes_for_same_pair,
    test_dump_params_reports_itp_defaults_comb_rule,
    test_itp_defaults_gen_pairs_and_fudge_are_visible_and_strict_fails,
    test_itp_defaults_comb_rule_geometric_sigma_mixing,
    test_itp_defaults_comb_rule_one_converts_c6_c12_and_overrides,
)

# ITP nonbonded insertion->deletion closure regression (1 function)
from cpu_cli.itp_nonbonded_energy_closure import (
    test_itp_atomtypes_insertion_then_deletion_deltaU_matches_analytic_lj_pair_energy,
)

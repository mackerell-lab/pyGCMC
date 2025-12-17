"""
INP compatibility tests for gcmc_gpu-style units and legacy keys.

Kept as a small compatibility shim for `tests/platform/test_cpu_cli.py`.
The actual test functions live in smaller modules under `tests/platform/cpu_cli/`.
"""

from __future__ import annotations

from .inp_units_compat_accept_seed import (
    test_dump_accept_log_consistent_with_final_pdb_counts,
    test_inp_random_seed_used_when_cli_missing,
)
from .inp_units_compat_box import (
    test_inp_units_gcmc_gpu_box_roundtrip_cryst1,
    test_inp_units_gcmc_gpu_fragmuex_kcal_to_kj_log,
)
from .inp_units_compat_dump_region import (
    test_inp_units_auto_defaults_to_gcmc_gpu_angstrom_for_inp_files,
    test_inp_units_gcmc_gpu_converts_grid_dx_and_cutoffs_and_target_volume,
    test_inp_units_gcmc_gpu_gcmc_region_numeric_conversion_affects_volume,
)
from .inp_units_compat_fragmuex import (
    test_fragmuex_scales_activity_and_acceptance_in_ideal_gas_limit,
)
from .inp_units_compat_nbar import (
    test_nbar_const_water_nbar_scales_all_fragments_by_fragconc,
    test_nbar_const_water_nbar_scales_multiple_fragments_in_one_run,
    test_nbar_number_water_nbar_updates_activity_after_first_insertion,
    test_nbar_volume_based_mode_sets_activity_from_concentration,
)
from .inp_units_compat_outputs import test_gcmc_gpu_active_and_muex_files_are_written

__all__ = [
    "test_inp_units_gcmc_gpu_box_roundtrip_cryst1",
    "test_inp_units_gcmc_gpu_fragmuex_kcal_to_kj_log",
    "test_fragmuex_scales_activity_and_acceptance_in_ideal_gas_limit",
    "test_dump_accept_log_consistent_with_final_pdb_counts",
    "test_inp_random_seed_used_when_cli_missing",
    "test_inp_units_auto_defaults_to_gcmc_gpu_angstrom_for_inp_files",
    "test_inp_units_gcmc_gpu_converts_grid_dx_and_cutoffs_and_target_volume",
    "test_inp_units_gcmc_gpu_gcmc_region_numeric_conversion_affects_volume",
    "test_nbar_volume_based_mode_sets_activity_from_concentration",
    "test_nbar_const_water_nbar_scales_all_fragments_by_fragconc",
    "test_nbar_const_water_nbar_scales_multiple_fragments_in_one_run",
    "test_nbar_number_water_nbar_updates_activity_after_first_insertion",
    "test_gcmc_gpu_active_and_muex_files_are_written",
]

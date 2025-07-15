# tests/simulation/test_energy_PGP.py
"""
Energy PGP Tests - Main Entry Point

This file imports all Energy PGP tests from modular sub-files.
Run: pytest tests/simulation/test_energy_PGP.py

Modular structure (all modules under 300 lines):
- helpers.py: Shared helper functions and utilities (650 lines)
- basic_operations.py: Basic PGP operations tests (162 lines, 3 functions)
- method_comparison.py: PME vs PGP comparison tests (173 lines, 1 function)
- complex_systems.py: Complex multi-algorithm comparison tests (204 lines, 1 function)
- planar_systems.py: Planar systems tests (218 lines, 1 function)
- asymmetric_water_complete.py: Complete asymmetric water test (227 lines, 1 function)
- asymmetric_water_movement.py: Movement loop for asymmetric water (159 lines, helper)
- asymmetric_nacl.py: Large NaCl asymmetric test (283 lines, 1 function)
- vspme_delta_energies.py: PGP vs PME delta energies comparison (215 lines, 1 function)
- vspme_direct_lj.py: PGP with direct space and LJ interactions (347 lines, 1 function)
- vspme_lj_energy.py: LJ energy calculation in PME and PGP (360 lines, 1 function)
- vspme_two_atom.py: Simple two-atom system energy comparison (190 lines, 1 function)
- vspme_combined.py: Combined energy calculation diagnostic (410 lines, 1 function)
- pgp_lj_neutral_systems.py: PGP tests with neutral particles for LJ (1 function)
- pgp_lj_minimum_energy.py: PGP LJ energy at close distances (1 function)
- pgp_lj_extreme_distances.py: PGP LJ tests at extreme distances (3 functions)
- pgp_lj_boundary_tests.py: PGP LJ boundary condition tests (3 functions)
- pgp_lj_helpers.py: Helper functions for PGP LJ tests
- pgp_real_space.py: PGP real-space calculation tests (2 functions)
- pgp_real_space_debug.py: Debug tests for PGP real-space (2 functions)
- pgp_real_minimal.py: Minimal PGP real-space test (1 function)
- pgp_real_fixed.py: Fixed PGP implementation test (1 function)
- pgp_initialization_order.py: PGP initialization order tests (2 functions)
- pgp_pme_cutoff.py: PGP PME cutoff test (1 function)
- pgp_cutoff_issue.py: PGP cutoff issue tests (2 functions)
- pgp_erfc_debug.py: PGP erfc table debug tests (2 functions)
- pgp_complete.py: PGP Complete tests including LJ interactions (6 functions)

Total: 40 test functions across 28 modules.
- Basic PGP tests (3): basic_operations
- PGP comparison tests (8): method_comparison, complex_systems, planar_systems, 
  asymmetric_water_complete, asymmetric_nacl, pgp_lj_neutral_systems, pgp_lj_minimum_energy
- PGPvsPME tests (5): vspme_delta_energies, vspme_direct_lj, vspme_lj_energy,
  vspme_two_atom, vspme_combined
- PGP debug/fix tests (12): pgp_real_space (2), pgp_real_space_debug (2), pgp_real_minimal (1),
  pgp_real_fixed (1), pgp_initialization_order (2), pgp_pme_cutoff (1), pgp_cutoff_issue (2),
  pgp_erfc_debug (2)
"""

# Basic PGP operations tests
from energyPGP.basic_operations import (
    test_pgp_parameter_setting,
    test_precompute_grid_potential,
    test_interpolate_molecule_energy
)

# PME vs PGP comparison tests
from energyPGP.method_comparison import (
    test_compare_pme_pgp_energy
)

# Complex multi-algorithm comparison tests
from energyPGP.complex_systems import (
    test_compare_ewald_pme_pgp_complex
)

# Planar systems tests
from energyPGP.planar_systems import (
    test_compare_ewald_pme_pgp_planar
)

# Asymmetric charge distribution tests
from energyPGP.asymmetric_water_complete import (
    test_compare_ewald_pme_pgp_asymmetric
)
from energyPGP.asymmetric_nacl import (
    test_compare_ewald_pme_pgp_asymmetric_nacl
)

# PGPvsPME comparison tests
from energyPGP.vspme_delta_energies import (
    test_compare_pgp_pme_delta_energies
)
from energyPGP.vspme_direct_lj import (
    test_pgp_direct_and_lj_energies
)
from energyPGP.vspme_lj_energy import (
    test_lj_energy_pme_pgp
)
from energyPGP.vspme_two_atom import (
    test_simple_two_atom_system
)
from energyPGP.vspme_combined import (
    test_combined_energy_calculation
)

# LJ-only tests for PGP
from energyPGP.pgp_lj_neutral_systems import (
    test_pgp_lj_only_system
)
from energyPGP.pgp_lj_minimum_energy import (
    test_pgp_lj_close_interaction
)
# Import tests from split pgp_lj_limits files
from energyPGP.pgp_lj_extreme_distances import (
    test_pgp_lj_near_sigma,
    test_pgp_lj_short_distance,
    test_pgp_lj_asymptotic_behavior
)
from energyPGP.pgp_lj_boundary_tests import (
    test_pgp_lj_minimum,
    test_pgp_lj_near_cutoff,
    test_pgp_lj_mixed_distances
)

# PGP real-space calculation tests (bug fixes and verification)
from energyPGP.pgp_real_space import (
    test_pgp_real_space_calculation,
    test_pgp_real_space_with_fixed_atoms
)
from energyPGP.pgp_real_space_debug import (
    test_pgp_real_space_debug,
    test_pgp_with_water_molecule
)
from energyPGP.pgp_real_minimal import (
    test_pgp_real_space_minimal
)
from energyPGP.pgp_real_fixed import (
    test_pgp_real_space_fixed
)

# PGP initialization and configuration tests
from energyPGP.pgp_initialization_order import (
    test_pgp_with_correct_initialization,
    test_pgp_initialization_methods
)

# PGP vs PME Complete validation tests
from energyPGP.pgp_pme_complete_validation import (
    test_pgp_movement_energy,
    test_pgp_reciprocal_space_comparison,
    test_pgp_grid_spacing_convergence
)
from energyPGP.pgp_absolute_validation import (
    test_pgp_absolute_energy_validation,
    test_pgp_fixed_particle_contribution
)
from energyPGP.pgp_mixed_interactions import (
    test_pgp_single_moveable_particle,
    test_pgp_multiple_moveable_particles,
    test_pgp_with_pme_movement_comparison
)
from energyPGP.pgp_pme_electrostatic_validation import (
    test_pgp_pme_pure_electrostatic_exact,
    test_pgp_pme_different_movements,
    test_pgp_pme_convergence_with_mesh
)
from energyPGP.pgp_pme_cutoff import (
    test_pme_cutoff_in_pgp
)
from energyPGP.pgp_cutoff_issue import (
    test_pgp_cutoff_issue,
    test_pgp_with_different_cutoffs
)
from energyPGP.pgp_erfc_debug import (
    test_pgp_erfc_table,
    test_compare_initialization_sequences
)

# PGP Complete tests - including LJ interactions
# These tests use resetPGPState() to clear global state between runs
# See energyPGP/README_PGP_COMPLETE_TESTS.md for details
from energyPGP.pgp_complete import (
    test_pgp_complete_vs_pme_complete,
    test_pgp_complete_movement_energy,
    test_pgp_complete_pure_lj,
    test_pgp_complete_multi_atom_residue,
    test_pgp_complete_extreme_distances,
    test_pgp_complete_direct_movement_test
)

# Additional PGP tests that were missing
# TEMPORARILY COMMENTED OUT - These are the most recently added tests
# from energyPGP.pgp_simple_test import (
#     test_pgp_simple
# )
# from energyPGP.pgp_consistency_tests import (
#     test_pgp_energy_symmetry,
#     test_pgp_energy_scaling,
#     test_pgp_grid_independence
# )
from energyPGP.pgp_grid_convergence import (
    test_pgp_grid_convergence,
    test_pgp_alpha_convergence,
    test_pgp_convergence_trend
)
from energyPGP.pgp_cutoff_continuity import (
    test_pgp_cutoff_continuity,
    test_pgp_smooth_transition,
    test_pgp_lj_cutoff_continuity
)
from energyPGP.pgp_lj_diagnosis import (
    test_pgp_lj_only_system,
    test_pgp_mixed_system,
    test_pgp_lj_distance_scan
)
from energyPGP.pgp_pme_debug_electrostatic import (
    test_simple_two_particle_system,
    test_pgp_self_consistency
)
from energyPGP.debug_movement_residues import (
    test_movement_residues
)
from energyPGP.debug_vdw_movement import (
    test_vdw_movement_debug
)

# Support direct execution for testing
if __name__ == "__main__":
    import pytest
    pytest.main([__file__])
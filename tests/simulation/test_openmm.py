# tests/simulation/test_openmm.py
"""
OpenMM Tests - Main Entry Point

Imports all OpenMM comparison tests from energyOpenmm/ subdirectory.
Total: 108 test functions including Drude tests.
"""

# Suppress SWIG-related deprecation warnings
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning, 
                       message="builtin type SwigPyPacked has no __module__ attribute")
warnings.filterwarnings("ignore", category=DeprecationWarning, 
                       message="builtin type SwigPyObject has no __module__ attribute")
warnings.filterwarnings("ignore", category=DeprecationWarning, 
                       message="builtin type swigvarlink has no __module__ attribute")

# Naive nonbonded comparison tests
from energyOpenmm.naive_energy_components import (
    test_openmm_energy_components,
    test_naive_energy_components,
    test_compare_openmm_naive_nonbonded
)

# Naive PBC and cutoff tests
from energyOpenmm.naive_pbc_cutoff import (
    test_compare_pbc_energies,
    test_compare_cutoff_effects
)

# Simple interaction tests
from energyOpenmm.simple_attractive_repulsive import (
    test_attractive_interaction,
    test_repulsive_interaction
)

# Simple PBC and cutoff tests
from energyOpenmm.simple_pbc_cutoff import (
    test_pbc_interaction,
    test_cutoff_effect
)

# Intermediate energy tests
from energyOpenmm.intermediate_energy_symmetry import test_energy_symmetry
from energyOpenmm.intermediate_custom_vs_standard import test_compare_custom_vs_standard_nonbonded

# Method comparison tests
from energyOpenmm.method_combined import test_compare_nonbonded_methods
from energyOpenmm.method_nocutoff import test_compare_nonbonded_nocutoff
from energyOpenmm.method_cutoff_nonperiodic import test_compare_nonbonded_cutoff_nonperiodic
from energyOpenmm.method_switching_functions import test_compare_switching_functions

# Energy analysis tests
from energyOpenmm.analysis_separate_terms import test_compare_separate_terms
from energyOpenmm.analysis_naive_vs_cutoff import test_compare_naive_vs_cutoff_energy
from energyOpenmm.analysis_detailed_comparison import test_detailed_energy_comparison_simple_vs_openmm_cutoff

# Periodic boundary tests
from energyOpenmm.periodic_self_energy import test_compare_all_methods_with_self_energy
from energyOpenmm.periodic_energy_terms import test_analyze_openmm_energy_terms
from energyOpenmm.periodic_cutoff_comparison import test_cutoff_periodic_comparison
from energyOpenmm.periodic_lj_coulomb import test_separate_lj_coulomb_periodic
from energyOpenmm.periodic_force_parameters import test_compare_force_parameters_and_energies

# File-based test
from energyOpenmm.file_based import test_verify_openmm_expressions

# NBFIX tests
from energyOpenmm.nbfix_parameters import (
    test_nbfix_overrides_lj_combination,
    test_multiple_nbfix_pairs
)
from energyOpenmm.nbfix_energy import test_ion_water_nbfix_energy
from energyOpenmm.nbfix_openmm import (
    test_nbfix_energy_vs_openmm,
    test_nbfix_energy_components_separately
)
from energyOpenmm.nbfix_pbc import test_nbfix_across_pbc

# PME Total bug fix tests
from energyOpenmm.pme_total_fix import (
    test_pme_total_single_residue,
    test_pme_total_multiple_residues,
    test_pme_total_vs_ewald
)

# PGP strict tolerance test
from energyOpenmm.pgp_strict_tolerance import test_pgp_with_strict_tolerance

# PME-Ewald consistency test
from energyOpenmm.pme_ewald_consistency import test_pme_ewald_consistency

# PME regression tests
from energyOpenmm.pme_total_lj_hypothesis import test_lj_hypothesis
from energyOpenmm.pme_residue_offset_pattern import test_residue_offset_pattern
from energyOpenmm.verify_pme_total_source import test_pme_total_sources
from energyOpenmm.analyze_pme_total_bug import test_pme_total_configurations

# PME ligand energy comparison
from energyOpenmm.pme_ligand_openmm_comparison import (
    test_pme_ligand_energy_vs_openmm,
    test_pme_movement_residues
)

# High precision PME comparison
from energyOpenmm.pme_high_precision_comparison import (
    test_pme_high_precision_simple_system,
    test_pme_high_precision_complex_system,
    test_pme_convergence_with_mesh_size
)

# PME ligand insertion tests
from energyOpenmm.pme_ligand_insertion import (
    test_ligand_insertion_energy,
    test_ligand_position_scan
)

# Movement residue PME tests
from energyOpenmm.movement_residue_pme import (
    test_movement_residue_pme,
    test_movement_residue_removal
)

# Residue activation demonstrations (not tests)

# Full nonbonded comparison tests
from energyOpenmm.full_nonbonded_pme import test_full_nonbonded_pme_lj
from energyOpenmm.full_nonbonded_cutoff import test_lj_only_comparison

# LJ double counting analysis tests
from energyOpenmm.lj_double_counting import (
    test_single_pair_lj,
    test_three_atom_system,
    test_residue_assignment
)

# LJ detailed comparison tests
from energyOpenmm.lj_detailed_comparison import (
    test_simple_two_atom_system,
    test_mixed_types
)

# PME LJ-only tests
from energyOpenmm.pme_lj_only_basic import (
    test_lj_only_cutoff,
    test_lj_only_pme
)
from energyOpenmm.pme_lj_only_analysis import (
    test_lj_distance_scan,
    test_lj_mixing_rules
)

# Debug distance tests
from energyOpenmm.debug_distances import (
    test_system_from_test_full_nonbonded,
    test_system_from_test_pme_lj_only
)

# LJ close distance tests
from energyOpenmm.lj_close_distance import (
    test_close_distances,
    test_energy_capping
)

# LJ double counting summary test
from energyOpenmm.lj_double_counting_summary import test_lj_double_counting_proof

# LJ multiparticle analysis tests
from energyOpenmm.lj_multiparticle_analysis import (
    test_progressive_atoms,
    test_cutoff_effects,
    test_pbc_boundary
)

# Fixed functions tests
from energyOpenmm.fixed_functions import (
    test_cutoff_fixed,
    test_pme_fixed,
    test_fixed_vs_openmm
)

# PME convergence and comparison tests
from energyOpenmm.pme_convergence import test_pme_parameter_sensitivity
from energyOpenmm.pme_medium_complexity import test_pme_medium_complexity
from energyOpenmm.pme_simple_comparison import test_simple_pme_comparison

# Intra-residue and NaCl tests
from energyOpenmm.intra_residue_difference import test_intra_residue_handling
from energyOpenmm.pme_nacl_correct import test_pme_vs_ewald_correct

# Intramolecular LJ tests
from energyOpenmm.intramolecular_lj import (
    test_single_molecule_lj,
    test_two_molecules_lj
)

# PME Complete tests
from energyOpenmm.pme_complete import (
    test_complete_cutoff_lj_only,
    test_complete_pme_lj_only,
    test_complete_pme_full_system,
    test_complete_vs_fixed_difference,
    test_complete_consistency
)

# PME with OpenMM parameters
from energyOpenmm.pme_with_openmm_params import (
    test_pme_with_openmm_params,
    test_simple_two_charge_system
)

# Drude-OpenMM comparison tests
from energyOpenmm.openmm_comparison_tests import (
    test_single_water_energy,
    # test_water_dimer_interaction,  # Skip: causes issues
    test_scf_convergence_tolerance,
    test_polarization_response,
    test_multiple_water_box,
    test_openmm_spring_constant_consistency,
    test_openmm_thole_screening_values,
    test_openmm_water_dipole_moment,
    test_openmm_scf_convergence_behavior
)

# Drude numerical comparison tests
from energyOpenmm.drude_numerical_comparison import (
    test_single_water_vacuum_exact,
    test_water_external_field_exact,
    test_water_dimer_exact,
    test_spring_constant_exact_values,
    test_thole_screening_exact
)

# Drude anisotropic polarizability tests (not yet implemented)
# from energyOpenmm.drude_anisotropic_tests import (
#     test_anisotropic_polarizability_setup,
#     test_anisotropic_response_to_field,
#     test_anisotropic_thole_screening,
#     test_anisotropic_energy_calculation,
#     test_anisotropic_water_model
# )

# Drude force verification tests
from energyOpenmm.drude_force_verification import (
    test_harmonic_spring_force,
    test_coulomb_force_consistency,
    test_drude_scf_force_balance,
    test_thole_screened_force,
    test_force_energy_consistency_complex
)

# Drude extreme condition tests
from energyOpenmm.drude_extreme_conditions import (
    test_scf_convergence_failure,
    test_hardwall_constraint_extreme,
    test_zero_polarizability,
    test_many_body_polarization,
    test_pathological_geometry
)

# Drude SCF comparison tests
from energyOpenmm.drude_scf_comparison import (
    test_single_drude_with_field
    # skip_test_water_dimer_scf  # Skip: segfault
)

if __name__ == "__main__":
    import pytest
    pytest.main([__file__])
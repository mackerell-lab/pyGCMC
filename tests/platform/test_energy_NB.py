# tests/simulation/test_energy_NB.py
"""
Energy Non-Bonded Tests - Main Entry Point

This file imports all Energy Non-Bonded tests from modular sub-files.
Run: pytest tests/simulation/test_energy_NB.py

Modular structure:
- Original test_naive_nonbonded.py (1278 lines, 17 functions) split into nonbonded_*.py modules
- Original test_switching_function.py (744 lines, 9 functions) split into switching_*.py modules

Non-bonded tests (17 functions from test_naive_nonbonded.py):
- Basic interaction tests: attractive, repulsive, electrostatic, combined
- Movement molecule tests: three movement molecules
- Validation tests: invalid parameters, inactive residues, energy symmetry
- Distance handling: very close distance, zero distance
- System-wide tests: all residues nonbonded/inactive
- Cutoff tests: cutoff nonperiodic, PBC basic/invalid/comparison/cross boundary

Switching function tests (9 functions from test_switching_function.py):
- Core switching tests: switching function calculation, energy with switching
- Method comparison: with Ewald, Monte Carlo system integration
- Diagnostic tests: print values, energy comparison, internal function validation
- Energy calculation tests: MCS energy calculation with switching

Energy components tests (3 functions):
- getTotalEnergyComponents function tests for simple systems, VdW interactions, and multiple residues

Total: 29 test functions across modular files (all modules under 300 lines each).
"""

# Non-bonded basic interaction tests
from energyNB.nonbonded_attractive import (
    test_attractive_interaction,
    test_electrostatic_interaction
)

from energyNB.nonbonded_repulsive import (
    test_repulsive_interaction,
    test_combined_interaction
)

# Non-bonded movement and validation tests
from energyNB.nonbonded_movement import (
    test_three_movement_molecules,
    test_invalid_forcefield_params,
    test_inactive_residue
)

# Non-bonded energy and distance tests
from energyNB.nonbonded_energy import (
    test_energy_symmetry,
    test_very_close_distance,
    test_zero_distance_handling
)

# Analytic energy contracts (LJ/Coulomb)
from energyNB.analytic_energy_contracts import (
    test_system_energy_lj_only_analytic,
    test_system_energy_coulomb_only_analytic,
    test_system_energy_components_mixed_analytic
)

# Non-bonded system-wide tests
from energyNB.nonbonded_system import (
    test_all_residues_nonbonded,
    test_all_residues_inactive
)

# Non-bonded cutoff tests
from energyNB.nonbonded_cutoff import (
    test_cutoff_nonperiodic
)

# Non-bonded PBC tests
from energyNB.nonbonded_pbc_basic import (
    test_pbc_basic,
    test_pbc_invalid_box,
    test_pbc_vs_nopbc
)

from energyNB.nonbonded_pbc_advanced import (
    test_pbc_cross_boundary
)

# Cutoff/PBC contracts
from energyNB.cutoff_pbc_contracts import (
    test_cutoff_excludes_far_pair,
    test_pbc_minimum_image_distance
)

# Pairtypes strict (1-4) behavior for DIRECT cutoff
from energyNB.pairtypes_14_strict_direct import (
    test_pairtypes_14_override_applies_in_direct_cutoff
)

# Switching function core tests
from energyNB.switching_core import (
    test_switching_function,
    test_energy_with_switching
)

# Switching function method tests
from energyNB.switching_methods import (
    test_with_ewald,
    test_monte_carlo_system_switching_function
)

# Switching function diagnostic tests
from energyNB.switching_print import (
    test_print_switching_values
)

from energyNB.switching_comparison import (
    test_compare_energy_with_without_switching
)

from energyNB.switching_internal import (
    test_internal_switching_function,
    test_internal_switching_function_ewald,
    test_mcs_energy_calculation_with_switching
)

# Energy components tests (getTotalEnergyComponents function)
from energyNB.energy_components_tests import (
    test_get_total_energy_components_simple,
    test_get_total_energy_components_with_vdw,
    test_get_total_energy_components_multiple_residues
)

# Support direct execution for testing
if __name__ == "__main__":
    import pytest
    pytest.main([__file__])

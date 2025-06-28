# tests/system/test_MonteCarloSystem.py
"""
MonteCarloSystem Tests - Main Entry Point

This file imports all MonteCarloSystem tests from modular sub-files.
Run: pytest tests/system/test_MonteCarloSystem.py

Modular structure (all files under 250 lines):
- Original test_MonteCarloSystem.py (1443 lines, 14 functions) split into:
  * MonteCarloSystem/initialization_tests.py: Basic initialization tests (3 functions)
  * MonteCarloSystem/mapping_tests.py: Type mapping and property tests (3 functions)  
  * MonteCarloSystem/movement_single.py: Single movement molecule test (1 function)
  * MonteCarloSystem/movement_two.py: Two movement molecules test (1 function)
  * MonteCarloSystem/movement_three.py: Three movement molecules test (1 function)
  * MonteCarloSystem/type_mapping_after_movement.py: Type mapping after movement test (1 function)
  * MonteCarloSystem/molecule_addition_comparison.py: Addition comparison test (1 function)
  * MonteCarloSystem/forcefield_tests.py: Force field initialization test (1 function)
- Helper fixtures and common code: MonteCarloSystem/helpers.py

Total: 14 test functions across modular files (all modules under 250 lines each).
Original total: 1443 lines → New structure: 1481 lines (+3% for better organization)
"""

# Import fixtures from helpers
from MonteCarloSystem.helpers import molecular_system, assert_arrays_almost_equal, charmm_ff

# Basic initialization tests (3 functions)
from MonteCarloSystem.initialization_tests import (
    test_initialize_from_molecular,
    test_initialize_from_molecular_empty,
    test_initialize_from_molecular_large_system,
    test_type_mapping
)

# Type mapping and property tests (3 functions)
from MonteCarloSystem.mapping_tests import (
    test_residue_atom_properties,
    test_residue_atom_properties_empty,
    test_type_mapping_from_molecular,
    test_residue_type_mapping_from_molecular
)

# Movement molecule tests (3 functions)
from MonteCarloSystem.movement_single import test_add_movement_molecules
from MonteCarloSystem.movement_two import test_add_two_movement_molecules  
from MonteCarloSystem.movement_three import test_add_three_movement_molecules

# Comparison and mapping tests (2 functions)
from MonteCarloSystem.type_mapping_after_movement import test_atom_type_maps_after_movement_molecules
from MonteCarloSystem.molecule_addition_comparison import test_compare_add_molecules_together_vs_separate

# Force field tests (1 function)
from MonteCarloSystem.forcefield_tests import test_initialize_force_field

# Support direct execution for testing
if __name__ == "__main__":
    import pytest
    pytest.main([__file__])
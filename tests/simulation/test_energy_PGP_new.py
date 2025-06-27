# tests/simulation/test_energy_PGP_new.py
"""
Energy PGP Tests - Main Entry Point

This file imports all Energy PGP tests from modular sub-files.
Run: pytest tests/simulation/test_energy_PGP_new.py

Modular structure (all modules under 300 lines):
- helpers.py: Shared helper functions and utilities (267 lines)
- basic_operations.py: Basic PGP operations tests (162 lines, 3 functions)
- method_comparison.py: PME vs PGP comparison tests (173 lines, 1 function)
- complex_systems.py: Complex multi-algorithm comparison tests (204 lines, 1 function)
- planar_systems.py: Planar systems tests (218 lines, 1 function)
- asymmetric_water_complete.py: Complete asymmetric water test (227 lines, 1 function)
- asymmetric_water_movement.py: Movement loop for asymmetric water (159 lines, helper)
- asymmetric_nacl.py: Large NaCl asymmetric test (283 lines, 1 function)

Total: 8 test functions across 8 modules, all under 300 lines each.
"""

# Basic PGP operations tests
from pgp.basic_operations import (
    test_pgp_parameter_setting,
    test_precompute_grid_potential,
    test_interpolate_molecule_energy
)

# PME vs PGP comparison tests
from pgp.method_comparison import (
    test_compare_pme_pgp_energy
)

# Complex multi-algorithm comparison tests
from pgp.complex_systems import (
    test_compare_ewald_pme_pgp_complex
)

# Planar systems tests
from pgp.planar_systems import (
    test_compare_ewald_pme_pgp_planar
)

# Asymmetric charge distribution tests
from pgp.asymmetric_water_complete import (
    test_compare_ewald_pme_pgp_asymmetric
)
from pgp.asymmetric_nacl import (
    test_compare_ewald_pme_pgp_asymmetric_nacl
)

# Support direct execution for testing
if __name__ == "__main__":
    import pytest
    pytest.main([__file__])
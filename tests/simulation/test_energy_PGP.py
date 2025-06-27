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

Total: 13 test functions across 13 modules (8 PGP + 5 PGPvsPME).
- PGP tests (8): basic_operations (3), method_comparison (1), complex_systems (1), 
  planar_systems (1), asymmetric_water_complete (1), asymmetric_nacl (1)
- PGPvsPME tests (5): vspme_delta_energies (1), vspme_direct_lj (1), vspme_lj_energy (1),
  vspme_two_atom (1), vspme_combined (1)
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

# Support direct execution for testing
if __name__ == "__main__":
    import pytest
    pytest.main([__file__])
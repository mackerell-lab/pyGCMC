# tests/simulation/test_energy_PME_new.py
"""
Energy PME Tests - Main Entry Point

This file imports all Energy PME tests from modular sub-files.
Run: pytest tests/simulation/test_energy_PME_new.py

Consolidated modular structure (all modules under 300 lines):
- helpers.py: Shared NaCl crystal creation helper (87 lines)
- crystal_file_helpers.py: File reading helper functions (133 lines)
- basic_tests.py: Basic PME tests - initialization, PME vs Ewald, movement (171 lines, 3 functions)
- parameter_optimization.py: Parameter optimization tests - spline order, error tolerance, mesh accuracy (225 lines, 3 functions)
- advanced_parameters.py: Advanced parameter tests - alpha dependency, LJ energy (185 lines, 2 functions)
- comparison_tests.py: PME comparison tests - Ewald vs PME with different systems (293 lines, 2 functions)
- small_mesh.py: PME small mesh accuracy test (130 lines, 1 function)
- ewald_exact.py: Exact Ewald test replicating pme.cpp (268 lines, 1 function)
- grid_operations.py: PME grid operations diagnostic test (192 lines, 1 function)
- pme_parameters.py: PME parameters influence test (195 lines, 1 function)
- cutoff_dependence.py: PME cutoff dependence test (177 lines, 1 function)

Total: 15 test functions across 11 modules (reduced from 17 modules).
Original file: 2050 lines → New consolidated structure: 1,856 lines total (-9% reduction)
All functions maintain complete compatibility with original implementation.

Consolidated grouping:
- Basic Tests (3): Initialization, PME vs Ewald comparison, Movement energy
- Parameter Optimization (3): Spline order effects, Error tolerance, Mesh accuracy
- Advanced Parameters (2): Alpha dependency, LJ energy calculation
- Comparison Tests (2): Ewald vs PME with random and amorphous systems
- Specialized Tests (5): Small mesh, Exact Ewald, Grid operations, PME parameters, Cutoff dependence
"""

# Basic PME tests (3 functions)
from energyPME.basic_tests import (
    test_pme_initialization,
    test_pme_vs_ewald,
    test_pme_movement_energy
)

# Parameter optimization tests (3 functions)
from energyPME.parameter_optimization import (
    test_pme_spline_order,
    test_pme_error_tolerance,
    test_pme_mesh_accuracy
)

# Advanced parameter tests (2 functions)
from energyPME.advanced_parameters import (
    test_pme_alpha_dependency,
    test_pme_lj_energy
)

# PME comparison tests (2 functions)
from energyPME.comparison_tests import (
    test_ewald_vs_pme_comparison,
    test_ewald_vs_pme_random
)

# PME small mesh accuracy test
from energyPME.small_mesh import (
    test_pme_small_mesh
)

# Exact Ewald test replicating pme.cpp
from energyPME.ewald_exact import (
    test_ewald_exact
)

# PME grid operations diagnostic test
from energyPME.grid_operations import (
    test_pme_grid_operations
)

# PME parameters influence test
from energyPME.pme_parameters import (
    test_pme_parameters
)

# PME cutoff dependence test
from energyPME.cutoff_dependence import (
    test_cutoff_dependence
)

# Support direct execution for testing
if __name__ == "__main__":
    import pytest
    pytest.main([__file__])
# tests/simulation/test_energy_drude.py
"""
Energy Drude Tests - Main Entry Point

This file imports all Energy Drude tests from modular sub-files.
Run: pytest tests/simulation/test_energy_drude.py

Test structure:
- helpers.py: Helper functions for creating Drude systems
- basic_tests.py: Basic Drude oscillator tests
- scf_tests.py: SCF convergence and optimization tests
- water_tests.py: SWM4-NDP water model tests
- screened_tests.py: Thole screening interaction tests
"""

# Basic Drude tests
from energyDrude.basic_tests import (
    test_drude_initialization,
    test_drude_harmonic_energy,
    test_drude_equilibrium_position,
    test_drude_anisotropic_polarizability,
    test_drude_multiple_particles
)

# SCF convergence tests
from energyDrude.scf_tests import (
    test_drude_scf_convergence,
    test_drude_scf_tolerance,
    test_drude_scf_damping,
    test_drude_scf_multiple_particles,
    test_drude_scf_with_external_field,
    test_drude_scf_iteration_limit
)

# Water model tests
from energyDrude.water_tests import (
    test_swm4_water_single,
    test_swm4_water_dimer,
    test_swm4_water_box,
    test_swm4_water_polarization
)

# Screened interaction tests
from energyDrude.screened_tests import (
    test_thole_screening,
    test_screened_pair_energy,
    test_multiple_screened_pairs,
    test_thole_parameter_effect
)

# Support direct execution for testing
if __name__ == "__main__":
    import pytest
    pytest.main([__file__])
#!/usr/bin/env python3
"""
Main test file for Drude oscillator functionality.
Imports tests from organized submodules.
"""

import pytest

# Import all tests from submodules
from energyDrude.basic_tests import (
    test_drude_imports,
    test_drude_particle_creation,
    test_screened_pair_creation,
    test_drude_scf_params,
    test_drude_algorithm_enum,
    test_thole_screening_function,
    test_drude_complete_basic
)

from energyDrude.energy_tests import (
    test_drude_with_mcstate,
    test_harmonic_energy,
    test_zero_displacement_energy
)

from energyDrude.water_tests import (
    test_swm4_ndp_water_setup,
    test_two_water_thole_interaction,
    test_water_residue_setup
)

if __name__ == "__main__":
    pytest.main([__file__, "-v"])
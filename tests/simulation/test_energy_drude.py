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
    test_drude_complete_basic,
    test_simple_drude_force
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

from energyDrude.thole_tests import (
    test_thole_screening_function as test_thole_screening_detailed,
    test_thole_screening_limits,
    test_thole_dipole_dipole_interaction,
    test_thole_parameter_sensitivity,
    test_multiple_thole_pairs
)

from energyDrude.hardwall_basic_tests import (
    test_hardwall_constraint_basic,
    test_hardwall_different_limits,
    test_hardwall_energy_consistency
)

from energyDrude.hardwall_advanced_tests import (
    test_hardwall_with_multiple_drudes,
    test_hardwall_anisotropic_force
)

# OpenMM comparison tests moved to test_openmm.py

from energyDrude.scf_basic_convergence_tests import (
    test_scf_basic_convergence,
    test_scf_damping_factor_effect,
    test_scf_iteration_limit
)

from energyDrude.scf_advanced_convergence_tests import (
    test_scf_tolerance_scaling,
    test_scf_adaptive_damping
)

from energyDrude.charmm_basic_validation_tests import (
    test_drude_charge_relationship,
    test_induced_dipole_in_uniform_field,
    test_thole_screening_values,
    test_swm4_ndp_geometry
)

from energyDrude.charmm_advanced_validation_tests import (
    test_drude_mass_redistribution,
    test_polarization_catastrophe_prevention,
    test_anisotropic_polarizability
)

from energyDrude.omm_scf_consistency_tests import test_scf_repeat_consistency
from energyDrude.omm_force_validation_tests import test_numerical_force_validation
from energyDrude.omm_thole_validation_tests import test_thole_screening_validation
from energyDrude.omm_swm4ndp_water_tests import test_swm4_ndp_water_system

if __name__ == "__main__":
    pytest.main([__file__, "-v"])
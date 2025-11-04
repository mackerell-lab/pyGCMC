# tests/simulation/test_movement_gcmc.py
"""
GCMC Movement Module Tests - Main Entry Point

This file imports all GCMC movement tests from modular sub-files.
Run: pytest tests/simulation/test_movement_gcmc.py

Test categories:
1. Basic GCMC operations - Insertion/deletion/translation/rotation
2. Acceptance criteria - μVT ensemble acceptance rates
3. Energy calculations - System and fragment energies
4. Cavity bias - Cavity detection and biased insertion
5. Fragment reservoir - Template and instance management
6. Scaling tests - Volume and chemical potential scaling
"""

import pytest
import pygcmc

# Import fixtures from movementGCMC/conftest.py for all tests in this file
from movementGCMC.conftest import (
    setup_gcmc_state, 
    gcmc_params,
    water_template,
    gcmc_mover,
    small_state,
    populated_state
)

# Basic GCMC operation tests (4 functions)
from movementGCMC.basic_operations import (
    test_insertion_attempt,
    test_deletion_attempt,
    test_translation_attempt,
    test_rotation_attempt
)

# Acceptance criteria tests (3 functions)
from movementGCMC.acceptance_tests import (
    test_acceptance_rate_range,
    test_chemical_potential_effect,
    test_temperature_effect
)

from movementGCMC.acceptance_formula_engine import (
    test_insertion_probability_detailed_matches_formula,
    test_deletion_probability_detailed_matches_formula
)

# Energy calculation tests (3 functions)
from movementGCMC.energy_tests import (
    test_energy_conservation,
    test_fragment_energy_calculation,
    test_system_energy_consistency
)

# Cavity bias tests (4 functions)
from movementGCMC.cavity_bias_tests import (
    test_cavity_detection,
    test_cavity_bias_insertion,
    test_cavity_grid_parameters,
    test_cavity_vs_uniform_insertion
)

# Fragment reservoir tests (4 functions)
from movementGCMC.reservoir_tests import (
    test_template_management,
    test_instance_creation_deletion,
    test_ghost_fragment_recycling,
    test_reservoir_statistics
)

# Scaling and ensemble tests (3 functions)
from movementGCMC.scaling_tests import (
    test_volume_scaling,
    test_chemical_potential_scaling,
    test_detailed_balance
)

# Configuration tests (3 functions)
from movementGCMC.configuration_tests import (
    test_statistics_sampling_occurs,
    test_max_translation_applied_to_moves,
    test_cavity_bias_toggle_bias_value
)

# =====================================================
# Additional tests migrated from platform/cpu/
# =====================================================

# Acceptance probability tests (functions)
from movementGCMC.acceptance_probability import (
    test_acceptance_probability_with_env,
    test_acceptance_probability_without_env,
    test_seed_chain
)

# Test classes from migrated files
from movementGCMC.acceptance_pure import TestAcceptancePure
from movementGCMC.detailed_balance import TestDetailedBalance
from movementGCMC.detailed_balance_fix import TestDetailedBalanceFix
from movementGCMC.insertion_deletion_balance import TestInsertionDeletionBalance
from movementGCMC.ideal_gas_validation import TestIdealGasValidation
from movementGCMC.gcmc_engine_config import TestGCMCEngineConfig
from movementGCMC.cavity_bias_cpu import TestCavityBias as TestCavityBiasCPU

# CLI/JSONL integration suites
from movementGCMC.cavity_bias_cli import TestCavityBias as TestCavityBiasCLI
from movementGCMC.cbmc_sanity_cli import TestCBMCSanity

# Temperature scaling tests (functions)
from movementGCMC.temperature_scaling import (
    test_boltzmann_distribution,
    test_metropolis_criterion,
    test_temperature_effect_on_acceptance,
    test_zero_temperature_limit
)

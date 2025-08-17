# tests/simulation/test_movement_cpp.py
"""
C++ Movement Module Bindings Tests - Main Entry Point

This file imports all C++ movement binding tests from modular sub-files.
Run: pytest tests/simulation/test_movement_cpp.py

Test categories:
1. Parameter tests - MovementParams configuration
2. Module tests - MovementModule functionality  
3. Active pool tests - Memory management
4. Result tests - MovementResult structure
5. Integration tests - Complete GCMC workflows
"""

# Parameter configuration tests (3 functions)
from movementCPP.params_tests import (
    test_default_params,
    test_custom_temperature,
    test_update_derived_parameters
)

# Main module functionality tests (7 functions)
from movementCPP.module_tests import (
    test_module_creation,
    test_insertion_attempt,
    test_deletion_attempt,
    test_translation_attempt,
    test_rotation_attempt,
    test_cavity_finding,
    test_statistics_tracking
)

# Active pool memory management tests (5 functions)
from movementCPP.active_pool_tests import (
    test_pool_creation,
    test_insert_molecule,
    test_delete_residue,
    test_fragmentation,
    test_capacity_check
)

# Movement result tests (2 functions)
from movementCPP.result_tests import (
    test_result_creation,
    test_result_summary
)

# Integration and workflow tests (2 functions)
from movementCPP.integration_tests import (
    test_gcmc_equilibration,
    test_acceptance_rates
)
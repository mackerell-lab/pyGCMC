# tests/simulation/test_movement_python.py
"""
Movement-based insertion tests - Main Entry Point

This file imports all movement insertion tests from modular sub-files.
Run: pytest tests/simulation/test_movement_python.py

Test categories:
1. Basic insertion tests - fundamental molecule insertion operations
2. Advanced insertion tests - cavity detection and energy-guided insertion
3. Residue activation tests - GCMC residue activation/deactivation
4. Insertion/deletion tests - GCMC insertion and deletion energy calculations
5. GCMC algorithm tests - Detailed balance and acceptance criteria
6. Practical GCMC tests - Real energy calculations with PyGCMC
"""

# Basic insertion tests (4 functions) - test_insert_with_existing_molecules moved to tmp/skipped_tests
from movementPython.basic_insertion_tests import (
    test_insert_single_molecule,
    test_insert_multiple_water_molecules,
    test_insert_ion_pair,
    test_insert_random_positions
)

# Advanced insertion tests (6 functions)
from movementPython.advanced_insertion_basic_tests import (
    test_simple_two_particle_insertion,
    test_cavity_detection
)
from movementPython.advanced_insertion_biased_tests import (
    test_biased_cavity_insertion,
    test_energy_guided_insertion
)
from movementPython.advanced_insertion_sequential_tests import (
    test_sequential_cavity_filling,
    test_water_cluster_formation
)

# Residue activation tests (4 functions) - moved from energyGCMC
from movementPython.residue_activation_tests import (
    test_residue_addition_energy,
    test_movement_energy_calculation,
    test_residue_activation_deactivation,
    test_multiple_residue_types
)

# GCMC insertion/deletion tests (6 functions) - moved from energyGCMC
from movementPython.insertion_deletion_tests import (
    test_single_molecule_insertion,
    test_molecule_insertion_with_interactions,
    test_molecule_deletion,
    test_multiple_insertions,
    test_biased_insertion,
    test_movement_only_insertion
)

# GCMC algorithm tests (3 functions)
from movementPython.gcmc_insertion_test_functions import (
    test_gcmc_benzene_insertion_in_protein,
    test_gcmc_water_insertion_simple,
    test_verify_detailed_balance
)

# Practical GCMC tests (3 functions)
from movementPython.practical_gcmc_benzene_test import test_practical_benzene_insertion
from movementPython.practical_gcmc_water_test import test_water_insertion_with_real_energy
from movementPython.practical_gcmc_cavity_bias_test import test_cavity_bias_effect

# Numerical stability tests (7 functions) - NEW
from movementPython.numerical_stability_tests import (
    test_extreme_energy_numerical_stability,
    test_configurational_bias_numerical_stability,
    test_temperature_extremes,
    test_gcmc_insertion_acceptance_stability,
    test_parallel_rosenbluth_calculation,
    test_insertion_with_pbc_wrapping
)

# Cavity bias tests (7 functions) - NEW (split into two files)
from movementPython.cavity_bias_basic_tests import (
    test_cavity_grid_creation,
    test_cavity_detection_with_molecules,
    test_biased_insertion_acceptance,
    test_cavity_bias_with_different_grid_spacings
)
from movementPython.cavity_bias_advanced_tests import (
    test_cavity_bias_energy_calculation,
    test_adaptive_cavity_grid,
    test_cavity_bias_detailed_balance
)

# Log-space stability tests (9 functions) - NEW (split into two files)
from movementPython.logspace_basic_tests import (
    test_logsumexp_extreme_values,
    test_prob_from_log_saturation,
    test_insertion_acceptance_extreme_parameters,
    test_deletion_acceptance_extreme_parameters,
    test_rosenbluth_weight_extreme_energies
)
from movementPython.logspace_advanced_tests import (
    test_configurational_bias_logspace,
    test_detailed_balance_logspace,
    test_extreme_temperature_logspace,
    test_parallel_logspace_calculations
)

# Integrated GCMC patterns tests (6 functions) - NEW (split into two files)
from movementPython.integrated_gcmc_basic import (
    test_config_driven_gcmc,
    test_gcmc_md_integration_pattern,
    test_multiple_fragment_types
)
from movementPython.integrated_gcmc_advanced import (
    test_convergence_monitoring,
    test_acceptance_rate_adaptation,
    test_realistic_waterbox_setup
)

# Translation and Rotation tests (5 functions) - NEW (split into two files)
from movementPython.translation_tests import (
    test_translation_move,
    test_rotation_move
)
from movementPython.rotation_tests import (
    test_combined_translation_rotation,
    test_translation_with_pbc,
    test_rotation_preserves_bond_lengths
)

# Pool memory management tests (8 functions) - NEW
from movementPython.pool_memory_management_tests import (
    test_active_pool_basic_operations,
    test_active_pool_energy_conservation,
    test_batched_operations,
    test_fragmentation_management,
    test_gpu_ready_capacity_limits,
    test_energy_component_separation,
    test_concurrent_insert_delete_pattern
)

# Pool simple test
def test_pool_simple():
    """Test simple pool operations"""
    from movementPython.active_pool import ActivePool
    import pygcmc
    
    # Simple water molecule 
    def create_water():
        atoms = []
        # O
        o = pygcmc.MCAtom()
        o.x, o.y, o.z = 0, 0, 0
        o.charge = -0.834
        o.type = 0
        atoms.append(o)
        # H1
        h1 = pygcmc.MCAtom()
        h1.x, h1.y, h1.z = 0.1, 0, 0
        h1.charge = 0.417
        h1.type = 1
        atoms.append(h1)
        # H2
        h2 = pygcmc.MCAtom()
        h2.x, h2.y, h2.z = 0, 0.1, 0
        h2.charge = 0.417
        h2.type = 1
        atoms.append(h2)
        return atoms
    
    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 2
    ff.ljSigma = [0.315, 0.0]
    ff.ljEps = [0.636, 0.0]
    
    pool = ActivePool(box=[5.0, 5.0, 5.0], cutoff=1.2, forcefield=ff)
    
    # Insert 3 waters
    res0 = pool.insert_molecule(create_water())
    assert len(pool.state.atoms) == 3
    assert len(pool.residue_metadata) == 1
    
    res1 = pool.insert_molecule(create_water())
    assert len(pool.state.atoms) == 6
    assert len(pool.residue_metadata) == 2
    
    res2 = pool.insert_molecule(create_water())
    assert len(pool.state.atoms) == 9
    assert len(pool.residue_metadata) == 3
    
    # Delete res1
    pool.delete_residue(res1)
    assert pool.residue_metadata[res1].active == False
    
    # Compact
    compacted = pool.compact(force=True)
    assert compacted == 1
    assert len(pool.residue_metadata) == 2
    assert len(pool.state.atoms) == 6

# Total: 70 test functions covering comprehensive GCMC insertion operations 
# (including integrated patterns, translation/rotation, cavity bias, log-space stability, and pool management)

# Support direct execution
if __name__ == "__main__":
    import pytest
    pytest.main([__file__, "-v"])
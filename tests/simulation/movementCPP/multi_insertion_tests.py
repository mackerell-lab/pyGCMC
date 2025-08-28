# tests/simulation/movementCPP/multi_insertion_tests.py
"""
Multi-insertion CBMC tests collection
Imports all multi-insertion test functions from extracted modules
"""

import pytest

# Import basic multi-insertion tests
from .multi_insertion_basic_funcs import (
    test_multi_insertion_basic,
    test_mproposal_scaling,
    test_parallel_insertion_count
)

from .multi_insertion_calculation_funcs import (
    test_region_volume_calculation,
    test_cbmc_weight_calculation
)

from .multi_insertion_advanced_funcs import (
    test_multi_insertion_with_cavity_bias,
    test_multi_insertion_performance
)

from .multi_insertion_balance_funcs import (
    test_multi_insertion_detailed_balance,
    test_multi_insertion_with_types
)

from .multi_insertion_thread_funcs import (
    test_multi_insertion_thread_safety
)

# Import improved multi-insertion tests
from .multi_insertion_improved_core_funcs import (
    test_multi_insertion_basic_strict,
    test_acceptance_rate_statistics,
    test_parallel_insertion_performance
)

from .multi_insertion_improved_balance_funcs import (
    test_detailed_balance_strict,
    test_cbmc_weight_validity
)

from .multi_insertion_improved_validation_funcs import (
    test_reproducibility_with_seed,
    test_cavity_bias_with_multi_insertion
)

# Import robustness tests
from .multi_insertion_robustness_monotonicity_funcs import (
    test_acceptance_rate_monotonicity_chemical_potential,
    test_acceptance_rate_monotonicity_volume
)

from .multi_insertion_robustness_energy_funcs import (
    test_energy_consistency_single_vs_system,
    test_region_independence_boundary
)

from .multi_insertion_robustness_determinism_funcs import (
    test_deterministic_config_selection,
    test_fallback_path_unit_consistency,
    test_energy_range_parameterization,
    test_openmp_determinism
)


# Optional: Organize tests into classes for better test discovery
class TestMultiInsertionBasic:
    """Basic multi-insertion tests"""
    def test_basic(self):
        test_multi_insertion_basic()
    
    def test_mproposal(self):
        test_mproposal_scaling()
    
    def test_parallel_count(self):
        test_parallel_insertion_count()
    
    def test_region_volume(self):
        test_region_volume_calculation()
    
    def test_cbmc_weight(self):
        test_cbmc_weight_calculation()


class TestMultiInsertionAdvanced:
    """Advanced multi-insertion tests"""
    def test_cavity_bias(self):
        test_multi_insertion_with_cavity_bias()
    
    def test_performance(self):
        test_multi_insertion_performance()
    
    def test_detailed_balance(self):
        test_multi_insertion_detailed_balance()
    
    def test_with_types(self):
        test_multi_insertion_with_types()
    
    def test_thread_safety(self):
        test_multi_insertion_thread_safety()


class TestMultiInsertionImproved:
    """Improved multi-insertion tests with strict validation"""
    def test_basic_strict(self):
        test_multi_insertion_basic_strict()
    
    def test_acceptance_statistics(self):
        test_acceptance_rate_statistics()
    
    def test_parallel_performance(self):
        test_parallel_insertion_performance()
    
    def test_detailed_balance_strict(self):
        test_detailed_balance_strict()
    
    def test_cbmc_validity(self):
        test_cbmc_weight_validity()
    
    def test_reproducibility(self):
        test_reproducibility_with_seed()
    
    def test_cavity_integration(self):
        test_cavity_bias_with_multi_insertion()


class TestMultiInsertionRobustness:
    """Robustness tests for multi-insertion"""
    def test_monotonicity_chemical_potential(self):
        test_acceptance_rate_monotonicity_chemical_potential()
    
    def test_monotonicity_volume(self):
        test_acceptance_rate_monotonicity_volume()
    
    def test_energy_consistency(self):
        test_energy_consistency_single_vs_system()
    
    def test_region_independence(self):
        test_region_independence_boundary()
    
    def test_deterministic_selection(self):
        test_deterministic_config_selection()
    
    def test_fallback_consistency(self):
        test_fallback_path_unit_consistency()
    
    def test_energy_parameterization(self):
        test_energy_range_parameterization()
    
    def test_openmp(self):
        test_openmp_determinism()
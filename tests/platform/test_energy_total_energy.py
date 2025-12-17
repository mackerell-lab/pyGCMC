# tests/platform/test_energy_total_energy.py
"""
Energy total helpers - Main Entry Point

This file imports total-energy convention tests from modular sub-files.
"""

from energy.total_energy_unique_pairs import (
    test_total_energy_unique_pairs_direct_matches_half_residue_sum,
    test_total_energy_unique_pairs_ewald_uses_ewald_total_plus_half_vdw,
)


if __name__ == "__main__":
    import pytest

    pytest.main([__file__])


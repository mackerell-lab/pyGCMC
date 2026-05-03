# tests/io/test_topParser_new.py
"""
TOP Parser Tests - Main Entry Point

This file imports all TOP parser tests from modular sub-files.
Run: pytest tests/io/test_topParser.py

Modular structure:
- basic_parsing.py: Basic TOP parsing tests (2 functions)
- molecular_structures.py: Molecular structure tests (3 functions)
- topology_terms.py: Topology terms tests (3 functions)
- error_handling.py: Error handling tests (2 functions)
- gcmc_systems.py: GCMC system topology tests (4 functions)
- enhanced_4wp7_molecular.py: Enhanced 4wp7 molecular tests (3 functions)
- enhanced_4wp7_advanced.py: Enhanced 4wp7 advanced tests (2 functions)
"""

import os
import pytest

@pytest.fixture
def test_data_dir():
    """Get the path to the test data directory."""
    current_dir = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(current_dir, "..", "data")

# Basic parsing functionality
from topParser.basic_parsing import (
    test_parse_protein_top,
    test_parse_step1_top,
    test_parse_top_string_with_isolated_temp_file
)

# Molecular structures
from topParser.molecular_structures import (
    test_parse_solvent_top,
    test_parse_benzene_top,
    test_parse_propane_top
)

# Topology terms
from topParser.topology_terms import (
    test_parse_bonds,
    test_parse_angles,
    test_parse_dihedrals
)

# Error handling
from topParser.error_handling import (
    test_parse_nonexistent_file,
    test_parse_invalid_top
)

# GCMC systems testing
from topParser.gcmc_systems import (
    test_parse_4wp7_gcmc_topology,
    test_4wp7_topology_atom_types,
    test_4wp7_force_field_includes,
    test_4wp7_parser_limitations
)

# Enhanced 4wp7 tests - molecular structure
from topParser.enhanced_4wp7_molecular import (
    test_4wp7_gcmc_molecules_detailed,
    test_4wp7_topology_connectivity,
    test_4wp7_water_structure_validation
)

# Enhanced 4wp7 tests - advanced features
from topParser.enhanced_4wp7_advanced import (
    test_4wp7_preprocessor_directives,
    test_4wp7_error_handling
)

# Support direct execution for testing
if __name__ == "__main__":
    import pytest
    pytest.main([__file__])

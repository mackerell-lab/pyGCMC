# tests/io/test_topParser_new.py
"""
TOP Parser Tests - Main Entry Point

This file imports all TOP parser tests from modular sub-files.
Run: pytest tests/io/test_topParser_new.py

Modular structure:
- basic_parsing.py: Basic TOP parsing tests (2 functions)
- molecular_structures.py: Molecular structure tests (3 functions)
- topology_terms.py: Topology terms tests (3 functions)
- error_handling.py: Error handling tests (2 functions)
"""

# Basic parsing functionality
from topParser.basic_parsing import (
    test_parse_protein_top,
    test_parse_step1_top
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

# Support direct execution for testing
if __name__ == "__main__":
    import pytest
    pytest.main([__file__])

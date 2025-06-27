# tests/io/test_psfParser_new.py
"""
PSF Parser Tests - Main Entry Point

This file imports all PSF parser tests from modular sub-files.
Run: pytest tests/io/test_psfParser_new.py

Modular structure:
- basic_parsing.py: Basic PSF parsing tests (4 functions)
- molecular_structures.py: Molecular structure tests (4 functions)
- topology_terms.py: Topology terms tests (4 functions)
- advanced_features.py: Advanced features tests (4 functions)
"""

# Basic parsing functionality
from psfParser.basic_parsing import (
    test_parse_protein_psf,
    test_parse_nonexistent_file,
    test_parse_invalid_psf,
    test_parse_out_of_order_psf
)

# Molecular structures
from psfParser.molecular_structures import (
    test_parse_solvent_psf,
    test_parse_benzene_psf,
    test_parse_propane_psf,
    test_parse_step1_psf
)

# Topology terms
from psfParser.topology_terms import (
    test_parse_bonds,
    test_parse_angles,
    test_parse_dihedrals,
    test_parse_impropers
)

# Advanced features
from psfParser.advanced_features import (
    test_parse_donors_acceptors,
    test_parse_cmap,
    test_parse_groups,
    test_parse_all_cmaps
)

# Support direct execution for testing
if __name__ == "__main__":
    import pytest
    pytest.main([__file__])

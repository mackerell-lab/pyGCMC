# tests/io/test_psfParser.py
"""
PSF Parser Tests - Main Entry Point

This file imports all PSF parser tests from modular sub-files.
Run: pytest tests/io/test_psfParser.py

Modular structure:
- basic_parsing.py: Basic PSF parsing tests (4 functions)
- molecular_structures.py: Molecular structure tests (4 functions)
- topology_terms.py: Topology terms tests (4 functions)
- advanced_features.py: Advanced features tests (4 functions)
- drude_particles.py: Drude particle tests (4 functions)
- drude_validation.py: Drude validation tests (4 functions)
- drude_extended_format.py: Drude extended format tests (4 functions)
- drude_polarizability_analysis.py: Drude polarizability tests (3 functions)
- drude_comprehensive_analysis.py: Drude comprehensive analysis tests (5 functions)
"""

import os
import pytest

@pytest.fixture
def test_data_dir():
    """Get the path to the test data directory."""
    current_dir = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(current_dir, "..", "data")

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

# Drude particles and lone pairs
from psfParser.drude_particles import (
    test_parse_drude_psf_basic,
    test_parse_drude_particles,
    test_parse_lone_pairs,
    test_parse_drude_charge_neutrality
)

# Drude system validation
from psfParser.drude_validation import (
    test_parse_drude_connectivity,
    test_parse_drude_hydrogen_bonding,
    test_parse_drude_atom_types,
    test_parse_drude_residue_composition
)

# Drude extended format
from psfParser.drude_extended_format import (
    test_drude_psf_extended_format_detection,
    test_drude_alpha_parameters_parsing,
    test_drude_thole_parameters_parsing,
    test_drude_extended_atom_line_format
)

# Drude polarizability analysis
from psfParser.drude_polarizability_analysis import (
    test_drude_polarizability_distribution,
    test_drude_residue_polarizability,
    test_drude_anisotropy_information
)

# Drude comprehensive analysis
from psfParser.drude_comprehensive_analysis import (
    test_enhanced_atom_type_distribution,
    test_thole_screening_analysis,
    test_enhanced_connectivity_analysis,
    test_hydrogen_bonding_analysis,
    test_psf_validation_summary
)

# Support direct execution for testing
if __name__ == "__main__":
    import pytest
    pytest.main([__file__])

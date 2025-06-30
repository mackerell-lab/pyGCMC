# tests/io/test_pdbParser.py
"""
PDB Parser Tests - Main Entry Point

This file imports all PDB parser tests from modular sub-files.
Run: pytest tests/io/test_pdbParser.py

Modular structure:
- basic_parsing.py: Basic PDB parsing tests (3 functions)
- structure_parsing.py: Special structure parsing tests (2 functions)  
- data_processing.py: Data processing and validation tests (4 functions)
- complex_structures.py: Complex structure parsing tests (4 functions)
- calculations.py: Calculation and analysis tests (1 function)
- gcmc_systems.py: GCMC system parsing tests (3 functions)
"""

# Basic parsing functionality
from pdbParser.basic_parsing import (
    test_parse_simple_pdb,
    test_parse_hetatm,
    test_parse_ter
)

# Special structure parsing
from pdbParser.structure_parsing import (
    test_parse_secondary_structure,
    test_parse_ssbond
)

# Data processing and validation
from pdbParser.data_processing import (
    test_parse_invalid_pdb,
    test_residue_atom_association,
    test_coordinate_parsing,
    test_occupancy_and_tempfactor
)

# Complex structure parsing
from pdbParser.complex_structures import (
    test_parse_protein_fragment,
    test_parse_crystal_info,
    test_parse_solvent_and_ligands,
    test_hydrogen_atoms
)

# Calculation and analysis
from pdbParser.calculations import (
    test_center_of_mass
)

# GCMC systems testing
from pdbParser.gcmc_systems import (
    test_parse_4wp7_gcmc_system,
    test_4wp7_system_residue_distribution,
    test_4wp7_coordinate_validation
)

# Support direct execution for testing
if __name__ == "__main__":
    import pytest
    pytest.main([__file__])
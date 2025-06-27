# tests/io/test_prmParser_new.py
"""
PRM Parser Tests - Main Entry Point

This file imports all PRM parser tests from modular sub-files.
Run: pytest tests/io/test_prmParser_new.py

Modular structure:
- basic_parsing.py: Basic PRM parsing tests (6 functions)
- error_handling.py: Error handling tests (4 functions)  
- multiple_files_basic.py: Basic multiple file tests (6 functions)
- multiple_files_complex.py: Complex multiple file tests (2 functions)
- charmm_protein.py: CHARMM protein force field tests (2 functions)
- charmm_forcefields.py: CHARMM force field interaction tests (4 functions)
"""

# Basic parsing functionality
from prmParser.basic_parsing import (
    test_parse_file_direct,
    test_parse_nonbonded_from_string,
    test_parse_nbfix_from_string,
    test_parse_from_file,
    test_parse_comments_and_empty_lines,
    test_special_formatting
)

# Error handling
from prmParser.error_handling import (
    test_invalid_file,
    test_invalid_file_old_api,
    test_invalid_atom_type,
    test_malformed_parameters
)

# Multiple file operations - basic
from prmParser.multiple_files_basic import (
    test_multiple_nbfix_combinations,
    test_multiple_file_parsing,
    test_random_parameter_combinations,
    test_parse_multiple_files,
    test_parse_multiple_files_with_invalid
)

# Multiple file operations - complex
from prmParser.multiple_files_complex import (
    test_multiple_parameter_files,
    test_prm_and_str_files
)

# CHARMM protein force field tests
from prmParser.charmm_protein import (
    test_charmm_prm_files,
    test_cgenff_prm_file
)

# CHARMM force field interaction tests
from prmParser.charmm_forcefields import (
    test_ion_ligand_nbfix,
    test_heterocyclic_parameters,
    test_nucleic_parameters,
    test_cross_forcefield_compatibility
)

# Support direct execution for testing
if __name__ == "__main__":
    import pytest
    pytest.main([__file__])
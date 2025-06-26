#!/usr/bin/env python
# tests/io/test_inpParser.py
"""
INP Parser Tests - Main Entry Point

This file imports all INP parser tests from modular sub-files.
Run: pytest tests/io/test_inpParser.py

Modular structure:
- file_parsing.py: File-based parsing tests (2 functions)
- string_parsing.py: String-based parsing tests (1 function)  
- validation.py: Validation and error handling tests (1 function)
"""

# File parsing functionality
from inp.file_parsing import (
    test_read_gcmc_inp,
    test_inp_parser_file
)

# String parsing functionality
from inp.string_parsing import (
    test_inp_parser_string
)

# Validation and error handling
from inp.validation import (
    test_inp_parser_validation
)

# Support direct execution for testing
if __name__ == "__main__":
    import pytest
    pytest.main([__file__])
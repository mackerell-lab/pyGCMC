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
# These tests read INP data from files on disk and verify the parsed
# Parameter object matches the expected reference values.
from inpParser.file_parsing import (
    test_read_gcmc_inp,
    test_inp_parser_file
)

# String parsing functionality
# Verifies that INP content supplied as a raw string is parsed correctly
# without relying on the filesystem.
from inpParser.string_parsing import (
    test_inp_parser_string
)

# Validation and error handling
# Ensures malformed or incomplete INP content raises the expected RuntimeError
# with a clear, descriptive message.
from inpParser.validation import (
    test_inp_parser_validation
)

# Support direct execution for testing
if __name__ == "__main__":
    import pytest
    pytest.main([__file__])
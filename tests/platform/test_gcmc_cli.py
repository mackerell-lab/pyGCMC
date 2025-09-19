"""
GCMC CLI Tests - Main Entry Point

This file imports all CLI tests from modular sub-files.
Run: pytest tests/platform/test_gcmc_cli.py

Test categories:
1. Basic tests - Help, missing input, basic run
2. Output tests - Verbose output, file generation
3. Validation tests - Deterministic seed, invalid input, parameter validation
"""

import pytest

# Import fixtures from cli/conftest.py for all tests in this file
from cli.conftest import (
    gcmc_cpu,
    test_data_dir,
    temp_dir
)

# Basic CLI tests (3 functions)
from cli.basic_tests import (
    test_gcmc_cpu_help,
    test_gcmc_cpu_missing_inp,
    test_gcmc_cpu_basic_run
)

# Output and file generation tests (2 functions)
from cli.output_tests import (
    test_gcmc_cpu_verbose_output,
    test_gcmc_cpu_output_files
)

# Validation and error handling tests (3 functions)
from cli.validation_tests import (
    test_gcmc_cpu_deterministic_seed,
    test_gcmc_cpu_invalid_inp,
    test_gcmc_cpu_parameter_validation
)
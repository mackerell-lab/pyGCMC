# tests/io/inp/test_validation.py
"""INP Parser validation and error handling tests."""

import pytest
import pygcmc


def test_inp_parser_validation():
    """Test validation and error handling for INP parser."""
    # Test missing required parameters
    with pytest.raises(RuntimeError, match="Missing required parameter: top"):
        pygcmc.io.INPParser.parse_string("pdb:test.pdb")

    with pytest.raises(RuntimeError, match="Missing required parameter: pdb"):
        pygcmc.io.INPParser.parse_string("top:test.top")

    # Test inconsistent fragment parameters
    with pytest.raises(RuntimeError, match="Inconsistent fragment parameters"):
        pygcmc.io.INPParser.parse_string("""
top:test.top
pdb:test.pdb
fragname:benx prpx
fragconc:0.25
fragmuex:-0.79 1.96
""")

    # Test invalid space parameters
    with pytest.raises(RuntimeError, match="Invalid grid_dx"):
        pygcmc.io.INPParser.parse_string("""
top:test.top
pdb:test.pdb
grid_dx:-1.0
""")

    with pytest.raises(RuntimeError, match="Invalid cutoff"):
        pygcmc.io.INPParser.parse_string("""
top:test.top
pdb:test.pdb
cutoff:-12.0
""")
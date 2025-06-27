# tests/io/topParser/error_handling.py
"""TOP Parser error handling tests."""

import os
import pytest
from pygcmc.io import TOPParser
from pygcmc.model import Topology, TopologyResidue, TopologyAtom


def test_parse_nonexistent_file(test_data_dir):
    """Test parsing a non-existent topology file."""
    top_file = os.path.join(test_data_dir, "nonexistent.top")
    parser = TOPParser()
    topology = Topology()
    
    assert not parser.parse_to_topology(top_file, topology), "Should fail for non-existent file"


def test_parse_invalid_top(test_data_dir, tmp_path):
    """Test parsing an invalid topology file."""
    # Create an invalid topology file
    invalid_top = tmp_path / "invalid.top"
    with open(invalid_top, "w") as f:
        f.write("This is not a topology file\n")
    
    parser = TOPParser()
    topology = Topology()
    assert not parser.parse_to_topology(str(invalid_top), topology), "Should fail for invalid topology file"

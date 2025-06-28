# tests/system/molecularSystem/helpers.py

import pytest
import os
import pygcmc
import math
from pygcmc.model import Topology
from pygcmc.io import TOPParser
import random

# Get the directory containing test data files
TEST_DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "data")

@pytest.fixture
def test_structure():
    """Fixture to create a test Structure object."""
    pdb_path = os.path.join(TEST_DATA_DIR, "test.pdb")
    return pygcmc.PDBParser.parse_file(pdb_path)

@pytest.fixture
def test_topology():
    """Fixture to create a test Topology object."""
    top_path = os.path.join(TEST_DATA_DIR, "test.top")
    return TOPParser.parse_file(top_path)
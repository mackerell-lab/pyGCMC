# tests/system/MonteCarloSystem/helpers.py

import pytest
import pygcmc
import os
import math
from pygcmc.model import Molecular
from pygcmc.io import PDBParser, TOPParser

# Get the directory containing test data files
TEST_DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "data")

@pytest.fixture
def molecular_system():
    """Create a test molecular system from PDB and TOP files."""
    pdb_path = os.path.join(TEST_DATA_DIR, "test.pdb")
    top_path = os.path.join(TEST_DATA_DIR, "test.top")
    
    structure = pygcmc.PDBParser.parse_file(pdb_path)
    topology = pygcmc.TOPParser.parse_file(top_path)
    
    mol_system = pygcmc.MolecularSystem()
    return mol_system.combine(structure, topology)

def assert_arrays_almost_equal(arr1, arr2, tol=1e-6):
    """Compare two arrays for approximate equality."""
    if len(arr1) != len(arr2):
        return False
    return all(abs(a - b) < tol for a, b in zip(arr1, arr2))

@pytest.fixture
def charmm_ff():
    """Create a test CHARMM force field."""
    ff = pygcmc.ForceField()
    
    # Load CHARMM force field files using the correct method
    pygcmc.PRMParser.parse_file_to_forcefield(os.path.join(TEST_DATA_DIR, "par_all36_cgenff.prm"), ff)
    pygcmc.PRMParser.parse_file_to_forcefield(os.path.join(TEST_DATA_DIR, "par_all36m_prot.prm"), ff)
    pygcmc.PRMParser.parse_file_to_forcefield(os.path.join(TEST_DATA_DIR, "silcs.str"), ff)
    pygcmc.PRMParser.parse_file_to_forcefield(os.path.join(TEST_DATA_DIR, "toppar_water_ions.str"), ff)
    
    return ff
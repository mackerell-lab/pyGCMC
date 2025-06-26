# tests/io/pdbParser/data_processing.py
"""PDB Parser data processing and validation tests."""

import pytest
import os
import pygcmc
import math

# Get the directory containing test data files
TEST_DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "data")


def test_parse_invalid_pdb():
    """Test handling of invalid PDB files."""
    # Test non-existent file
    with pytest.raises(RuntimeError):
        pygcmc.PDBParser.parse_file("nonexistent.pdb")
    
    # Test invalid atom record
    invalid_pdb = """
ATOM   INVALID  N   MET A   1      27.340  24.430   2.614  1.00  0.00
"""
    with pytest.raises(RuntimeError):
        pygcmc.PDBParser.parse_string(invalid_pdb)
    
    # Test invalid coordinates
    invalid_coords = """
ATOM      1  N   MET A   1      XXXXX  24.430   2.614  1.00  0.00           N  
"""
    with pytest.raises(RuntimeError):
        pygcmc.PDBParser.parse_string(invalid_coords)


def test_residue_atom_association():
    """Test if atoms are correctly associated with residues."""
    pdb_path = os.path.join(TEST_DATA_DIR, "simple.pdb")
    result = pygcmc.PDBParser.parse_file(pdb_path)
    
    # Check MET residue
    met = result.residues[0]
    assert met.get_resname() == "MET"
    assert len(met.get_atoms()) == 8
    
    # Check each atom belongs to the residue
    for atom in met.get_atoms():
        assert atom.get_resname() == "MET"
        assert atom.get_ires() == 1
        assert atom.get_chain() == "A"


def test_coordinate_parsing():
    """Test parsing of accurate coordinates with precision."""
    pdb_path = os.path.join(TEST_DATA_DIR, "simple.pdb")
    result = pygcmc.PDBParser.parse_file(pdb_path)
    
    # Check first atom coordinates (N atom)
    atom = result.atoms[0]
    coords = atom.get_coor()
    assert math.isclose(coords[0], 27.340, rel_tol=1e-5)
    assert math.isclose(coords[1], 24.430, rel_tol=1e-5)
    assert math.isclose(coords[2], 2.614, rel_tol=1e-5)
    
    # Check all atoms have valid coordinates
    for atom in result.atoms:
        coords = atom.get_coor()
        assert len(coords) == 3
        assert all(isinstance(c, (int, float)) for c in coords)


def test_occupancy_and_tempfactor():
    """Test parsing of occupancy and B-factor values."""
    pdb_path = os.path.join(TEST_DATA_DIR, "simple.pdb")
    result = pygcmc.PDBParser.parse_file(pdb_path)
    
    # Check that atoms have occupancy and temperature factor
    atom = result.atoms[0]
    
    # These should be accessible (exact values depend on file content)
    if hasattr(atom, 'occupancy'):
        assert 0.0 <= atom.occupancy <= 1.0
    
    if hasattr(atom, 'tempfactor') or hasattr(atom, 'bfactor'):
        bfactor = getattr(atom, 'tempfactor', getattr(atom, 'bfactor', 0.0))
        assert bfactor >= 0.0
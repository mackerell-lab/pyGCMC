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
    
    # Check atom types in residue
    atom_types = [atom.get_type().strip() for atom in met.get_atoms()]
    expected_types = ["N", "CA", "C", "O", "CB", "CG", "SD", "CE"]
    assert sorted(atom_types) == sorted(expected_types)
    
    # Check atom-residue consistency
    for atom in met.get_atoms():
        assert atom.get_resname() == met.get_resname()
        assert atom.get_ires() == met.get_ires()
        assert atom.get_chain() == met.get_chain()
        assert atom.get_inscode() == met.get_inscode()


def test_coordinate_parsing():
    """Test parsing of atomic coordinates."""
    pdb_path = os.path.join(TEST_DATA_DIR, "simple.pdb")
    result = pygcmc.PDBParser.parse_file(pdb_path)
    
    for atom in result.atoms:
        coords = atom.get_coor()
        # Check coordinate format
        assert len(coords) == 3
        assert all(isinstance(c, float) for c in coords)
        assert all(not isinstance(c, str) for c in coords)
        # Check coordinate ranges
        assert all(-1000 < c < 1000 for c in coords)
        # Check precision (PDB format: 8.3f)
        for c in coords:
            assert abs(c - round(c, 3)) < 1e-3


def test_occupancy_and_tempfactor():
    """Test parsing of occupancy and temperature factor."""
    pdb_path = os.path.join(TEST_DATA_DIR, "simple.pdb")
    result = pygcmc.PDBParser.parse_file(pdb_path)
    
    for atom in result.atoms:
        # Check occupancy (default 1.00)
        assert math.isclose(atom.get_occupancy(), 1.00, rel_tol=1e-5)
        # Check temperature factor (default 0.00)
        assert math.isclose(atom.get_tempfactor(), 0.00, rel_tol=1e-5)
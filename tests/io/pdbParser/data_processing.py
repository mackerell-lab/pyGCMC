# tests/io/pdb/data_processing.py
"""PDB Parser data processing and validation tests."""

import pytest
import os
import pygcmc
import math

# Get the directory containing test data files
TEST_DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "data")


def test_parse_invalid_pdb():
    """Test handling of invalid or malformed PDB files."""
    # Test non-existent file
    with pytest.raises((FileNotFoundError, IOError)):
        pygcmc.PDBParser.parse_file("nonexistent.pdb")
    
    # Test empty file
    import tempfile
    with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f:
        f.write("")
        temp_path = f.name
    
    try:
        result = pygcmc.PDBParser.parse_file(temp_path)
        # Empty file should result in empty structure
        assert len(result.atoms) == 0
        assert len(result.residues) == 0
    finally:
        os.unlink(temp_path)
    
    # Test malformed ATOM line
    with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f:
        f.write("ATOM      1  N   ALA A   1      invalid coordinates\n")
        temp_path = f.name
    
    try:
        # Parser should handle malformed lines gracefully
        result = pygcmc.PDBParser.parse_file(temp_path)
        # Depending on implementation, might skip invalid lines or raise exception
    except (ValueError, RuntimeError):
        # This is acceptable behavior for malformed input
        pass
    finally:
        os.unlink(temp_path)


def test_residue_atom_association():
    """Test correct association between residues and atoms."""
    pdb_path = os.path.join(TEST_DATA_DIR, "simple.pdb")
    result = pygcmc.PDBParser.parse_file(pdb_path)
    
    # Check residue-atom relationships
    assert len(result.residues) == 1
    residue = result.residues[0]
    
    # Check residue properties
    assert residue.get_name() == "MET"
    assert residue.get_chain() == "A"
    assert residue.get_ires() == 1
    
    # Check that all atoms belong to the residue
    residue_atoms = residue.get_atoms()
    assert len(residue_atoms) == 8  # MET has 8 atoms
    
    # Verify atom-residue consistency
    for atom in residue_atoms:
        assert atom.get_resname() == residue.get_name()
        assert atom.get_chain() == residue.get_chain()
        assert atom.get_ires() == residue.get_ires()


def test_coordinate_parsing():
    """Test accurate parsing of atomic coordinates."""
    pdb_path = os.path.join(TEST_DATA_DIR, "simple.pdb")
    result = pygcmc.PDBParser.parse_file(pdb_path)
    
    # Test specific coordinates from simple.pdb
    # These values should match the actual file content
    atom = result.atoms[0]  # First atom (N)
    coords = atom.get_coor()
    
    # Check coordinate precision
    assert math.isclose(coords[0], 27.340, abs_tol=1e-3)
    assert math.isclose(coords[1], 24.430, abs_tol=1e-3)
    assert math.isclose(coords[2], 2.614, abs_tol=1e-3)
    
    # Check that all atoms have valid coordinates
    for atom in result.atoms:
        coords = atom.get_coor()
        assert len(coords) == 3
        # Coordinates should be finite numbers
        assert all(math.isfinite(c) for c in coords)


def test_occupancy_and_tempfactor():
    """Test parsing of occupancy and temperature factor values."""
    pdb_path = os.path.join(TEST_DATA_DIR, "simple.pdb")
    result = pygcmc.PDBParser.parse_file(pdb_path)
    
    # Check occupancy and B-factor for first atom
    atom = result.atoms[0]
    
    # Default values if not specified in PDB
    occupancy = getattr(atom, 'occupancy', 1.0)
    bfactor = getattr(atom, 'tempfactor', 0.0)
    
    # Occupancy should be between 0 and 1
    assert 0.0 <= occupancy <= 1.0
    
    # B-factor should be non-negative
    assert bfactor >= 0.0
    
    # Check that all atoms have valid occupancy/B-factor values
    for atom in result.atoms:
        if hasattr(atom, 'occupancy'):
            assert 0.0 <= atom.occupancy <= 1.0
        if hasattr(atom, 'tempfactor'):
            assert atom.tempfactor >= 0.0
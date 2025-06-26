# tests/io/pdb/basic_parsing.py
"""PDB Parser basic parsing tests."""

import pytest
import os
import pygcmc
import math

# Get the directory containing test data files
TEST_DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "data")


def test_parse_simple_pdb():
    """Test parsing a simple PDB file with basic ATOM records."""
    pdb_path = os.path.join(TEST_DATA_DIR, "simple.pdb")
    result = pygcmc.PDBParser.parse_file(pdb_path)
    
    # Check if atoms were parsed
    assert len(result.atoms) == 8  # MET has 8 atoms
    assert len(result.residues) == 1  # One MET residue
    
    # Check first atom's properties (N)
    atom = result.atoms[0]
    assert atom.get_bynu() == 1
    assert atom.get_type() == "N"
    assert atom.get_resname() == "MET"
    assert atom.get_chain() == "A"
    assert atom.get_ires() == 1
    
    # Check coordinates
    coords = atom.get_coor()
    assert len(coords) == 3
    assert math.isclose(coords[0], 27.340, rel_tol=1e-5)
    assert math.isclose(coords[1], 24.430, rel_tol=1e-5)
    assert math.isclose(coords[2], 2.614, rel_tol=1e-5)


def test_parse_hetatm():
    """Test parsing HETATM records."""
    pdb_path = os.path.join(TEST_DATA_DIR, "hetatm.pdb")
    result = pygcmc.PDBParser.parse_file(pdb_path)
    
    # Check counts
    protein_atoms = [atom for atom in result.atoms if not atom.is_hetatm()]
    water_atoms = [atom for atom in result.atoms if atom.is_hetatm()]
    assert len(protein_atoms) == 5  # ALA has 5 atoms
    assert len(water_atoms) == 3    # 3 water molecules
    
    # Check water properties
    water = water_atoms[0]
    assert water.is_hetatm()
    assert water.get_resname() == "HOH"
    assert water.get_type() == "O"
    
    # Check water coordinates (from original test)
    coords = water.get_coor()
    assert math.isclose(coords[0], 15.168, rel_tol=1e-5)
    assert math.isclose(coords[1], 20.391, rel_tol=1e-5)
    assert math.isclose(coords[2], 18.649, rel_tol=1e-5)


def test_parse_ter():
    """Test parsing TER records and chain termination."""
    pdb_path = os.path.join(TEST_DATA_DIR, "multichain.pdb")
    result = pygcmc.PDBParser.parse_file(pdb_path)
    
    # Check chain counts
    chains = set(atom.get_chain() for atom in result.atoms)
    assert chains == {"A", "B", "C"}
    
    # Check residues in different chains
    chain_residues = {}
    for res in result.residues:
        chain = res.get_chain()
        if chain not in chain_residues:
            chain_residues[chain] = []
        chain_residues[chain].append(res)
    
    # Verify chain contents
    assert len(chain_residues["A"]) == 1  # GLY
    assert len(chain_residues["B"]) == 1  # ALA
    assert len(chain_residues["C"]) == 1  # VAL
    
    # Check residue types
    assert chain_residues["A"][0].get_resname() == "GLY"
    assert chain_residues["B"][0].get_resname() == "ALA"
    assert chain_residues["C"][0].get_resname() == "VAL"
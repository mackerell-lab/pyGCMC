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
    
    # Check residue indices
    water_residues = [res for res in result.residues if res.get_name() == "HOH"]
    assert len(water_residues) == 1
    
    water_res = water_residues[0]
    assert water_res.get_name() == "HOH"
    assert len(water_res.get_atoms()) == 3


def test_parse_ter():
    """Test parsing TER records that terminate chains."""
    pdb_path = os.path.join(TEST_DATA_DIR, "multichain.pdb")
    result = pygcmc.PDBParser.parse_file(pdb_path)
    
    # Check that we have multiple chains
    chains = set(atom.get_chain() for atom in result.atoms)
    assert len(chains) >= 2  # At least two chains
    
    # Check that TER records properly terminate chains
    chain_a_atoms = [atom for atom in result.atoms if atom.get_chain() == "A"]
    chain_b_atoms = [atom for atom in result.atoms if atom.get_chain() == "B"]
    
    assert len(chain_a_atoms) > 0
    assert len(chain_b_atoms) > 0
    
    # Verify that chain termination is handled correctly
    # The exact test depends on the multichain.pdb content
    assert len(result.atoms) == len(chain_a_atoms) + len(chain_b_atoms)
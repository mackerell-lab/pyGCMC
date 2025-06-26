# tests/io/pdbParser/complex_structures.py
"""PDB Parser complex structure parsing tests."""

import pytest
import os
import pygcmc
import math

# Get the directory containing test data files
TEST_DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "data")


def test_parse_protein_fragment():
    """Test parsing protein fragment."""
    pdb_path = os.path.join(TEST_DATA_DIR, "test.pdb")
    result = pygcmc.PDBParser.parse_file(pdb_path)
    
    # Check total number of protein atoms
    protein_atoms = [atom for atom in result.atoms if not atom.is_hetatm()]
    assert len(protein_atoms) == 196, f"Expected 196 protein atoms, got {len(protein_atoms)}"
    
    # Check number of residues (only protein residues)
    protein_residues = set((atom.get_resname(), atom.get_ires()) 
                          for atom in protein_atoms 
                          if atom.get_resname() in ["ALA", "VAL", "PRO", "ASN", "GLN"])
    assert len(protein_residues) == 9, f"Expected 9 protein residues, got {len(protein_residues)}"
    
    # Check first residue (ALA 7)
    ala_atoms = [atom for atom in protein_atoms if atom.get_resname() == "ALA" and atom.get_ires() == 7]
    assert len(ala_atoms) == 12, f"Expected 12 atoms in ALA 7, got {len(ala_atoms)}"
    
    # Check first atom properties
    first_atom = protein_atoms[0]
    assert first_atom.get_type() == "N"
    assert first_atom.get_resname() == "ALA"
    assert first_atom.get_ires() == 7
    assert first_atom.get_chain() == " "
    assert not first_atom.is_hetatm()


def test_parse_crystal_info():
    """Test parsing crystallographic information."""
    pdb_path = os.path.join(TEST_DATA_DIR, "test.pdb")
    result = pygcmc.PDBParser.parse_file(pdb_path)
    
    # Check box dimensions
    assert hasattr(result, "box_dimensions")
    assert len(result.box_dimensions) == 6
    assert math.isclose(result.box_dimensions[0], 127.022, rel_tol=1e-5)  # a
    assert math.isclose(result.box_dimensions[1], 133.419, rel_tol=1e-5)  # b
    assert math.isclose(result.box_dimensions[2], 132.854, rel_tol=1e-5)  # c
    assert math.isclose(result.box_dimensions[3], 90.0, rel_tol=1e-5)     # alpha
    assert math.isclose(result.box_dimensions[4], 90.0, rel_tol=1e-5)     # beta
    assert math.isclose(result.box_dimensions[5], 90.0, rel_tol=1e-5)     # gamma


def test_parse_solvent_and_ligands():
    result = pygcmc.PDBParser.parse_file(os.path.join(TEST_DATA_DIR, "test.pdb"))
    
    # Check water molecules (SOL)
    water_atoms = [atom for atom in result.atoms if atom.get_resname() == "SOL"]
    assert len(water_atoms) == 30, f"Expected 30 water atoms (10 molecules), got {len(water_atoms)}"
    
    # Check BENX ligand (marked as ATOM in test.pdb)
    benx_atoms = [atom for atom in result.atoms if atom.get_resname() == "BENX"]
    assert len(benx_atoms) == 13, f"Expected 13 BENX atoms, got {len(benx_atoms)}"
    
    # Check PRPX ligand (marked as ATOM in test.pdb)
    prpx_atoms = [atom for atom in result.atoms if atom.get_resname() == "PRPX"]
    assert len(prpx_atoms) == 24, f"Expected 24 PRPX atoms (2 molecules), got {len(prpx_atoms)}"
    
    # Check first water molecule
    first_water = water_atoms[0]
    assert first_water.get_type() == "OW"
    assert first_water.get_resname() == "SOL"
    assert not first_water.is_hetatm()  # In test.pdb, these are ATOM records
    
    # Check first BENX molecule
    first_benx = benx_atoms[0]
    assert first_benx.get_type() == "CG"
    assert first_benx.get_resname() == "BENX"
    assert not first_benx.is_hetatm()  # In test.pdb, these are ATOM records
    
    # Check first PRPX molecule
    first_prpx = prpx_atoms[0]
    assert first_prpx.get_type() == "H11"
    assert first_prpx.get_resname() == "PRPX"
    assert not first_prpx.is_hetatm()  # In test.pdb, these are ATOM records


def test_hydrogen_atoms():
    pdb_path = os.path.join(TEST_DATA_DIR, "test.pdb")
    result = pygcmc.PDBParser.parse_file(pdb_path)
    
    # Get all atoms from first ALA residue (ALA 7)
    ala_atoms = [atom for atom in result.atoms if atom.get_resname() == "ALA" and atom.get_ires() == 7]
    h_atoms = [atom for atom in ala_atoms if atom.get_type().strip()[0] == "H"]
    
    assert len(h_atoms) == 7, f"Expected 7 hydrogen atoms in ALA 7, got {len(h_atoms)}"
    
    # Check types of hydrogen atoms
    h_types = set(atom.get_type().strip() for atom in h_atoms)
    expected_types = {"H1", "H2", "H3", "HA", "HB1", "HB2", "HB3"}
    assert h_types == expected_types, f"Expected H types {expected_types}, got {h_types}"
    
    # Check coordinates of first hydrogen
    h1_atom = next(atom for atom in h_atoms if atom.get_type().strip() == "H1")
    coords = h1_atom.get_coor()
    assert len(coords) == 3
    assert abs(coords[0] - 76.044) < 0.001
    assert abs(coords[1] - 92.324) < 0.001
    assert abs(coords[2] - 93.379) < 0.001
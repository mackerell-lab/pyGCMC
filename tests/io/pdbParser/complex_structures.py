# tests/io/pdb/complex_structures.py
"""PDB Parser complex structure parsing tests."""

import pytest
import os
import pygcmc

# Get the directory containing test data files
TEST_DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "data")


def test_parse_protein_fragment():
    """Test parsing protein fragments with multiple residues."""
    pdb_path = os.path.join(TEST_DATA_DIR, "test.pdb")
    result = pygcmc.PDBParser.parse_file(pdb_path)
    
    # Check that multiple residues were parsed
    assert len(result.residues) > 1
    assert len(result.atoms) > 8  # More than a single residue
    
    # Check residue variety
    residue_names = set(res.get_resname() for res in result.residues)
    assert len(residue_names) > 1  # Multiple residue types
    
    # Check chain continuity
    chains = set(atom.get_chain() for atom in result.atoms)
    assert len(chains) >= 1
    
    # Verify backbone atoms are present
    backbone_atoms = ["N", "CA", "C", "O"]
    found_backbone = set()
    
    for atom in result.atoms:
        if atom.get_type() in backbone_atoms:
            found_backbone.add(atom.get_type())
    
    # Should find at least some backbone atoms
    assert len(found_backbone) > 0
    
    # Check residue numbering
    residue_numbers = [res.get_ires() for res in result.residues]
    assert len(set(residue_numbers)) == len(residue_numbers)  # Unique numbering


def test_parse_solvent_and_ligands():
    """Test parsing files containing protein, solvent, and ligands."""
    pdb_path = os.path.join(TEST_DATA_DIR, "water.pdb")
    result = pygcmc.PDBParser.parse_file(pdb_path)
    
    # Separate protein and non-protein atoms
    protein_atoms = []
    water_atoms = []
    other_atoms = []
    
    for atom in result.atoms:
        if atom.get_resname() in ["HOH", "WAT", "TIP", "TIP3", "SPC"]:
            water_atoms.append(atom)
        elif atom.is_hetatm():
            other_atoms.append(atom)
        else:
            protein_atoms.append(atom)
    
    # Should have some water molecules
    assert len(water_atoms) > 0
    
    # Check water molecule structure
    water_residues = [res for res in result.residues 
                     if res.get_resname() in ["HOH", "WAT", "TIP", "TIP3", "SPC"]]
    
    if water_residues:
        water_res = water_residues[0]
        water_res_atoms = water_res.get_atoms()
        # Water should have 1-3 atoms (O, or O+H, or O+H+H)
        assert 1 <= len(water_res_atoms) <= 3
        
        # Check for oxygen atom
        oxygen_found = any(atom.get_type().startswith('O') for atom in water_res_atoms)
        assert oxygen_found
    
    # Total atoms should equal sum of categories
    total_atoms = len(protein_atoms) + len(water_atoms) + len(other_atoms)
    assert total_atoms == len(result.atoms)


def test_hydrogen_atoms():
    """Test proper handling of hydrogen atoms."""
    pdb_path = os.path.join(TEST_DATA_DIR, "simple.pdb")
    result = pygcmc.PDBParser.parse_file(pdb_path)
    
    # Look for hydrogen atoms
    hydrogen_atoms = []
    heavy_atoms = []
    
    for atom in result.atoms:
        atom_type = atom.get_type()
        if atom_type.startswith('H') or atom_type.startswith('h'):
            hydrogen_atoms.append(atom)
        else:
            heavy_atoms.append(atom)
    
    # Check that heavy atoms are present
    assert len(heavy_atoms) > 0
    
    # If hydrogens are present, verify their properties
    for h_atom in hydrogen_atoms:
        atom_type = h_atom.get_type()
        assert atom_type.startswith('H') or atom_type.startswith('h')
        
        # Hydrogen coordinates should be valid
        coords = h_atom.get_coor()
        assert len(coords) == 3
        assert all(isinstance(c, (int, float)) for c in coords)
        
        # Hydrogen should belong to a residue
        assert h_atom.get_resname() is not None
        assert h_atom.get_chain() is not None
    
    # Verify atom count consistency
    assert len(hydrogen_atoms) + len(heavy_atoms) == len(result.atoms)
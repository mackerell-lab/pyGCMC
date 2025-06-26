# tests/io/pdb/calculations.py
"""PDB Parser calculation and analysis tests."""

import pytest
import os
import pygcmc
import math

# Get the directory containing test data files
TEST_DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "data")


def test_center_of_mass():
    """Test center of mass calculation for parsed structures."""
    pdb_path = os.path.join(TEST_DATA_DIR, "simple.pdb")
    result = pygcmc.PDBParser.parse_file(pdb_path)
    
    # Check that we have atoms to work with
    assert len(result.atoms) > 0
    
    # Calculate center of mass manually for verification
    total_mass = 0.0
    com_x, com_y, com_z = 0.0, 0.0, 0.0
    
    for atom in result.atoms:
        # Use atomic mass approximations
        atom_type = atom.get_type()
        if atom_type.startswith('C'):
            mass = 12.01
        elif atom_type.startswith('N'):
            mass = 14.01
        elif atom_type.startswith('O'):
            mass = 16.00
        elif atom_type.startswith('S'):
            mass = 32.06
        elif atom_type.startswith('H'):
            mass = 1.008
        else:
            mass = 12.01  # Default to carbon
        
        coords = atom.get_coor()
        com_x += mass * coords[0]
        com_y += mass * coords[1]
        com_z += mass * coords[2]
        total_mass += mass
    
    if total_mass > 0:
        com_x /= total_mass
        com_y /= total_mass
        com_z /= total_mass
        
        # Center of mass should be within reasonable bounds
        # (roughly centered around the atomic coordinates)
        min_x = min(atom.get_coor()[0] for atom in result.atoms)
        max_x = max(atom.get_coor()[0] for atom in result.atoms)
        min_y = min(atom.get_coor()[1] for atom in result.atoms)
        max_y = max(atom.get_coor()[1] for atom in result.atoms)
        min_z = min(atom.get_coor()[2] for atom in result.atoms)
        max_z = max(atom.get_coor()[2] for atom in result.atoms)
        
        assert min_x <= com_x <= max_x
        assert min_y <= com_y <= max_y
        assert min_z <= com_z <= max_z
    
    # Test geometric center as simpler alternative
    geom_center_x = sum(atom.get_coor()[0] for atom in result.atoms) / len(result.atoms)
    geom_center_y = sum(atom.get_coor()[1] for atom in result.atoms) / len(result.atoms)
    geom_center_z = sum(atom.get_coor()[2] for atom in result.atoms) / len(result.atoms)
    
    # Geometric center should be finite
    assert math.isfinite(geom_center_x)
    assert math.isfinite(geom_center_y)
    assert math.isfinite(geom_center_z)
    
    # For a single residue, center should be reasonably close to first atom
    first_atom_coords = result.atoms[0].get_coor()
    distance = math.sqrt(
        (geom_center_x - first_atom_coords[0])**2 +
        (geom_center_y - first_atom_coords[1])**2 +
        (geom_center_z - first_atom_coords[2])**2
    )
    
    # Distance should be reasonable (less than 100 Angstroms for typical residues)
    assert distance < 100.0
# tests/io/pdbParser/enhanced_4wp7_tests.py
"""
4wp7 Structure Tests - Molecular integrity and crystal structure.

This file tests molecular integrity and crystal structure aspects of 4wp7 PDB file
including:

Tests include:
- Molecular integrity and structure validation
- Crystal cell and spatial distribution analysis
"""

import pytest
import os
import time
import math
from pygcmc.io import PDBParser

# Get the directory containing test data files
TEST_DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "data")


def test_4wp7_molecule_integrity():
    """Test molecular integrity and structure validation for GCMC molecules and water."""
    pdb_path = os.path.join(TEST_DATA_DIR, "4wp7", "4wp7_fixed_with_5l13_silcs.1.prod.74.rec.pdb")
    
    # Skip test if file doesn't exist
    if not os.path.exists(pdb_path):
        pytest.skip(f"Test file {pdb_path} not found")
    
    result = PDBParser.parse_file(pdb_path)
    assert result is not None, "Failed to parse 4wp7 PDB file"
    
    # Expected molecular structures based on raw PDB file analysis
    # Atom counts verified from: grep "^ATOM" *.pdb | awk '{print $4}' | sort | uniq -c
    # Molecules per type calculated by dividing total atoms by atoms per molecule
    expected_gcmc_molecules = {
        "BENX": {"atoms_per_molecule": 13, "description": "Benzene with hydrogens"},    # 3731 atoms total
        "PRPX": {"atoms_per_molecule": 12, "description": "Propane with explicit hydrogens"},  # 3504 atoms total
        "IMIA": {"atoms_per_molecule": 9, "description": "Imidazole"},                 # 2601 atoms total  
        "MAMY": {"atoms_per_molecule": 9, "description": "Methylammonium"},            # 2313 atoms total
        "DMEE": {"atoms_per_molecule": 9, "description": "Dimethyl ether"},            # 2169 atoms total
        "ACEY": {"atoms_per_molecule": 8, "description": "Acetate"},                   # 2048 atoms total
        "MEOH": {"atoms_per_molecule": 6, "description": "Methanol"},                  # 1662 atoms total
        "FORM": {"atoms_per_molecule": 6, "description": "Formamide"}                  # 1662 atoms total
    }
    
    # Count atoms per residue - now using actual residue objects from parser
    residue_atom_counts = {}
    residue_instances = {}
    
    # Use the residue objects from the parser instead of counting atoms manually
    for residue in result.residues:
        res_name = residue.get_resname()
        if res_name in expected_gcmc_molecules or res_name == "SOL":
            atom_count = residue.atom_count()
            res_id = f"{res_name}_{residue.get_ires()}_{id(residue)}"  # Include object ID to make unique
            
            residue_atom_counts[res_id] = atom_count
            if res_name not in residue_instances:
                residue_instances[res_name] = []
            residue_instances[res_name].append(res_id)
    
    # Verify GCMC molecule structures
    for mol_name, expected_info in expected_gcmc_molecules.items():
        if mol_name in residue_instances:
            # Check a sample of molecules of this type
            sample_molecules = residue_instances[mol_name][:5]  # Check first 5 instances
            
            for res_id in sample_molecules:
                actual_atoms = residue_atom_counts[res_id]
                expected_atoms = expected_info["atoms_per_molecule"]
                
                assert actual_atoms == expected_atoms, \
                    f"{mol_name} molecule {res_id} ({expected_info['description']}): " \
                    f"expected {expected_atoms} atoms, got {actual_atoms}"
    
    # Verify TIP3P water structure (SOL residues should have exactly 3 atoms)
    if "SOL" in residue_instances:
        water_sample = residue_instances["SOL"][:10]  # Check first 10 water molecules
        
        for res_id in water_sample:
            actual_atoms = residue_atom_counts[res_id]
            assert actual_atoms == 3, \
                f"Water molecule {res_id} should have 3 atoms (OW + 2×HW), got {actual_atoms}"
    
    print(f"✓ Molecular integrity validated: GCMC molecules and water structures correct")


def test_4wp7_crystal_cell_parsing():
    """Test crystal cell information parsing and coordinate validation."""
    pdb_path = os.path.join(TEST_DATA_DIR, "4wp7", "4wp7_fixed_with_5l13_silcs.1.prod.74.rec.pdb")
    
    # Skip test if file doesn't exist
    if not os.path.exists(pdb_path):
        pytest.skip(f"Test file {pdb_path} not found")
    
    result = PDBParser.parse_file(pdb_path)
    assert result is not None, "Failed to parse 4wp7 PDB file"
    
    # Read raw file to extract CRYST1 record
    cryst1_info = None
    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith("CRYST1"):
                # Parse CRYST1 record: a, b, c, alpha, beta, gamma
                parts = line.split()
                cryst1_info = {
                    'a': float(parts[1]),
                    'b': float(parts[2]), 
                    'c': float(parts[3]),
                    'alpha': float(parts[4]),
                    'beta': float(parts[5]),
                    'gamma': float(parts[6])
                }
                break
    
    assert cryst1_info is not None, "CRYST1 record not found in PDB file"
    
    # Verify expected crystal cell dimensions (from file header)
    expected_dims = {'a': 132.126, 'b': 132.213, 'c': 122.070}
    for dim, expected_val in expected_dims.items():
        actual_val = cryst1_info[dim]
        assert abs(actual_val - expected_val) < 0.01, \
            f"Crystal cell {dim}: expected {expected_val}, got {actual_val}"
    
    # Verify orthogonal cell (90 degree angles)
    for angle in ['alpha', 'beta', 'gamma']:
        assert abs(cryst1_info[angle] - 90.0) < 0.01, \
            f"Expected orthogonal cell, {angle} = {cryst1_info[angle]}"
    
    # Verify atom coordinates are within reasonable bounds
    x_coords = [atom.get_x() for atom in result.atoms]
    y_coords = [atom.get_y() for atom in result.atoms]
    z_coords = [atom.get_z() for atom in result.atoms]
    
    # Coordinates should generally be within cell dimensions (allowing for some flexibility)
    x_range = max(x_coords) - min(x_coords)
    y_range = max(y_coords) - min(y_coords)
    z_range = max(z_coords) - min(z_coords)
    
    assert x_range <= cryst1_info['a'] + 10, \
        f"X coordinate range {x_range:.1f} exceeds cell dimension a={cryst1_info['a']}"
    assert y_range <= cryst1_info['b'] + 10, \
        f"Y coordinate range {y_range:.1f} exceeds cell dimension b={cryst1_info['b']}"
    assert z_range <= cryst1_info['c'] + 10, \
        f"Z coordinate range {z_range:.1f} exceeds cell dimension c={cryst1_info['c']}"
    
    # Verify reasonable coordinate distribution (no clustering at origin)
    center_x = sum(x_coords) / len(x_coords)
    center_y = sum(y_coords) / len(y_coords)
    center_z = sum(z_coords) / len(z_coords)
    
    # Center should be roughly in middle of cell
    assert 20 <= center_x <= cryst1_info['a'] - 20, \
        f"X center {center_x:.1f} suggests poor distribution in cell"
    assert 20 <= center_y <= cryst1_info['b'] - 20, \
        f"Y center {center_y:.1f} suggests poor distribution in cell"
    assert 20 <= center_z <= cryst1_info['c'] - 20, \
        f"Z center {center_z:.1f} suggests poor distribution in cell"
    
    print(f"✓ Crystal cell validated: {cryst1_info['a']:.1f}×{cryst1_info['b']:.1f}×{cryst1_info['c']:.1f} Å³")



# tests/io/pdbParser/comparison_tests.py
"""
Comparison tests between direct Python PDB parsing and pygcmc parser.
Tests basic atom-level and structure-level comparisons.
"""

import pytest
import os
import math
import pygcmc
from .direct_parser_core import (
    DirectPDBAtom, read_pdb_direct, compare_atom_data, TEST_DATA_DIR
)


def test_4wp7_direct_vs_pygcmc_parsing():
    """Compare direct Python PDB parsing with pygcmc parser for 4wp7 system."""
    pdb_path = os.path.join(TEST_DATA_DIR, "4wp7", "4wp7_fixed_with_5l13_silcs.1.prod.74.rec.pdb")
    
    # Skip test if file doesn't exist
    if not os.path.exists(pdb_path):
        pytest.skip(f"Test file {pdb_path} not found")
    
    # Parse with both methods
    direct_atoms, direct_metadata = read_pdb_direct(pdb_path)
    pygcmc_result = pygcmc.PDBParser.parse_file(pdb_path)
    
    # Basic comparison
    assert len(direct_atoms) == len(pygcmc_result.atoms), \
        f"Atom count mismatch: direct={len(direct_atoms)}, pygcmc={len(pygcmc_result.atoms)}"
    
    print(f"Both parsers found {len(direct_atoms)} atoms")
    
    # Sample comparison (first 100 atoms for performance)
    sample_size = min(100, len(direct_atoms))
    mismatches = 0
    
    for i in range(sample_size):
        direct_atom = direct_atoms[i]
        pygcmc_atom = pygcmc_result.atoms[i]
        
        differences = compare_atom_data(direct_atom, pygcmc_atom)
        if differences:
            mismatches += 1
            if mismatches <= 5:  # Only report first 5 mismatches
                print(f"Atom {i+1} differences: {'; '.join(differences)}")
    
    # Allow small number of mismatches due to different parsing approaches
    mismatch_rate = mismatches / sample_size
    assert mismatch_rate < 0.05, f"Too many mismatches: {mismatches}/{sample_size} ({mismatch_rate:.1%})"
    
    print(f"Comparison: {mismatches}/{sample_size} atoms had differences ({mismatch_rate:.1%})")


def test_4wp7_coordinate_ranges_comparison():
    """Compare coordinate ranges and statistics between parsers."""
    pdb_path = os.path.join(TEST_DATA_DIR, "4wp7", "4wp7_fixed_with_5l13_silcs.1.prod.74.rec.pdb")
    
    # Skip test if file doesn't exist
    if not os.path.exists(pdb_path):
        pytest.skip(f"Test file {pdb_path} not found")
    
    # Parse with both methods
    direct_atoms, _ = read_pdb_direct(pdb_path)
    pygcmc_result = pygcmc.PDBParser.parse_file(pdb_path)
    
    # Calculate coordinate ranges for direct parsing
    if direct_atoms:
        direct_x_coords = [atom.x for atom in direct_atoms]
        direct_y_coords = [atom.y for atom in direct_atoms]
        direct_z_coords = [atom.z for atom in direct_atoms]
        
        direct_ranges = {
            'x': (min(direct_x_coords), max(direct_x_coords)),
            'y': (min(direct_y_coords), max(direct_y_coords)),
            'z': (min(direct_z_coords), max(direct_z_coords))
        }
    
    # Calculate coordinate ranges for pygcmc parsing
    if pygcmc_result.atoms:
        pygcmc_x_coords = [atom.get_coor()[0] for atom in pygcmc_result.atoms]
        pygcmc_y_coords = [atom.get_coor()[1] for atom in pygcmc_result.atoms]
        pygcmc_z_coords = [atom.get_coor()[2] for atom in pygcmc_result.atoms]
        
        pygcmc_ranges = {
            'x': (min(pygcmc_x_coords), max(pygcmc_x_coords)),
            'y': (min(pygcmc_y_coords), max(pygcmc_y_coords)),
            'z': (min(pygcmc_z_coords), max(pygcmc_z_coords))
        }
    
    # Compare coordinate ranges (should be identical or very close)
    tolerance = 1e-3
    for coord in ['x', 'y', 'z']:
        direct_min, direct_max = direct_ranges[coord]
        pygcmc_min, pygcmc_max = pygcmc_ranges[coord]
        
        assert math.isclose(direct_min, pygcmc_min, abs_tol=tolerance), \
            f"{coord} min mismatch: direct={direct_min:.3f}, pygcmc={pygcmc_min:.3f}"
        
        assert math.isclose(direct_max, pygcmc_max, abs_tol=tolerance), \
            f"{coord} max mismatch: direct={direct_max:.3f}, pygcmc={pygcmc_max:.3f}"
    
    print(f"Coordinate ranges match: X[{direct_ranges['x'][0]:.1f}, {direct_ranges['x'][1]:.1f}], "
          f"Y[{direct_ranges['y'][0]:.1f}, {direct_ranges['y'][1]:.1f}], "
          f"Z[{direct_ranges['z'][0]:.1f}, {direct_ranges['z'][1]:.1f}]")


def test_specific_atom_verification():
    """Verify specific known atoms between direct and pygcmc parsing."""
    pdb_path = os.path.join(TEST_DATA_DIR, "4wp7", "4wp7_fixed_with_5l13_silcs.1.prod.74.rec.pdb")
    
    # Skip test if file doesn't exist
    if not os.path.exists(pdb_path):
        pytest.skip(f"Test file {pdb_path} not found")
    
    # Parse with both methods
    direct_atoms, _ = read_pdb_direct(pdb_path)
    pygcmc_result = pygcmc.PDBParser.parse_file(pdb_path)
    
    # Find first ASP N atom (from line 6 in PDB file)
    # ATOM      1  N   ASP     8      42.526  94.418  83.282  1.00  0.00           N
    
    direct_asp_n = None
    for atom in direct_atoms:
        if (atom.atom_num == 1 and atom.atom_name == "N" and 
            atom.res_name == "ASP" and atom.res_seq == 8):
            direct_asp_n = atom
            break
    
    pygcmc_asp_n = None
    for atom in pygcmc_result.atoms:
        if (atom.get_bynu() == 1 and atom.get_type() == "N" and 
            atom.get_resname() == "ASP" and atom.get_ires() == 8):
            pygcmc_asp_n = atom
            break
    
    assert direct_asp_n is not None, "Could not find ASP N atom in direct parsing"
    assert pygcmc_asp_n is not None, "Could not find ASP N atom in pygcmc parsing"
    
    # Verify coordinates match exactly
    tolerance = 1e-3
    assert math.isclose(direct_asp_n.x, 42.526, abs_tol=tolerance), f"Direct X coordinate wrong: {direct_asp_n.x}"
    assert math.isclose(direct_asp_n.y, 94.418, abs_tol=tolerance), f"Direct Y coordinate wrong: {direct_asp_n.y}"
    assert math.isclose(direct_asp_n.z, 83.282, abs_tol=tolerance), f"Direct Z coordinate wrong: {direct_asp_n.z}"
    
    pygcmc_coords = pygcmc_asp_n.get_coor()
    assert math.isclose(pygcmc_coords[0], 42.526, abs_tol=tolerance), f"Pygcmc X coordinate wrong: {pygcmc_coords[0]}"
    assert math.isclose(pygcmc_coords[1], 94.418, abs_tol=tolerance), f"Pygcmc Y coordinate wrong: {pygcmc_coords[1]}"
    assert math.isclose(pygcmc_coords[2], 83.282, abs_tol=tolerance), f"Pygcmc Z coordinate wrong: {pygcmc_coords[2]}"
    
    # Compare the two atoms directly
    differences = compare_atom_data(direct_asp_n, pygcmc_asp_n)
    assert len(differences) == 0, f"ASP N atom differences: {'; '.join(differences)}"
    
    print("First ASP N atom verification passed - both parsers read identical data")
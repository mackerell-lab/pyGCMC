# tests/io/pdbParser/calculations.py
"""PDB Parser calculation and analysis tests."""

import pytest
import os
import pygcmc
import math

# Get the directory containing test data files
TEST_DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "data")


def test_center_of_mass():
    """Test center of mass calculation for residues."""
    pdb_path = os.path.join(TEST_DATA_DIR, "test.pdb")
    result = pygcmc.PDBParser.parse_file(pdb_path)
    
    # Test ALA 7 center of mass
    ala7 = next(res for res in result.residues if res.get_resname() == "ALA" and res.get_ires() == 7)
    
    # Debug information
    print("\nALA 7 atoms:")
    total_mass = 0.0
    weighted_pos = [0.0, 0.0, 0.0]
    for atom in ala7.get_atoms():
        coords = atom.get_coor()
        mass = atom.get_mass()
        print(f"Atom {atom.get_type()}: mass={mass:.3f}, coords=({coords[0]:.3f}, {coords[1]:.3f}, {coords[2]:.3f})")
        total_mass += mass
        for i in range(3):
            weighted_pos[i] += mass * coords[i]
    
    if total_mass > 0:
        for i in range(3):
            weighted_pos[i] /= total_mass
    print(f"\nCalculated COM: ({weighted_pos[0]:.3f}, {weighted_pos[1]:.3f}, {weighted_pos[2]:.3f})")
    
    # Get the COM from the residue
    com = ala7.get_center_of_mass()
    print(f"Residue COM: ({com[0]:.3f}, {com[1]:.3f}, {com[2]:.3f})")
    
    # Verify COM coordinates
    assert len(com) == 3
    # The COM should be somewhere in the middle of the residue
    assert 77.0 < com[0] < 79.0, f"COM x-coordinate {com[0]} not in range (77.0, 79.0)"
    assert 91.0 < com[1] < 93.0, f"COM y-coordinate {com[1]} not in range (91.0, 93.0)"
    assert 92.0 < com[2] < 94.0, f"COM z-coordinate {com[2]} not in range (92.0, 94.0)"
    
    # Test VAL 8 center of mass
    val8 = next(res for res in result.residues if res.get_resname() == "VAL" and res.get_ires() == 8)
    com = val8.get_center_of_mass()
    
    # Verify COM coordinates
    assert len(com) == 3
    assert 78.0 < com[0] < 80.0, f"VAL8 COM x-coordinate {com[0]} not in range (78.0, 80.0)"
    assert 94.0 < com[1] < 96.0, f"VAL8 COM y-coordinate {com[1]} not in range (94.0, 96.0)"
    assert 89.0 < com[2] < 91.0, f"VAL8 COM z-coordinate {com[2]} not in range (89.0, 91.0)"
    
    # Test water molecule (SOL)
    sol = next(res for res in result.residues if res.get_resname() == "SOL")
    com = sol.get_center_of_mass()
    
    # Verify COM coordinates for water
    assert len(com) == 3
    # The COM should be close to the oxygen atom position for water
    assert all(isinstance(x, float) for x in com)
    assert all(not math.isnan(x) for x in com)
    assert all(not math.isinf(x) for x in com)
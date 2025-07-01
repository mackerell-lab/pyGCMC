# tests/io/psfParser/drude_extended_format.py
"""
Drude PSF Parser Tests - Extended Format Verification

This file tests whether the PSF parser correctly reads extended Drude PSF format
including alpha (polarizability) and Thole parameters.

Based on analysis from analyze_drude_psf.sh:
- Total polarizable atoms: 3814
- Average polarizability: 0.8018 Å³
- Median polarizability: 0.8005 Å³
- Min polarizability: 0.1260 Å³
- Max polarizability: 2.3990 Å³
- Average Thole parameter: -1.220
- Thole range: -2.180 to -0.467

NOTE: These tests use pytest.xfail to document expected behavior
that is not yet implemented in the current PSF parser.
"""

import pytest
import os
from pygcmc.io import PSFParser
from pygcmc.model import Topology

# Get the directory containing test data files
TEST_DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "data")


def test_drude_psf_extended_format_detection():
    """Test if PSF parser detects extended Drude format from header."""
    psf_path = os.path.join(TEST_DATA_DIR, "4wp7", "4wp7_drude.psf")
    
    if not os.path.exists(psf_path):
        pytest.skip(f"Test file {psf_path} not found")
    
    # Read the first line to check format
    with open(psf_path, 'r') as f:
        header_line = f.readline().strip()
    
    # Verify the header contains DRUDE and EXT flags
    assert "PSF" in header_line, "PSF header not found"
    assert "EXT" in header_line, "EXT (extended) format flag not found"
    assert "DRUDE" in header_line, "DRUDE format flag not found"
    
    # The header should be: PSF EXT CMAP DRUDE XPLOR AUTOG
    expected_flags = ["PSF", "EXT", "CMAP", "DRUDE", "XPLOR"]
    for flag in expected_flags:
        assert flag in header_line, f"Expected flag '{flag}' not found in header"


def test_drude_alpha_parameters_parsing():
    """Test if PSF parser reads alpha (polarizability) parameters.
    
    Expected values from analyze_drude_psf.sh:
    - Total polarizable atoms: 3814
    - Average polarizability: 0.8018 Å³
    - Median polarizability: 0.8005 Å³
    - Min polarizability: 0.1260 Å³
    - Max polarizability: 2.3990 Å³
    
    Top atom types by average polarizability:
    - OE2: 2.3990 Å³ (n=32)
    - OD2: 2.3990 Å³ (n=26)
    - OE1: 1.9771 Å³ (n=50)
    - OD1: 1.9042 Å³ (n=45)
    """
    psf_path = os.path.join(TEST_DATA_DIR, "4wp7", "4wp7_drude.psf")
    
    if not os.path.exists(psf_path):
        pytest.skip(f"Test file {psf_path} not found")
    
    parser = PSFParser()
    topology = Topology()
    success = parser.parse_to_topology(psf_path, topology)
    assert success, "Failed to parse 4wp7_drude.psf"
    
    # Check if atoms have alpha values
    atoms_with_alpha = 0
    alpha_values = []
    alpha_by_type = {}
    
    for i in range(topology.get_num_atoms()):
        atom = topology.get_atom(i)
        
        # The atom should have alpha property or get_alpha() method
        if hasattr(atom, 'alpha'):
            alpha = atom.alpha
        elif hasattr(atom, 'get_alpha'):
            alpha = atom.get_alpha()
        else:
            # This is where current implementation fails
            continue
            
        if alpha > 0:
            atoms_with_alpha += 1
            alpha_values.append(alpha)
            atom_type = atom.type
            if atom_type not in alpha_by_type:
                alpha_by_type[atom_type] = []
            alpha_by_type[atom_type].append(alpha)
    
    # Verify total count
    assert atoms_with_alpha == 3814, f"Expected 3814 atoms with alpha values, got {atoms_with_alpha}"
    
    # Verify alpha value statistics
    assert len(alpha_values) == 3814
    avg_alpha = sum(alpha_values) / len(alpha_values)
    assert abs(avg_alpha - 0.8018) < 0.001, f"Expected average alpha 0.8018, got {avg_alpha:.4f}"
    
    # Verify range
    assert abs(min(alpha_values) - 0.1260) < 0.001, f"Expected min alpha 0.1260, got {min(alpha_values)}"
    assert abs(max(alpha_values) - 2.3990) < 0.001, f"Expected max alpha 2.3990, got {max(alpha_values)}"
    
    # Verify specific atom types have expected alpha values
    if "OE2" in alpha_by_type:
        avg_oe2 = sum(alpha_by_type["OE2"]) / len(alpha_by_type["OE2"])
        assert abs(avg_oe2 - 2.3990) < 0.001, f"Expected OE2 alpha 2.3990, got {avg_oe2}"


def test_drude_thole_parameters_parsing():
    """Test if PSF parser reads Thole screening parameters.
    
    Expected values from analyze_drude_psf.sh:
    - Parent atoms with Thole parameters: 3814
    - Average Thole parameter: -1.220
    - Thole parameter range: Min: -2.180, Max: -0.467
    """
    psf_path = os.path.join(TEST_DATA_DIR, "4wp7", "4wp7_drude.psf")
    
    if not os.path.exists(psf_path):
        pytest.skip(f"Test file {psf_path} not found")
    
    parser = PSFParser()
    topology = Topology()
    success = parser.parse_to_topology(psf_path, topology)
    assert success, "Failed to parse 4wp7_drude.psf"
    
    # Check if atoms have Thole values
    atoms_with_thole = 0
    thole_values = []
    
    for i in range(topology.get_num_atoms()):
        atom = topology.get_atom(i)
        
        # The atom should have thole property or get_thole() method
        # Note: The Atom class may need a thole field added
        if hasattr(atom, 'thole'):
            thole = atom.thole
        elif hasattr(atom, 'get_thole'):
            thole = atom.get_thole()
        else:
            # This is where current implementation fails
            continue
            
        if thole != 0:
            atoms_with_thole += 1
            thole_values.append(thole)
    
    # Verify count
    assert atoms_with_thole == 3814, f"Expected 3814 atoms with Thole values, got {atoms_with_thole}"
    
    # Verify Thole statistics
    avg_thole = sum(thole_values) / len(thole_values)
    assert abs(avg_thole - (-1.220)) < 0.001, f"Expected average Thole -1.220, got {avg_thole:.3f}"
    
    # Verify range
    assert abs(min(thole_values) - (-2.180)) < 0.001, f"Expected min Thole -2.180, got {min(thole_values)}"
    assert abs(max(thole_values) - (-0.467)) < 0.001, f"Expected max Thole -0.467, got {max(thole_values)}"


def test_drude_extended_atom_line_format():
    """Test the extended PSF atom line format by examining raw file."""
    psf_path = os.path.join(TEST_DATA_DIR, "4wp7", "4wp7_drude.psf")
    
    if not os.path.exists(psf_path):
        pytest.skip(f"Test file {psf_path} not found")
    
    # Find and examine a sample atom line with Drude parameters
    with open(psf_path, 'r') as f:
        in_atom_section = False
        sample_parent_line = None
        sample_drude_line = None
        
        for line in f:
            if "!NATOM" in line:
                in_atom_section = True
                continue
            elif "!NBOND" in line:
                break
            elif in_atom_section and line.strip():
                # Look for a parent atom with alpha (column 11 > 0)
                parts = line.split()
                if len(parts) >= 11:
                    atom_name = parts[4] if len(parts) > 4 else ""
                    atom_type = parts[5] if len(parts) > 5 else ""
                    
                    # Find a parent atom (not Drude, not lone pair)
                    if not atom_name.startswith("D") and not atom_name.startswith("LP"):
                        if len(parts) >= 11 and float(parts[10]) > 0:  # Has alpha
                            sample_parent_line = line.strip()
                    
                    # Find a Drude particle
                    if atom_type == "DRUD" and atom_name.startswith("D"):
                        sample_drude_line = line.strip()
                    
                    if sample_parent_line and sample_drude_line:
                        break
    
    assert sample_parent_line is not None, "Could not find parent atom with alpha in PSF"
    assert sample_drude_line is not None, "Could not find Drude particle in PSF"
    
    # Parse the parent atom line
    parts = sample_parent_line.split()
    assert len(parts) >= 11, f"Parent atom line has only {len(parts)} fields, expected at least 11 for extended format"
    
    # Extended format should have:
    # 1:atom_id 2:segname 3:resid 4:resname 5:atomname 6:atomtype 
    # 7:charge 8:mass 9:unused 10:thole 11:alpha
    alpha_value = float(parts[10])  # Column 11 (0-indexed: 10)
    thole_value = float(parts[9])   # Column 10 (0-indexed: 9)
    
    assert alpha_value > 0, f"Expected positive alpha value, got {alpha_value}"
    assert thole_value < 0, f"Expected negative Thole value, got {thole_value}"
    
    print(f"Sample parent atom line: {sample_parent_line}")
    print(f"  Alpha (polarizability): {alpha_value}")
    print(f"  Thole parameter: {thole_value}")
    
    # Parse the Drude particle line
    parts = sample_drude_line.split()
    drude_mass = float(parts[7])
    assert abs(drude_mass - 0.4) < 0.001, f"Drude particle mass should be 0.4, got {drude_mass}"
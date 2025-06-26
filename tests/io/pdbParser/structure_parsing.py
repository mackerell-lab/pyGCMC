# tests/io/pdbParser/structure_parsing.py
"""PDB Parser special structure parsing tests."""

import pytest
import os
import pygcmc

# Get the directory containing test data files
TEST_DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "data")


def test_parse_secondary_structure():
    """Test parsing HELIX and SHEET records."""
    pdb_path = os.path.join(TEST_DATA_DIR, "secondary.pdb")
    result = pygcmc.PDBParser.parse_file(pdb_path)
    
    # Check HELIX records
    assert len(result.helices) == 1  # Chain A has helices
    assert len(result.helices["A"]) == 2  # Two helices in chain A
    
    # Check helix classes
    helices = result.helices["A"]
    assert helices[0].helixClass == 1  # Right-handed alpha
    assert helices[1].helixClass == 1  # Right-handed alpha
    
    # Check SHEET records
    assert len(result.sheets) == 1  # Chain B has sheets
    assert len(result.sheets["B"]) == 2  # Two strands in chain B
    
    # Check sheet info format
    sheet_info = result.sheets["B"][1]  # Second strand
    parts = sheet_info.split(":")
    assert len(parts) == 3
    assert parts[0] == "S1"  # Sheet ID
    assert parts[1] == "2"   # Strand number
    assert parts[2] == "-1"  # Anti-parallel sense


def test_parse_ssbond():
    """Test parsing SSBOND records."""
    pdb_path = os.path.join(TEST_DATA_DIR, "ssbond.pdb")
    result = pygcmc.PDBParser.parse_file(pdb_path)

    # Check SSBOND records
    assert len(result.ssbonds) == 2

    # Check first bond (intra-chain)
    bond1 = result.ssbonds[0]
    parts = bond1.split("-")
    assert parts[0] == "A:3 "   # First CYS: chain A, residue 3
    assert parts[1] == "A:20 "  # Second CYS: chain A, residue 20

    # Check second bond (inter-chain)
    bond2 = result.ssbonds[1]
    parts = bond2.split("-")
    assert parts[0] == "A:15 "  # First CYS: chain A, residue 15
    assert parts[1] == "B:5 "   # Second CYS: chain B, residue 5

    # Check CYS residue structure
    cys_residues = [res for res in result.residues if res.get_resname() == "CYS"]
    assert len(cys_residues) > 0
    assert cys_residues[0].get_resname() == "CYS"



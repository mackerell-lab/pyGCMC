# tests/io/pdb/structure_parsing.py
"""PDB Parser special structure parsing tests."""

import pytest
import os
import pygcmc

# Get the directory containing test data files
TEST_DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "data")


def test_parse_secondary_structure():
    """Test parsing secondary structure information."""
    pdb_path = os.path.join(TEST_DATA_DIR, "secondary.pdb")
    result = pygcmc.PDBParser.parse_file(pdb_path)
    
    # Check that atoms were parsed
    assert len(result.atoms) > 0
    assert len(result.residues) > 0
    
    # Check for secondary structure elements
    # Note: This test depends on the specific content of secondary.pdb
    # We're checking for basic parsing capability rather than specific values
    helices = result.secondary_structures.get("HELIX", [])
    sheets = result.secondary_structures.get("SHEET", [])
    
    # At least one secondary structure element should be present
    assert len(helices) + len(sheets) > 0
    
    # Check helix properties if present
    if helices:
        helix = helices[0]
        assert hasattr(helix, 'initResName')
        assert hasattr(helix, 'initChainID')
        assert hasattr(helix, 'initSeqNum')
        assert hasattr(helix, 'endResName')
        assert hasattr(helix, 'endChainID')
        assert hasattr(helix, 'endSeqNum')


def test_parse_ssbond():
    """Test parsing disulfide bond information."""
    pdb_path = os.path.join(TEST_DATA_DIR, "ssbond.pdb")
    result = pygcmc.PDBParser.parse_file(pdb_path)
    
    # Check that atoms were parsed
    assert len(result.atoms) > 0
    assert len(result.residues) > 0
    
    # Check for disulfide bonds
    ssbonds = result.ssbonds
    assert len(ssbonds) > 0
    
    # Check disulfide bond properties
    ssbond = ssbonds[0]
    assert hasattr(ssbond, 'res1Name')
    assert hasattr(ssbond, 'chainID1') 
    assert hasattr(ssbond, 'seqNum1')
    assert hasattr(ssbond, 'res2Name')
    assert hasattr(ssbond, 'chainID2')
    assert hasattr(ssbond, 'seqNum2')
    
    # Both residues should be cysteine
    assert ssbond.res1Name == "CYS"
    assert ssbond.res2Name == "CYS"


def test_parse_crystal_info():
    """Test parsing crystallographic information."""
    pdb_path = os.path.join(TEST_DATA_DIR, "water.pdb")
    result = pygcmc.PDBParser.parse_file(pdb_path)
    
    # Check that crystal information was parsed if present
    if hasattr(result, 'crystal_info') and result.crystal_info:
        crystal = result.crystal_info
        
        # Check for unit cell parameters
        if hasattr(crystal, 'a'):
            assert crystal.a > 0
        if hasattr(crystal, 'b'):
            assert crystal.b > 0  
        if hasattr(crystal, 'c'):
            assert crystal.c > 0
        if hasattr(crystal, 'alpha'):
            assert 0 < crystal.alpha <= 180
        if hasattr(crystal, 'beta'):
            assert 0 < crystal.beta <= 180
        if hasattr(crystal, 'gamma'):
            assert 0 < crystal.gamma <= 180
    
    # At minimum, atoms should be parsed
    assert len(result.atoms) > 0
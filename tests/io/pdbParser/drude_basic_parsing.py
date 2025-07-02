#!/usr/bin/env python3
"""
Drude Force Field PDB Parser Tests - Basic Parsing

Tests for basic parsing and atom type recognition of Drude polarizable force field PDB files.
This test suite validates the parsing of 4wp7_drude.pdb which contains:
- Standard atoms (parent nuclear centers)  
- Drude particles (polarizable electron clouds, prefixed with 'D')
- Lone pairs (directional electron density, prefixed with 'LP')

Drude Force Field Background:
The Drude oscillator model treats polarization by attaching a virtual particle
(Drude particle) to each polarizable atom via a harmonic spring. This allows
explicit modeling of electronic polarization effects in molecular simulations.

Expected 4wp7_drude.pdb Structure (from raw file analysis):
- Total atoms: 12,979
- Parent atoms: ~7,630 (nuclear centers like N, CA, CB, etc.)  
- Drude particles: ~3,814 (virtual electrons like DN, DCA, DCB, etc.)
- Lone pairs: ~1,535 (electron density like LPOA, LPOB, etc.)
- Polarizable systems: 3,814 (parent-Drude pairs)
- Residue types: 20 amino acids with polarization
"""

import pytest
import os
from typing import Dict, List, Set
from collections import defaultdict

# Import pygcmc if available
try:
    import pygcmc
    from pygcmc.io import PDBParser
    PYGCMC_AVAILABLE = True
except ImportError:
    PYGCMC_AVAILABLE = False

# Test data directory
TEST_DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "data")


# Ground truth fixture is provided by conftest.py


@pytest.mark.skipif(not PYGCMC_AVAILABLE, reason="pygcmc not available")
def test_drude_pdb_file_exists():
    """Test that the Drude PDB file exists and is accessible"""
    pdb_path = os.path.join(TEST_DATA_DIR, "4wp7", "4wp7_drude.pdb")
    assert os.path.exists(pdb_path), f"Drude PDB file not found: {pdb_path}"
    assert os.path.getsize(pdb_path) > 0, "Drude PDB file is empty"


@pytest.mark.skipif(not PYGCMC_AVAILABLE, reason="pygcmc not available")
def test_drude_pdb_basic_parsing():
    """Test basic parsing of Drude PDB file"""
    pdb_path = os.path.join(TEST_DATA_DIR, "4wp7", "4wp7_drude.pdb")
    
    if not os.path.exists(pdb_path):
        pytest.skip(f"Drude PDB file not found: {pdb_path}")
    
    result = PDBParser.parse_file(pdb_path)
    assert result is not None, "Failed to parse Drude PDB file"
    assert len(result.atoms) > 0, "No atoms parsed from Drude PDB file"


@pytest.mark.skipif(not PYGCMC_AVAILABLE, reason="pygcmc not available") 
def test_drude_pdb_atom_counts(drude_ground_truth):
    """Test that parsed atom counts match raw file analysis"""
    if not drude_ground_truth:
        pytest.skip("Could not analyze raw Drude PDB file")
        
    pdb_path = os.path.join(TEST_DATA_DIR, "4wp7", "4wp7_drude.pdb")
    if not os.path.exists(pdb_path):
        pytest.skip(f"Drude PDB file not found: {pdb_path}")
    
    result = PDBParser.parse_file(pdb_path)
    
    # Verify total atom count matches raw file
    expected_total = drude_ground_truth['total_atoms']
    actual_total = len(result.atoms)
    assert actual_total == expected_total, \
        f"Total atom count mismatch: expected {expected_total}, got {actual_total}"


@pytest.mark.skipif(not PYGCMC_AVAILABLE, reason="pygcmc not available")
def test_drude_pdb_atom_type_recognition(drude_ground_truth):
    """Test recognition of different Drude atom types"""
    if not drude_ground_truth:
        pytest.skip("Could not analyze raw Drude PDB file")
        
    pdb_path = os.path.join(TEST_DATA_DIR, "4wp7", "4wp7_drude.pdb") 
    if not os.path.exists(pdb_path):
        pytest.skip(f"Drude PDB file not found: {pdb_path}")
    
    result = PDBParser.parse_file(pdb_path)
    
    # Categorize atoms by type
    parent_atoms = []
    drude_particles = []
    lone_pairs = []
    
    for atom in result.atoms:
        atom_name = atom.get_type()
        
        if atom_name.startswith('D') and len(atom_name) > 1:
            drude_particles.append(atom)
        elif atom_name.startswith('LP'):
            lone_pairs.append(atom)
        else:
            parent_atoms.append(atom)
    
    # Verify counts match ground truth (with tolerance for parser differences)
    expected_drude = drude_ground_truth['drude_particles']
    expected_lp = drude_ground_truth['lone_pairs'] 
    expected_parent = drude_ground_truth['parent_atoms']
    
    assert len(drude_particles) == expected_drude, \
        f"Drude particle count: expected {expected_drude}, got {len(drude_particles)}"
    
    assert len(lone_pairs) == expected_lp, \
        f"Lone pair count: expected {expected_lp}, got {len(lone_pairs)}"
    
    assert len(parent_atoms) == expected_parent, \
        f"Parent atom count: expected {expected_parent}, got {len(parent_atoms)}"


@pytest.mark.skipif(not PYGCMC_AVAILABLE, reason="pygcmc not available")
def test_drude_polarizable_residue_types(drude_ground_truth):
    """Test that all expected polarizable residue types are present"""
    if not drude_ground_truth:
        pytest.skip("Could not analyze raw Drude PDB file")
        
    pdb_path = os.path.join(TEST_DATA_DIR, "4wp7", "4wp7_drude.pdb")
    if not os.path.exists(pdb_path):
        pytest.skip(f"Drude PDB file not found: {pdb_path}")
    
    result = PDBParser.parse_file(pdb_path)
    
    # Find residue types with Drude particles
    polarizable_residues = set()
    for atom in result.atoms:
        if atom.get_type().startswith('D') and len(atom.get_type()) > 1:
            polarizable_residues.add(atom.get_resname())
    
    expected_polarizable = drude_ground_truth['polarizable_residues']
    
    # All expected polarizable residues should be found
    missing_residues = expected_polarizable - polarizable_residues
    assert not missing_residues, f"Missing polarizable residues: {missing_residues}"
    
    # Should have reasonable number of polarizable residue types (amino acids)
    assert len(polarizable_residues) >= 15, \
        f"Too few polarizable residue types: {len(polarizable_residues)}"



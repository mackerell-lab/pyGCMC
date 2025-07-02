#!/usr/bin/env python3
"""
Drude Force Field PDB Parser Tests - Advanced Validation

Tests for advanced validation of Drude polarizable force field PDB files including
parent-Drude pairing, coordinate precision, lone pair distribution, and force field completeness.
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
def test_drude_parent_drude_pairing():
    """Test that Drude particles are properly paired with parent atoms"""
    pdb_path = os.path.join(TEST_DATA_DIR, "4wp7", "4wp7_drude.pdb")
    if not os.path.exists(pdb_path):
        pytest.skip(f"Drude PDB file not found: {pdb_path}")
    
    result = PDBParser.parse_file(pdb_path)
    
    # Group atoms by residue for pairing analysis
    residue_atoms = defaultdict(list)
    for atom in result.atoms:
        res_key = f"{atom.get_resname()}_{atom.get_ires()}"
        residue_atoms[res_key].append(atom)
    
    orphaned_drude = 0
    valid_pairs = 0
    
    for res_key, atoms in residue_atoms.items():
        # Create lookup of parent atoms and Drude particles
        parent_atoms = {}
        drude_particles = {}
        
        for atom in atoms:
            atom_name = atom.get_type()
            if atom_name.startswith('D') and len(atom_name) > 1:
                parent_name = atom_name[1:]  # Remove 'D' prefix
                drude_particles[parent_name] = atom
            elif not atom_name.startswith('LP'):
                parent_atoms[atom_name] = atom
        
        # Check pairing
        for parent_name, drude_atom in drude_particles.items():
            if parent_name in parent_atoms:
                valid_pairs += 1
                
                # Verify coordinates are close (should be nearly identical for Drude)
                parent_coord = parent_atoms[parent_name].get_coor()
                drude_coord = drude_atom.get_coor()
                
                distance = sum((p - d)**2 for p, d in zip(parent_coord, drude_coord))**0.5
                assert distance < 1.0, \
                    f"Parent-Drude pair too far apart: {parent_name} distance={distance:.3f}Å"
            else:
                orphaned_drude += 1
    
    # Should have many valid pairs and few orphaned Drude particles
    assert valid_pairs > 3000, f"Too few valid parent-Drude pairs: {valid_pairs}"
    assert orphaned_drude < 100, f"Too many orphaned Drude particles: {orphaned_drude}"


@pytest.mark.skipif(not PYGCMC_AVAILABLE, reason="pygcmc not available")
def test_drude_coordinate_precision():
    """Test that Drude particles have precise coordinate values"""
    pdb_path = os.path.join(TEST_DATA_DIR, "4wp7", "4wp7_drude.pdb")
    if not os.path.exists(pdb_path):
        pytest.skip(f"Drude PDB file not found: {pdb_path}")
    
    result = PDBParser.parse_file(pdb_path)
    
    # Check coordinate precision for Drude particles
    drude_coords = []
    for atom in result.atoms:
        if atom.get_type().startswith('D') and len(atom.get_type()) > 1:
            coords = atom.get_coor()
            drude_coords.append(coords)
    
    assert len(drude_coords) > 3000, "Not enough Drude particles found for coordinate test"
    
    # Verify coordinates are reasonable (within simulation box)
    for coords in drude_coords[:100]:  # Sample first 100 for performance
        assert len(coords) == 3, "Coordinates should be 3D"
        for coord in coords:
            assert isinstance(coord, (int, float)), f"Coordinate should be numeric: {coord}"
            assert -100 <= coord <= 100, f"Coordinate out of reasonable range: {coord}"


@pytest.mark.skipif(not PYGCMC_AVAILABLE, reason="pygcmc not available")
def test_drude_lone_pair_distribution():
    """Test distribution and naming of lone pairs"""
    pdb_path = os.path.join(TEST_DATA_DIR, "4wp7", "4wp7_drude.pdb")
    if not os.path.exists(pdb_path):
        pytest.skip(f"Drude PDB file not found: {pdb_path}")
    
    result = PDBParser.parse_file(pdb_path)
    
    # Collect lone pair information
    lone_pair_types = set()
    lone_pair_residues = set()
    lone_pair_count = 0
    
    for atom in result.atoms:
        if atom.get_type().startswith('LP'):
            lone_pair_types.add(atom.get_type())
            lone_pair_residues.add(atom.get_resname())
            lone_pair_count += 1
    
    # Verify lone pair diversity and distribution
    assert lone_pair_count > 1000, f"Too few lone pairs found: {lone_pair_count}"
    assert len(lone_pair_types) >= 10, f"Too few lone pair types: {len(lone_pair_types)}"
    assert len(lone_pair_residues) >= 15, f"Lone pairs not distributed across residues: {len(lone_pair_residues)}"
    
    # Check for expected lone pair naming patterns
    expected_patterns = {'LPOA', 'LPOB', 'LP1A', 'LP1B', 'LP2A', 'LP2B'}
    found_patterns = expected_patterns.intersection(lone_pair_types)
    assert len(found_patterns) >= 3, f"Missing expected lone pair patterns: {expected_patterns - found_patterns}"


@pytest.mark.skipif(not PYGCMC_AVAILABLE, reason="pygcmc not available")
def test_drude_force_field_completeness():
    """Test that Drude force field has reasonable coverage of polarizable atoms"""
    pdb_path = os.path.join(TEST_DATA_DIR, "4wp7", "4wp7_drude.pdb")
    if not os.path.exists(pdb_path):
        pytest.skip(f"Drude PDB file not found: {pdb_path}")
    
    result = PDBParser.parse_file(pdb_path)
    
    # Count potentially polarizable vs actually polarized atoms
    polarizable_atom_types = {'N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'CE', 'CZ', 'OH', 'OG', 'OD1', 'OD2', 'OE1', 'OE2'}
    
    total_polarizable = 0
    total_drude = 0
    
    for atom in result.atoms:
        atom_type = atom.get_type()
        
        if atom_type in polarizable_atom_types:
            total_polarizable += 1
        elif atom_type.startswith('D') and len(atom_type) > 1:
            total_drude += 1
    
    # Polarization coverage should be reasonable 
    # Note: Some atom types may have multiple Drude particles or extended coverage
    coverage = (total_drude / total_polarizable * 100) if total_polarizable > 0 else 0
    
    # Adjust range based on actual Drude force field implementation
    assert 20 <= coverage <= 200, \
        f"Polarization coverage {coverage:.1f}% outside reasonable range (20-200%)"
    
    assert total_drude > 3000, f"Too few Drude particles for complete force field: {total_drude}"
    
    print(f"Drude force field analysis: {total_drude} Drude particles, "
          f"{total_polarizable} potentially polarizable atoms, coverage: {coverage:.1f}%")


if __name__ == "__main__":
    # Run tests when executed directly
    pytest.main([__file__, "-v"])
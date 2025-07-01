# tests/io/psfParser/drude_particles.py
"""
Drude PSF Parser Tests - Drude Particles and Lone Pairs

This file contains tests for Drude particle and lone pair parsing.
Based on analysis of 4wp7_drude.psf file.

Tests include:
- Basic Drude PSF parsing
- Drude particle identification and properties
- Lone pair parsing and classification
- Charge distribution analysis
"""

import pytest
import os
from pygcmc.io import PSFParser
from pygcmc.model import Topology

# Get the directory containing test data files
TEST_DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "data")


def test_parse_drude_psf_basic():
    """Test basic parsing of Drude PSF file with all sections."""
    psf_path = os.path.join(TEST_DATA_DIR, "4wp7", "4wp7_drude.psf")
    
    # Skip test if file doesn't exist
    if not os.path.exists(psf_path):
        pytest.skip(f"Test file {psf_path} not found")
    
    parser = PSFParser()
    topology = Topology()
    
    # Parse the PSF file
    success = parser.parse_to_topology(psf_path, topology)
    assert success, "Failed to parse 4wp7_drude.psf"
    
    # Verify basic counts from analysis
    assert topology.get_num_atoms() == 12979, f"Expected 12979 atoms, got {topology.get_num_atoms()}"
    assert topology.get_num_bonds() == 13008, f"Expected 13008 bonds, got {topology.get_num_bonds()}"
    assert topology.get_num_angles() == 13981, f"Expected 13981 angles, got {topology.get_num_angles()}"
    assert topology.get_num_dihedrals() == 20482, f"Expected 20482 dihedrals, got {topology.get_num_dihedrals()}"
    assert topology.get_num_impropers() == 1160, f"Expected 1160 impropers, got {topology.get_num_impropers()}"
    
    # Verify CMAP cross-terms
    assert topology.get_num_cmaps() == 492, f"Expected 492 CMAP terms, got {topology.get_num_cmaps()}"


def test_parse_drude_particles():
    """Test parsing and validation of Drude particles."""
    psf_path = os.path.join(TEST_DATA_DIR, "4wp7", "4wp7_drude.psf")
    
    if not os.path.exists(psf_path):
        pytest.skip(f"Test file {psf_path} not found")
    
    parser = PSFParser()
    topology = Topology()
    success = parser.parse_to_topology(psf_path, topology)
    assert success, "Failed to parse 4wp7_drude.psf"
    
    # Count Drude particles and verify properties
    drude_count = 0
    drude_masses = []
    drude_charges = []
    parent_atoms = 0
    
    for i in range(topology.get_num_atoms()):
        atom = topology.get_atom(i)
        atom_name = atom.name
        atom_type = atom.type
        
        # Drude particles have type "DRUD" and names starting with "D"
        if atom_type == "DRUD" and atom_name.startswith("D"):
            drude_count += 1
            drude_masses.append(atom.mass)
            drude_charges.append(atom.charge)
        elif not atom_name.startswith("LP"):  # Not a lone pair
            parent_atoms += 1
    
    # Verify counts match analysis
    assert drude_count == 3814, f"Expected 3814 Drude particles, got {drude_count}"
    assert parent_atoms == 7630, f"Expected 7630 parent atoms, got {parent_atoms}"
    
    # Verify Drude particle masses (should all be 0.4 amu)
    assert all(abs(mass - 0.4) < 1e-4 for mass in drude_masses), \
        "All Drude particles should have mass 0.4 amu"
    
    # Verify total Drude charge
    total_drude_charge = sum(drude_charges)
    assert abs(total_drude_charge - (-7143.307)) < 0.1, \
        f"Expected total Drude charge -7143.307, got {total_drude_charge:.3f}"
    
    # Verify polarization coverage
    polarization_ratio = drude_count / parent_atoms
    assert 0.499 <= polarization_ratio <= 0.500, \
        f"Expected ~50% polarization coverage, got {polarization_ratio:.3f}"


def test_parse_lone_pairs():
    """Test parsing of lone pairs in Drude PSF."""
    psf_path = os.path.join(TEST_DATA_DIR, "4wp7", "4wp7_drude.psf")
    
    if not os.path.exists(psf_path):
        pytest.skip(f"Test file {psf_path} not found")
    
    parser = PSFParser()
    topology = Topology()
    success = parser.parse_to_topology(psf_path, topology)
    assert success, "Failed to parse 4wp7_drude.psf"
    
    # Count lone pairs by type
    lone_pair_types = {}
    total_lone_pairs = 0
    
    for i in range(topology.get_num_atoms()):
        atom = topology.get_atom(i)
        if atom.name.startswith("LP"):
            total_lone_pairs += 1
            lp_type = atom.name
            lone_pair_types[lp_type] = lone_pair_types.get(lp_type, 0) + 1
    
    # Verify total lone pair count
    assert total_lone_pairs == 1535, f"Expected 1535 lone pairs, got {total_lone_pairs}"
    
    # Verify specific lone pair types from analysis
    expected_lp_counts = {
        "LPOB": 493,
        "LPOA": 493,
        "LP1B": 82,
        "LP1A": 82,
        "LPGB": 65,
        "LPGA": 65,
        "LP2B": 58,
        "LP2A": 58,
        "LP": 45
    }
    
    for lp_type, expected_count in expected_lp_counts.items():
        actual_count = lone_pair_types.get(lp_type, 0)
        assert actual_count == expected_count, \
            f"Expected {expected_count} {lp_type} lone pairs, got {actual_count}"


def test_parse_drude_charge_neutrality():
    """Test charge neutrality in Drude PSF system."""
    psf_path = os.path.join(TEST_DATA_DIR, "4wp7", "4wp7_drude.psf")
    
    if not os.path.exists(psf_path):
        pytest.skip(f"Test file {psf_path} not found")
    
    parser = PSFParser()
    topology = Topology()
    success = parser.parse_to_topology(psf_path, topology)
    assert success, "Failed to parse 4wp7_drude.psf"
    
    # Calculate charges by category
    parent_charge = 0.0
    drude_charge = 0.0
    lone_pair_charge = 0.0
    total_charge = 0.0
    
    for i in range(topology.get_num_atoms()):
        atom = topology.get_atom(i)
        charge = atom.charge
        total_charge += charge
        
        if atom.type == "DRUD" and atom.name.startswith("D"):
            drude_charge += charge
        elif atom.name.startswith("LP"):
            lone_pair_charge += charge
        else:
            parent_charge += charge
    
    # Verify charge values match analysis
    assert abs(total_charge - (-3.0)) < 0.1, \
        f"Expected total charge -3.0, got {total_charge:.3f}"
    assert abs(parent_charge - 7561.173) < 0.1, \
        f"Expected parent charge 7561.173, got {parent_charge:.3f}"
    assert abs(drude_charge - (-7143.307)) < 0.1, \
        f"Expected Drude charge -7143.307, got {drude_charge:.3f}"
    assert abs(lone_pair_charge - (-420.866)) < 0.1, \
        f"Expected lone pair charge -420.866, got {lone_pair_charge:.3f}"
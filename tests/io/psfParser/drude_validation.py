# tests/io/psfParser/drude_validation.py
"""
Drude PSF Parser Tests - System Validation

This file contains validation tests for Drude PSF systems.
Based on analysis of 4wp7_drude.psf file.

Tests include:
- Connectivity validation
- Hydrogen bonding analysis
- Atom type distribution
- Residue composition validation
"""

import pytest
import os
from pygcmc.io import PSFParser
from pygcmc.model import Topology

# Get the directory containing test data files
TEST_DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "data")


def test_parse_drude_connectivity():
    """Test connectivity statistics in Drude PSF."""
    psf_path = os.path.join(TEST_DATA_DIR, "4wp7", "4wp7_drude.psf")
    
    if not os.path.exists(psf_path):
        pytest.skip(f"Test file {psf_path} not found")
    
    parser = PSFParser()
    topology = Topology()
    success = parser.parse_to_topology(psf_path, topology)
    assert success, "Failed to parse 4wp7_drude.psf"
    
    num_atoms = topology.get_num_atoms()
    num_bonds = topology.get_num_bonds()
    num_angles = topology.get_num_angles()
    num_dihedrals = topology.get_num_dihedrals()
    num_impropers = topology.get_num_impropers()
    
    # Calculate connectivity metrics
    bond_density = num_bonds / num_atoms
    angle_density = num_angles / num_atoms
    dihedral_density = num_dihedrals / num_atoms
    improper_density = num_impropers / num_atoms
    
    # Verify densities match analysis
    assert abs(bond_density - 1.002) < 0.01, \
        f"Expected bond density 1.002, got {bond_density:.3f}"
    assert abs(angle_density - 1.077) < 0.01, \
        f"Expected angle density 1.077, got {angle_density:.3f}"
    assert abs(dihedral_density - 1.578) < 0.01, \
        f"Expected dihedral density 1.578, got {dihedral_density:.3f}"
    assert abs(improper_density - 0.089) < 0.01, \
        f"Expected improper density 0.089, got {improper_density:.3f}"
    
    # Verify topology richness
    topology_richness = (num_bonds + num_angles + num_dihedrals) / num_atoms
    assert abs(topology_richness - 3.657) < 0.01, \
        f"Expected topology richness 3.657, got {topology_richness:.3f}"


def test_parse_drude_hydrogen_bonding():
    """Test hydrogen bonding information in Drude PSF."""
    psf_path = os.path.join(TEST_DATA_DIR, "4wp7", "4wp7_drude.psf")
    
    if not os.path.exists(psf_path):
        pytest.skip(f"Test file {psf_path} not found")
    
    parser = PSFParser()
    topology = Topology()
    success = parser.parse_to_topology(psf_path, topology)
    assert success, "Failed to parse 4wp7_drude.psf"
    
    # Get donor and acceptor information
    num_donors = topology.get_num_donors()
    num_acceptors = topology.get_num_acceptors()
    total_hbond_sites = num_donors + num_acceptors
    hbond_density = total_hbond_sites / topology.get_num_atoms()
    
    # Verify hydrogen bonding statistics
    assert num_donors == 823, f"Expected 823 donors, got {num_donors}"
    assert num_acceptors == 727, f"Expected 727 acceptors, got {num_acceptors}"
    assert total_hbond_sites == 1550, f"Expected 1550 H-bond sites, got {total_hbond_sites}"
    assert abs(hbond_density - 0.1194) < 0.001, \
        f"Expected H-bond density 0.1194, got {hbond_density:.4f}"


def test_parse_drude_atom_types():
    """Test atom type distribution in Drude PSF."""
    psf_path = os.path.join(TEST_DATA_DIR, "4wp7", "4wp7_drude.psf")
    
    if not os.path.exists(psf_path):
        pytest.skip(f"Test file {psf_path} not found")
    
    parser = PSFParser()
    topology = Topology()
    success = parser.parse_to_topology(psf_path, topology)
    assert success, "Failed to parse 4wp7_drude.psf"
    
    # Count atoms by element category
    element_counts = {
        "Carbon": 0,
        "Hydrogen": 0,
        "Oxygen": 0,
        "Nitrogen": 0,
        "Sulfur": 0,
        "Drude": 0,
        "LonePair": 0,
        "Other": 0
    }
    
    for i in range(topology.get_num_atoms()):
        atom = topology.get_atom(i)
        atom_name = atom.name
        atom_type = atom.type
        
        if atom_type == "DRUD" and atom_name.startswith("D"):
            element_counts["Drude"] += 1
        elif atom_name.startswith("LP"):
            element_counts["LonePair"] += 1
        elif atom_name.startswith("H"):
            element_counts["Hydrogen"] += 1
        elif atom_name.startswith("C"):
            element_counts["Carbon"] += 1
        elif atom_name.startswith("N"):
            element_counts["Nitrogen"] += 1
        elif atom_name.startswith("O"):
            element_counts["Oxygen"] += 1
        elif atom_name.startswith("S"):
            element_counts["Sulfur"] += 1
        else:
            element_counts["Other"] += 1
    
    # Verify element counts match analysis
    assert element_counts["Drude"] == 3814, f"Expected 3814 Drude atoms, got {element_counts['Drude']}"
    assert element_counts["LonePair"] == 1535, f"Expected 1535 lone pairs, got {element_counts['LonePair']}"
    assert element_counts["Hydrogen"] == 3816, f"Expected 3816 hydrogen atoms, got {element_counts['Hydrogen']}"
    assert element_counts["Carbon"] == 2434, f"Expected 2434 carbon atoms, got {element_counts['Carbon']}"
    assert element_counts["Oxygen"] == 718, f"Expected 718 oxygen atoms, got {element_counts['Oxygen']}"
    assert element_counts["Nitrogen"] == 643, f"Expected 643 nitrogen atoms, got {element_counts['Nitrogen']}"
    assert element_counts["Sulfur"] == 19, f"Expected 19 sulfur atoms, got {element_counts['Sulfur']}"


def test_parse_drude_residue_composition():
    """Test residue composition in Drude PSF."""
    psf_path = os.path.join(TEST_DATA_DIR, "4wp7", "4wp7_drude.psf")
    
    if not os.path.exists(psf_path):
        pytest.skip(f"Test file {psf_path} not found")
    
    parser = PSFParser()
    topology = Topology()
    success = parser.parse_to_topology(psf_path, topology)
    assert success, "Failed to parse 4wp7_drude.psf"
    
    # Count residues by type
    residue_counts = {}
    
    for i in range(topology.get_num_residues()):
        residue = topology.get_residue(i)
        res_name = residue.name
        residue_counts[res_name] = residue_counts.get(res_name, 0) + 1
    
    # Verify we have protein residues
    protein_residues = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", 
                       "GLY", "HIS", "HSD", "ILE", "LEU", "LYS", "MET", 
                       "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]
    
    total_protein_residues = sum(residue_counts.get(res, 0) for res in protein_residues)
    assert total_protein_residues > 0, "No protein residues found in Drude PSF"
    
    # Verify we have the expected number of unique residues
    # From analysis, we know there are 494 residues with proper Drude pairing
    assert total_protein_residues >= 494, \
        f"Expected at least 494 protein residues, got {total_protein_residues}"
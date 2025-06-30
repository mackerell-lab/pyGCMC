# tests/io/topParser/gcmc_systems.py
"""TOP Parser tests for GCMC systems - focusing on 4wp7 topology parsing."""

import pytest
import os
import pygcmc
import re
from pygcmc.io import TOPParser
from pygcmc.model import Topology

# Get the directory containing test data files
TEST_DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "data")


def test_parse_4wp7_gcmc_topology():
    """Test parsing 4wp7 GCMC system TOP file with multiple protein chains."""
    top_path = os.path.join(TEST_DATA_DIR, "4wp7", "4wp7_fixed_with_5l13_silcs.1.gc.74.top")
    
    # Skip test if file doesn't exist
    if not os.path.exists(top_path):
        pytest.skip(f"Test file {top_path} not found")
    
    parser = TOPParser()
    topology = Topology()
    
    # Parse the topology file
    success = parser.parse_to_topology(top_path, topology)
    assert success, "Failed to parse 4wp7 TOP file"
    
    # Verify basic topology information
    assert topology.get_num_atoms() > 20000, f"Expected >20k atoms, got {topology.get_num_atoms()}"
    # Add stricter check for this specific large system
    assert topology.get_num_atoms() > 150000, f"Expected >150k atoms for 4wp7 system, got {topology.get_num_atoms()}"
    
    # Should have multiple residues due to multiple protein chains
    assert topology.get_num_residues() > 100, f"Expected >100 residues, got {topology.get_num_residues()}"
    # Add stricter check for large GCMC system
    assert topology.get_num_residues() > 5000, f"Expected >5k residues for 4wp7 GCMC system, got {topology.get_num_residues()}"
    
    # Multiple protein chains means multiple segments
    assert topology.get_num_segments() >= 1, f"Expected multiple segments, got {topology.get_num_segments()}"
    # Add stricter check for truly multiple segments
    assert topology.get_num_segments() > 1, f"Expected multiple segments for multi-chain system, got {topology.get_num_segments()}"
    
    # Should have substantial bonding information
    assert topology.get_num_bonds() > 10000, f"Expected >10k bonds, got {topology.get_num_bonds()}"
    assert topology.get_num_angles() > 15000, f"Expected >15k angles, got {topology.get_num_angles()}"
    assert topology.get_num_dihedrals() > 20000, f"Expected >20k dihedrals, got {topology.get_num_dihedrals()}"
    # Add stricter checks for this large system
    assert topology.get_num_bonds() > 100000, f"Expected >100k bonds for 4wp7 system, got {topology.get_num_bonds()}"
    assert topology.get_num_angles() > 200000, f"Expected >200k angles for 4wp7 system, got {topology.get_num_angles()}"
    assert topology.get_num_dihedrals() > 300000, f"Expected >300k dihedrals for 4wp7 system, got {topology.get_num_dihedrals()}"
    
    # Find ASP residue (should be residue 8 based on TOP file)
    asp_residues = []
    for i in range(topology.get_num_residues()):
        residue = topology.get_residue(i)
        if residue.name == "ASP" and residue.resid == 8:
            asp_residues.append(residue)
    
    assert len(asp_residues) >= 1, "Should find at least one ASP residue with resid 8"
    
    # Check first ASP residue atoms
    asp = asp_residues[0]
    asp_atom_names = []
    for atom_idx in asp.atoms:
        atom = topology.get_atom(atom_idx)
        asp_atom_names.append(atom.name)
    
    # ASP should have standard amino acid atoms
    expected_asp_atoms = {"N", "CA", "C", "O", "CB", "CG", "OD1", "OD2"}
    found_atoms = set(asp_atom_names)
    common_atoms = expected_asp_atoms.intersection(found_atoms)
    
    assert len(common_atoms) >= 6, f"ASP should have standard atoms, found {found_atoms}"


def test_4wp7_topology_atom_types():
    """Test atom types and parameters in 4wp7 topology.""" 
    top_path = os.path.join(TEST_DATA_DIR, "4wp7", "4wp7_fixed_with_5l13_silcs.1.gc.74.top")
    
    # Skip test if file doesn't exist
    if not os.path.exists(top_path):
        pytest.skip(f"Test file {top_path} not found")
    
    parser = TOPParser()
    topology = Topology()
    
    # Parse the topology file
    success = parser.parse_to_topology(top_path, topology)
    assert success, "Failed to parse 4wp7 TOP file"
    
    # Check that we have expected CHARMM36 atom types
    atom_types_found = set()
    charges_found = []
    masses_found = []
    
    for i in range(topology.get_num_atoms()):
        atom = topology.get_atom(i)
        atom_types_found.add(atom.type)
        charges_found.append(atom.charge)
        masses_found.append(atom.mass)
    
    # Should find standard CHARMM36 protein atom types
    expected_charmm_types = {"NH3", "HC", "CT1", "HB", "C", "O", "CT2", "CT3"}
    found_charmm_types = atom_types_found.intersection(expected_charmm_types)
    
    assert len(found_charmm_types) >= 4, f"Expected CHARMM36 types, found {found_charmm_types}"
    
    # Verify charges are reasonable (should be between -2 and +2 for most atoms)
    reasonable_charges = [c for c in charges_found if -2.0 <= c <= 2.0]
    assert len(reasonable_charges) > len(charges_found) * 0.95, "Most charges should be reasonable"
    
    # Verify masses are reasonable (should be > 0 and < 200 for standard atoms)
    reasonable_masses = [m for m in masses_found if 0.5 <= m <= 200.0]
    assert len(reasonable_masses) > len(masses_found) * 0.9, "Most masses should be reasonable"
    
    # Check specific N-terminal nitrogen (should be NH3 type with positive charge)
    first_n_atom = None
    for i in range(min(100, topology.get_num_atoms())):  # Check first 100 atoms
        atom = topology.get_atom(i)
        if atom.name == "N" and atom.type == "NH3":
            first_n_atom = atom
            break
    
    if first_n_atom:
        # N-terminal nitrogen should have negative charge (around -0.3)
        assert -0.5 <= first_n_atom.charge <= 0.0, f"N-terminal N charge should be negative, got {first_n_atom.charge}"
        # Standard nitrogen mass
        assert 13.0 <= first_n_atom.mass <= 15.0, f"Nitrogen mass should be ~14, got {first_n_atom.mass}"


def test_4wp7_force_field_includes():
    """Test force field include directives in 4wp7 topology."""
    top_path = os.path.join(TEST_DATA_DIR, "4wp7", "4wp7_fixed_with_5l13_silcs.1.gc.74.top")
    
    # Skip test if file doesn't exist
    if not os.path.exists(top_path):
        pytest.skip(f"Test file {top_path} not found")
    
    # Read file content to verify includes
    with open(top_path, 'r') as f:
        content = f.read()
    
    # Should include CHARMM36 force field
    assert "charmm36" in content.lower(), "Should reference CHARMM36 force field"
    assert "#include" in content, "Should have include directives"
    
    # Check for proper moleculetype definition
    assert "[moleculetype]" in content, "Should have moleculetype section"
    assert "Protein_chain_P" in content, "Should define protein chain molecule type"
    
    # Should have atoms section
    assert "[atoms]" in content, "Should have atoms section"
    
    # Verify file structure includes key sections
    required_sections = ["moleculetype", "atoms"]
    for section in required_sections:
        pattern = rf"\[\s*{section}\s*\]"
        assert re.search(pattern, content), f"Should have [{section}] section"
    
    # Check that includes reference proper file structure
    include_pattern = r"#include\s+[\"']([^\"']+)[\"']"
    includes = re.findall(include_pattern, content)
    
    assert len(includes) >= 1, "Should have at least one include directive"
    
    # First include should be the force field
    if includes:
        first_include = includes[0]
        assert "forcefield" in first_include.lower() or "charmm" in first_include.lower(), \
            f"First include should be force field file, got {first_include}"


def test_4wp7_multiple_protein_chains():
    """Test that 4wp7 topology correctly handles multiple identical protein chains."""
    top_path = os.path.join(TEST_DATA_DIR, "4wp7", "4wp7_fixed_with_5l13_silcs.1.gc.74.top")
    
    # Skip test if file doesn't exist
    if not os.path.exists(top_path):
        pytest.skip(f"Test file {top_path} not found")
    
    parser = TOPParser()
    topology = Topology()
    
    # Parse the topology file
    success = parser.parse_to_topology(top_path, topology)
    assert success, "Failed to parse 4wp7 TOP file"
    
    # Count atoms by residue type to detect multiple chains
    residue_counts = {}
    for i in range(topology.get_num_residues()):
        residue = topology.get_residue(i)
        key = (residue.name, residue.resid)
        residue_counts[key] = residue_counts.get(key, 0) + 1
    
    # Should have multiple copies of the same residues (indicating multiple protein chains)
    asp_8_count = residue_counts.get(("ASP", 8), 0)
    
    # From the TOP file analysis, we saw multiple ASP 8 entries, indicating multiple chains
    assert asp_8_count >= 1, f"Should have at least one ASP residue 8, found {asp_8_count}"
    # Add stricter check for multiple chains
    assert asp_8_count >= 2, f"Should have multiple ASP residue 8 for multi-chain system, found {asp_8_count}"
    
    # Count total protein-like residues
    protein_residue_names = {"ALA", "ARG", "ASN", "ASP", "GLN", "GLU", "GLY", "HIS", "ILE",
                            "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"}
    
    protein_residues = 0
    for i in range(topology.get_num_residues()):
        residue = topology.get_residue(i)
        if residue.name in protein_residue_names:
            protein_residues += 1
    
    # Should have substantial number of protein residues due to multiple chains
    assert protein_residues > 50, f"Expected many protein residues from multiple chains, got {protein_residues}"
    
    # Verify we have reasonable atom count for multiple protein chains
    total_atoms = topology.get_num_atoms()
    atoms_per_chain_estimate = total_atoms // max(1, asp_8_count) if asp_8_count > 0 else total_atoms
    
    # Each protein chain should have reasonable number of atoms (few hundred to few thousand)
    if asp_8_count > 0:
        assert 100 <= atoms_per_chain_estimate <= 10000, \
            f"Atoms per chain estimate ({atoms_per_chain_estimate}) seems unreasonable"
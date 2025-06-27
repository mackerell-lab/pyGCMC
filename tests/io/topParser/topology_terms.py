# tests/io/topParser/topology_terms.py
"""TOP Parser topology terms tests."""

import os
import pytest
from pygcmc.io import TOPParser
from pygcmc.model import Topology, TopologyResidue, TopologyAtom

# Get the directory containing test data files
TEST_DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "data")


def test_parse_bonds():
    """Test parsing bonds section from topology file."""
    top_file = os.path.join(TEST_DATA_DIR, "test.top")
    parser = TOPParser()
    topology = Topology()
    
    assert parser.parse_to_topology(top_file, topology), "Failed to parse topology file"
    
    # Check total number of bonds
    assert topology.get_num_bonds() == 163, "Wrong number of bonds (should be 163: protein + BENX + 2×PRPX)"
    
    # Check specific bonds in ALA-7 (N-terminal)
    ala_id = topology.find_residue("ALA", 7)
    ala = topology.get_residue(ala_id)
    
    # Find atom indices for N-terminal ALA
    n_idx = None
    ht_indices = []
    ca_idx = None
    
    for atom_idx in ala.atoms:
        atom = topology.get_atom(atom_idx)
        if atom.name == "N":
            n_idx = atom_idx
        elif atom.name in ["HT1", "HT2", "HT3"]:
            ht_indices.append(atom_idx)
        elif atom.name == "CA":
            ca_idx = atom_idx
    
    # Check N-H bonds
    for ht_idx in ht_indices:
        assert topology.has_bond(n_idx, ht_idx), f"Missing N-HT bond between atoms {n_idx} and {ht_idx}"
    
    # Check N-CA bond
    assert topology.has_bond(n_idx, ca_idx), f"Missing N-CA bond between atoms {n_idx} and {ca_idx}"


def test_parse_angles():
    """Test parsing angles section from topology file."""
    top_file = os.path.join(TEST_DATA_DIR, "test.top")
    parser = TOPParser()
    topology = Topology()
    
    assert parser.parse_to_topology(top_file, topology), "Failed to parse topology file"
    
    # Check total number of angles
    assert topology.get_num_angles() == 297, "Wrong number of angles (should be 297: protein[243] + BENX + 2×PRPX[18×2])"
    
    # Check specific angles in ALA-7 (N-terminal)
    ala_id = topology.find_residue("ALA", 7)
    ala = topology.get_residue(ala_id)
    
    # Find atom indices for N-terminal ALA
    n_idx = None
    ht_indices = []
    ca_idx = None
    
    for atom_idx in ala.atoms:
        atom = topology.get_atom(atom_idx)
        if atom.name == "N":
            n_idx = atom_idx
        elif atom.name in ["HT1", "HT2", "HT3"]:
            ht_indices.append(atom_idx)
        elif atom.name == "CA":
            ca_idx = atom_idx
    
    # Check HT-N-HT angles
    for i, ht1_idx in enumerate(ht_indices):
        for ht2_idx in ht_indices[i+1:]:
            assert topology.has_angle(ht1_idx, n_idx, ht2_idx), \
                f"Missing HT-N-HT angle between atoms {ht1_idx}, {n_idx}, and {ht2_idx}"
    
    # Check HT-N-CA angles
    for ht_idx in ht_indices:
        assert topology.has_angle(ht_idx, n_idx, ca_idx), \
            f"Missing HT-N-CA angle between atoms {ht_idx}, {n_idx}, and {ca_idx}"


def test_parse_dihedrals():
    """Test parsing dihedrals section from topology file."""
    top_file = os.path.join(TEST_DATA_DIR, "test.top")
    parser = TOPParser()
    topology = Topology()
    
    assert parser.parse_to_topology(top_file, topology), "Failed to parse topology file"
    
    # Check total number of dihedrals
    assert topology.get_num_dihedrals() == 422, "Wrong number of dihedrals (should be 422: protein[362] + BENX[24] + 2×PRPX[18×2])"
    
    # Check specific dihedrals in ALA-7 (N-terminal)
    ala_id = topology.find_residue("ALA", 7)
    ala = topology.get_residue(ala_id)
    
    # Find atom indices for N-terminal ALA
    n_idx = None
    ca_idx = None
    c_idx = None
    ha_idx = None  # Alpha hydrogen on CA
    
    # Print all atom names to help debug
    print("\nAtoms in N-terminal ALA:")
    for atom_idx in ala.atoms:
        atom = topology.get_atom(atom_idx)
        print(f"  {atom.name} (type: {atom.type})")
        if atom.name == "N":
            n_idx = atom_idx
        elif atom.name == "CA":
            ca_idx = atom_idx
        elif atom.name == "C":
            c_idx = atom_idx
        elif atom.name == "HA":  # Alpha hydrogen commonly exists in CHARMM
            ha_idx = atom_idx
    
    # Check if we found the backbone atoms
    assert n_idx is not None, "N atom not found in N-terminal ALA"
    assert ca_idx is not None, "CA atom not found in N-terminal ALA"
    assert c_idx is not None, "C atom not found in N-terminal ALA"
    assert ha_idx is not None, "HA atom not found in N-terminal ALA"
    
    # Check backbone dihedral that should exist in CHARMM
    # N-CA-C-O or N-CA-C-next_N are common backbone dihedrals
    found_backbone_dihedral = False
    
    # Use range() to iterate over atom indices
    for atom_idx in range(topology.get_num_atoms()):
        # Look for O or next residue's N that forms dihedral with N-CA-C
        atom = topology.get_atom(atom_idx)
        if atom_idx not in ala.atoms and (atom.name == "O" or atom.name == "N"):
            if topology.has_dihedral(n_idx, ca_idx, c_idx, atom_idx):
                found_backbone_dihedral = True
                print(f"Found backbone dihedral: N-CA-C-{atom.name}")
                break
    
    assert found_backbone_dihedral, "No backbone dihedral found for N-terminal ALA"

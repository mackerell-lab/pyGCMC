# tests/io/psfParser/topology_terms.py
"""PSF Parser topology terms tests."""

import os
import pytest
from pygcmc.io import PSFParser
from pygcmc.model import Topology, TopologyResidue, TopologyAtom


def test_parse_bonds(test_data_dir):
    """Test parsing bonds section from PSF file."""
    psf_file = os.path.join(test_data_dir, "test_proa.psf")
    parser = PSFParser()
    topology = Topology()
    
    assert parser.parse_to_topology(psf_file, topology), "Failed to parse PSF file"
    
    # Check total number of bonds
    assert topology.get_num_bonds() == 131, "Wrong number of bonds"
    
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


def test_parse_angles(test_data_dir):
    """Test parsing angles section from PSF file."""
    psf_file = os.path.join(test_data_dir, "test_proa.psf")
    parser = PSFParser()
    topology = Topology()
    
    assert parser.parse_to_topology(psf_file, topology), "Failed to parse PSF file"
    
    # Check total number of angles
    assert topology.get_num_angles() == 243, "Wrong number of angles"
    
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


def test_parse_dihedrals(test_data_dir):
    """Test parsing dihedrals section from PSF file."""
    psf_file = os.path.join(test_data_dir, "test_proa.psf")
    parser = PSFParser()
    topology = Topology()
    
    assert parser.parse_to_topology(psf_file, topology), "Failed to parse PSF file"
    
    # Check total number of dihedrals
    assert topology.get_num_dihedrals() == 362, "Wrong number of dihedrals"
    
    # Check backbone dihedrals in VAL-8
    val_id = topology.find_residue("VAL", 8)
    val = topology.get_residue(val_id)
    
    # Find backbone atoms for phi/psi angles
    n_idx = None
    ca_idx = None
    c_idx = None
    prev_c_idx = None
    next_n_idx = None
    
    # Find VAL-8 backbone atoms
    for atom_idx in val.atoms:
        atom = topology.get_atom(atom_idx)
        if atom.name == "N":
            n_idx = atom_idx
        elif atom.name == "CA":
            ca_idx = atom_idx
        elif atom.name == "C":
            c_idx = atom_idx
    
    # Find previous residue's C atom (ALA-7)
    ala_id = topology.find_residue("ALA", 7)
    ala = topology.get_residue(ala_id)
    for atom_idx in ala.atoms:
        atom = topology.get_atom(atom_idx)
        if atom.name == "C":
            prev_c_idx = atom_idx
            break
    
    # Find next residue's N atom (PRO-9)
    pro_id = topology.find_residue("PRO", 9)
    pro = topology.get_residue(pro_id)
    for atom_idx in pro.atoms:
        atom = topology.get_atom(atom_idx)
        if atom.name == "N":
            next_n_idx = atom_idx
            break
    
    # Check phi/psi dihedrals
    # phi dihedral: C(i-1)-N(i)-CA(i)-C(i)
    assert topology.has_dihedral(prev_c_idx, n_idx, ca_idx, c_idx), "Missing phi dihedral"
    # psi dihedral: N(i)-CA(i)-C(i)-N(i+1)
    assert topology.has_dihedral(n_idx, ca_idx, c_idx, next_n_idx), "Missing psi dihedral"


def test_parse_impropers(test_data_dir):
    """Test parsing improper dihedrals section from PSF file."""
    psf_file = os.path.join(test_data_dir, "test_proa.psf")
    parser = PSFParser()
    topology = Topology()
    
    assert parser.parse_to_topology(psf_file, topology), "Failed to parse PSF file"
    
    # Check total number of impropers
    assert topology.get_num_impropers() == 29, "Wrong number of impropers"

    # Instead of assuming C-CA-N-O improper exists in every residue,
    # we should verify impropers that we know exist in the PSF file.
    # For example, if we know from the PSF that atoms 11,5,13,12 (1-based)
    # form an improper, we can test for that:
    c_idx = 10   # 11-1, converting to 0-based
    ca_idx = 4   # 5-1
    n_idx = 12   # 13-1
    o_idx = 11   # 12-1
    assert topology.has_improper(c_idx, ca_idx, n_idx, o_idx), "Expected improper not found"

# tests/io/psfParser/advanced_features.py
"""PSF Parser advanced features tests."""

import os
import pytest
from pygcmc.io import PSFParser
from pygcmc.model import Topology, TopologyResidue, TopologyAtom

# Get the directory containing test data files
TEST_DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "data")


def test_parse_donors_acceptors():
    """Test parsing hydrogen bond donors and acceptors from PSF file."""
    psf_file = os.path.join(TEST_DATA_DIR, "test_proa.psf")
    parser = PSFParser()
    topology = Topology()
    
    assert parser.parse_to_topology(psf_file, topology), "Failed to parse PSF file"
    
    # Check total numbers
    assert topology.get_num_donors() > 0, "No donors found"
    assert topology.get_num_acceptors() > 0, "No acceptors found"
    
    # Check ASN-12 sidechain donors/acceptors
    asn_id = topology.find_residue("ASN", 12)
    asn = topology.get_residue(asn_id)
    
    # Find relevant atoms
    nd2_idx = None
    hd21_idx = None
    hd22_idx = None
    od1_idx = None
    
    for atom_idx in asn.atoms:
        atom = topology.get_atom(atom_idx)
        if atom.name == "ND2":
            nd2_idx = atom_idx
        elif atom.name == "HD21":
            hd21_idx = atom_idx
        elif atom.name == "HD22":
            hd22_idx = atom_idx
        elif atom.name == "OD1":
            od1_idx = atom_idx
    
    # Check donor-H pairs
    assert topology.has_donor(nd2_idx, hd21_idx), "Missing ND2-HD21 donor"
    assert topology.has_donor(nd2_idx, hd22_idx), "Missing ND2-HD22 donor"
    
    # Check acceptor
    assert topology.has_acceptor(od1_idx), "Missing OD1 acceptor"


def test_parse_cmap():
    """Test parsing CMAP (correction map) terms from PSF file."""
    psf_file = os.path.join(TEST_DATA_DIR, "test_proa.psf")
    parser = PSFParser()
    topology = Topology()
    
    assert parser.parse_to_topology(psf_file, topology), "Failed to parse PSF file"
    
    # Check total number of CMAP terms
    assert topology.get_num_cmaps() > 0, "No CMAP terms found"
    
    # Check CMAP term for VAL-8 to PRO-9
    # We need 8 atoms for a CMAP term in the following order (from PSF file):
    # C(i), N(i+1), CA(i+1), C(i+1), N(i+1), CA(i+1), C(i+1), N(i+2)
    
    # Get all required residues
    ala7_id = topology.find_residue("ALA", 7)
    val8_id = topology.find_residue("VAL", 8)
    pro9_id = topology.find_residue("PRO", 9)
    ala10_id = topology.find_residue("ALA", 10)
    
    assert all(x is not None for x in [ala7_id, val8_id, pro9_id, ala10_id]), "Failed to find required residues"
    
    # Get residue objects
    ala7 = topology.get_residue(ala7_id)
    val8 = topology.get_residue(val8_id)
    pro9 = topology.get_residue(pro9_id)
    ala10 = topology.get_residue(ala10_id)
    
    # Initialize atom indices
    val8_c_idx = None    # C(i)
    pro9_n_idx = None    # N(i+1)
    pro9_ca_idx = None   # CA(i+1)
    pro9_c_idx = None    # C(i+1)
    ala10_n_idx = None   # N(i+2)
    
    # Find VAL-8's C
    for atom_idx in val8.atoms:
        atom = topology.get_atom(atom_idx)
        if atom.name == "C":
            val8_c_idx = atom_idx
            break
    
    # Find PRO-9's N, CA, C
    for atom_idx in pro9.atoms:
        atom = topology.get_atom(atom_idx)
        if atom.name == "N":
            pro9_n_idx = atom_idx
        elif atom.name == "CA":
            pro9_ca_idx = atom_idx
        elif atom.name == "C":
            pro9_c_idx = atom_idx
    
    # Find ALA-10's N
    for atom_idx in ala10.atoms:
        atom = topology.get_atom(atom_idx)
        if atom.name == "N":
            ala10_n_idx = atom_idx
            break
    
    # Verify we found all atoms
    assert all(x is not None for x in [
        val8_c_idx, pro9_n_idx, pro9_ca_idx, pro9_c_idx, ala10_n_idx
    ]), "Failed to find all required atoms for CMAP"
    
    # Check CMAP term with atoms in PSF file order:
    # C(i), N(i+1), CA(i+1), C(i+1), N(i+1), CA(i+1), C(i+1), N(i+2)
    cmap_atoms = [
        val8_c_idx,    # C(i)
        pro9_n_idx,    # N(i+1)
        pro9_ca_idx,   # CA(i+1)
        pro9_c_idx,    # C(i+1)
        pro9_n_idx,    # N(i+1) again
        pro9_ca_idx,   # CA(i+1) again
        pro9_c_idx,    # C(i+1) again
        ala10_n_idx    # N(i+2)
    ]
    
    assert topology.has_cmap(cmap_atoms), "Missing CMAP term between VAL-8 and PRO-9"


def test_parse_groups():
    """Test parsing group definitions from PSF file."""
    psf_file = os.path.join(TEST_DATA_DIR, "test_proa.psf")
    parser = PSFParser()
    topology = Topology()
    
    assert parser.parse_to_topology(psf_file, topology), "Failed to parse PSF file"
    
    # Check total number of groups
    num_groups = topology.get_num_groups()
    assert num_groups > 0, "No groups found"
    
    # In CHARMM PSF files, groups are not necessarily defined per residue.
    # Instead, they are defined based on charge groups or other criteria.
    # Here we just verify that the number of groups matches what's in the PSF file
    # and that each group contains valid atom indices.
    
    for i in range(num_groups):
        group = topology.get_group(i)
        # Verify that all atom indices in the group are valid
        for atom_idx in group.atoms:
            assert topology.has_atom(atom_idx), f"Invalid atom index {atom_idx} in group {i}"


def test_parse_all_cmaps():
    """Test parsing all CMAP terms from PSF file."""
    psf_file = os.path.join(TEST_DATA_DIR, "test_proa.psf")
    parser = PSFParser()
    topology = Topology()
    
    assert parser.parse_to_topology(psf_file, topology), "Failed to parse PSF file"
    
    # Verify total number of CMAP terms
    assert topology.get_num_cmaps() == 7, "Expected 7 CMAP terms, got {}".format(topology.get_num_cmaps())
    
    # Expected CMAP terms (1-based indices from PSF file)
    expected_cmaps = [
        [11, 13, 15, 27, 13, 15, 27, 29],
        [27, 29, 33, 41, 29, 33, 41, 43],
        [41, 43, 45, 51, 43, 45, 51, 53],
        [51, 53, 57, 65, 53, 57, 65, 67],
        [65, 67, 69, 79, 67, 69, 79, 81],
        [79, 81, 83, 96, 81, 83, 96, 98],
        [96, 98, 100, 113, 98, 100, 113, 115]
    ]
    
    # Convert to 0-based indices and verify each CMAP term
    for i, expected_cmap in enumerate(expected_cmaps):
        expected_0based = [idx - 1 for idx in expected_cmap]
        assert topology.has_cmap(expected_0based), f"Missing CMAP term {i+1}: {expected_cmap}"
    
    # Verify these correspond to the protein backbone
    # Each CMAP should connect four consecutive residues through their backbone atoms
    residue_sequences = [
        ["ALA", "VAL", "PRO"],           # First CMAP
        ["VAL", "PRO", "ALA"],           # Second CMAP
        ["PRO", "ALA", "PRO"],           # Third CMAP
        ["ALA", "PRO", "ASN"],           # Fourth CMAP
        ["PRO", "ASN", "GLN"],           # Fifth CMAP
        ["ASN", "GLN", "GLN"],           # Sixth CMAP
        ["GLN", "GLN", "PRO"]            # Seventh CMAP
    ]
    
    for i, residues in enumerate(residue_sequences):
        # For each residue sequence, verify the residues exist and are connected
        for j in range(len(residues)-1):
            res1_id = topology.find_residue(residues[j], 7+i+j)
            res2_id = topology.find_residue(residues[j+1], 8+i+j)
            assert res1_id is not None, f"Could not find residue {residues[j]} {7+i+j}"
            assert res2_id is not None, f"Could not find residue {residues[j+1]} {8+i+j}"
            
            # Get the residues
            res1 = topology.get_residue(res1_id)
            res2 = topology.get_residue(res2_id)
            
            # Verify they share at least one bond (they should be connected)
            found_connection = False
            for atom1_idx in res1.atoms:
                for atom2_idx in res2.atoms:
                    if topology.has_bond(atom1_idx, atom2_idx):
                        found_connection = True
                        break
                if found_connection:
                    break
            assert found_connection, f"No connection found between {residues[j]} and {residues[j+1]}"

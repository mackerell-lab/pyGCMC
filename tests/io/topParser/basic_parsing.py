# tests/io/topParser/basic_parsing.py
"""TOP Parser basic parsing tests."""

import os
import pytest
from pygcmc.io import TOPParser
from pygcmc.model import Topology, TopologyResidue, TopologyAtom


def test_parse_protein_top(test_data_dir):
    """Test parsing protein topology file (test.top)."""
    top_file = os.path.join(test_data_dir, "test.top")
    parser = TOPParser()
    topology = Topology()
    
    # Parse the topology file
    assert parser.parse_to_topology(top_file, topology), "Failed to parse protein topology file"
    
    # Verify basic topology information
    assert topology.get_num_atoms() == 196, "Wrong total number of atoms (should be 196: 129 protein + 13 BENX + 24 PRPX + 30 SOL)"
    assert topology.get_num_residues() > 0, "No residues found in topology"
    assert topology.get_num_segments() > 0, "No segments found in topology"
    assert topology.get_num_bonds() == 163, "Wrong number of bonds (should be 163: protein + BENX + 2×PRPX)"
    assert topology.get_num_angles() == 297, "Wrong number of angles (should be 297: protein[243] + BENX + 2×PRPX[18×2])"
    assert topology.get_num_dihedrals() == 422, "Wrong number of dihedrals (should be 422: protein[362] + BENX[24] + 2×PRPX[18×2])"
    assert topology.get_num_impropers() == 29, "Wrong number of impropers (all from protein)"
    assert topology.get_num_cmaps() == 7, "Wrong number of CMAPs (should be 7: one for each amino acid pair in protein except the last one)"
    
    # Check specific residues
    residues = {"ALA", "VAL", "PRO", "ASN", "GLN"}
    for residue in residues:
        found = False
        for i in range(topology.get_num_residues()):
            if topology.get_residue(i).name == residue:
                found = True
                break
        assert found, f"Residue {residue} not found in topology"
    
    # Check N-terminal ALA (residue 7)
    ala_id = topology.find_residue("ALA", 7)
    assert ala_id is not None, "N-terminal ALA not found"
    ala = topology.get_residue(ala_id)
    
    # Check N-terminal atoms
    n_atom = None
    h_atoms = []
    for atom_idx in ala.atoms:
        atom = topology.get_atom(atom_idx)
        if atom.name == "N":
            n_atom = atom
        elif atom.name in ["H1", "H2", "H3"]:
            h_atoms.append(atom)
    
    assert n_atom is not None, "N atom not found in N-terminal ALA"
    assert len(h_atoms) == 3, "Missing H atoms in N-terminal ALA"
    assert n_atom.type == "NH3", "Wrong type for N-terminal N atom"
    assert abs(n_atom.charge + 0.30) < 1e-6, "Wrong charge for N-terminal N atom"
    for h in h_atoms:
        assert h.type == "HC", f"Wrong type for {h.name}"
        assert abs(h.charge - 0.33) < 1e-6, f"Wrong charge for {h.name}"


def test_parse_step1_top(test_data_dir):
    """Test parsing step1_pdbreader.top file (which includes step1_pdbreader.itp)."""
    top_file = os.path.join(test_data_dir, "step1_pdbreader.top")
    parser = TOPParser()
    topology = Topology()
    
    # Parse the "step1_pdbreader.top" file
    assert parser.parse_to_topology(top_file, topology), "Failed to parse step1_pdbreader.top"

    # In the [molecules] section, there are 11 repeats of:
    #   rna1 1
    #   MG 5
    #   TIP3 50
    #
    # So there should be:
    #   - 11 RNA segments (we won't check residue counts within each segment, but you could expand to do so)
    #   - 11 × 5 = 55 total MG residues
    #   - 11 × 50 = 550 total TIP3 residues
    
    mg_count = 0
    tip3_count = 0
    rna_segments = 0

    # Count how many residues are named "MG" or "TIP3"
    for i in range(topology.get_num_residues()):
        res = topology.get_residue(i)
        if res.name == "MG":
            mg_count += 1
        elif res.name == "TIP3":
            tip3_count += 1

    # Count RNA segments by looking at residue segment names
    # We'll consider a new RNA segment starts when we see a residue with a different segment name
    current_segment = None
    for i in range(topology.get_num_residues()):
        res = topology.get_residue(i)
        if res.segment != current_segment and res.name not in ["MG", "TIP3"]:
            rna_segments += 1
            current_segment = res.segment

    assert mg_count == 60, f"Expected 60 MG residues, found {mg_count}"
    assert tip3_count == 600, f"Expected 600 TIP3 residues, found {tip3_count}"
    assert rna_segments == 12, f"Expected 12 RNA segments, found {rna_segments}"
    
    # You can add additional checks (bonds, angles, dihedrals, charges) as desired.
    # For example, check total atom count or partial checks:
    # assert topology.get_num_atoms() == <some_expected_number>
    # ...

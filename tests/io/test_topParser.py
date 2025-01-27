# tests/io/test_topParser.py

import os
import pytest
from pygcmc.io import TopParser
from pygcmc.model import Topology, TopologyResidue, TopologyAtom

@pytest.fixture
def test_data_dir():
    """Get the path to the test data directory."""
    current_dir = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(current_dir, "..", "data")

def test_parse_protein_top(test_data_dir):
    """Test parsing protein topology file (test.top)."""
    top_file = os.path.join(test_data_dir, "test.top")
    parser = TopParser()
    topology = Topology()
    
    # Parse the topology file
    assert parser.parse_to_topology(top_file, topology), "Failed to parse protein topology file"
    
    # Verify basic topology information
    assert topology.get_num_atoms() == 129, "Wrong number of atoms"
    assert topology.get_num_residues() > 0, "No residues found in topology"
    assert topology.get_num_segments() > 0, "No segments found in topology"
    assert topology.get_num_bonds() == 131, "Wrong number of bonds"
    assert topology.get_num_angles() == 243, "Wrong number of angles"
    assert topology.get_num_dihedrals() == 362, "Wrong number of dihedrals"
    assert topology.get_num_impropers() == 29, "Wrong number of impropers"
    
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

def test_parse_solvent_top(test_data_dir):
    """Test parsing solvent topology file (sol.itp)."""
    top_file = os.path.join(test_data_dir, "mols", "sol.itp")
    parser = TopParser()
    topology = Topology()
    
    # Parse the topology file
    assert parser.parse_to_topology(top_file, topology), "Failed to parse solvent topology file"
    
    # Find water residues
    water_residues = set()
    for i in range(topology.get_num_residues()):
        res = topology.get_residue(i)
        if res.name in {"SOL", "WAT", "HOH", "TIP3"}:
            water_residues.add(i)
    
    assert len(water_residues) > 0, "No water residues found"
    
    # Check each water residue
    for res_id in water_residues:
        res = topology.get_residue(res_id)
        assert len(res.atoms) == 3, f"Water residue {res.name} {res.number} should have 3 atoms"
        
        # Find OW and HW atoms
        ow = None
        hw = []
        for atom_idx in res.atoms:
            atom = topology.get_atom(atom_idx)
            if atom.name == "OW":
                ow = atom
            elif atom.name in ["HW1", "HW2"]:
                hw.append(atom)
        
        # Verify water atoms
        assert ow is not None, f"OW atom not found in {res.name} {res.number}"
        assert len(hw) == 2, f"Wrong number of HW atoms in {res.name} {res.number}"
        
        # Check topology properties
        assert ow.type == "OT", "Wrong type for OW atom"
        assert abs(ow.charge + 0.834) < 1e-6, "Wrong charge for OW atom"
        assert abs(ow.mass - 15.9994) < 1e-6, "Wrong mass for OW atom"
        
        for h in hw:
            assert h.type == "HT", f"Wrong type for {h.name}"
            assert abs(h.charge - 0.417) < 1e-6, f"Wrong charge for {h.name}"
            assert abs(h.mass - 1.0080) < 1e-6, f"Wrong mass for {h.name}"

def test_parse_benzene_top(test_data_dir):
    """Test parsing benzene topology file (benx.itp)."""
    top_file = os.path.join(test_data_dir, "mols", "benx.itp")
    parser = TopParser()
    topology = Topology()
    
    # Parse the topology file
    assert parser.parse_to_topology(top_file, topology), "Failed to parse benzene topology file"
    
    # Find benzene residues
    benx_residues = set()
    for i in range(topology.get_num_residues()):
        res = topology.get_residue(i)
        if res.name == "BENX":
            benx_residues.add(i)
    
    assert len(benx_residues) > 0, "No benzene residues found"
    
    # Check each benzene residue
    for res_id in benx_residues:
        res = topology.get_residue(res_id)
        expected_atoms = {
            "CG", "CD1", "CD2", "CE1", "CE2", "CZ",
            "HG", "HD1", "HD2", "HE1", "HE2", "HZ", "LPA"
        }
        found_atoms = set()
        
        for atom_idx in res.atoms:
            atom = topology.get_atom(atom_idx)
            found_atoms.add(atom.name)
            
            # Check specific atoms
            if atom.name == "CG":
                assert atom.type == "CG2R61", "Wrong type for CG atom"
                assert abs(atom.charge + 0.115) < 1e-6, "Wrong charge for CG atom"
            elif atom.name.startswith("CD"):
                assert atom.type == "CG2R61", "Wrong type for CD atom"
            elif atom.name.startswith("CE"):
                assert atom.type == "CG2R61", "Wrong type for CE atom"
            elif atom.name == "CZ":
                assert atom.type == "CG2R61", "Wrong type for CZ atom"
            elif atom.name.startswith("H"):
                assert atom.type == "HGR61", "Wrong type for H atom"
        
        assert found_atoms == expected_atoms, f"Missing or extra atoms in benzene residue"

def test_parse_propane_top(test_data_dir):
    """Test parsing propane topology file (prpx.itp)."""
    top_file = os.path.join(test_data_dir, "mols", "prpx.itp")
    parser = TopParser()
    topology = Topology()
    
    # Parse the topology file
    assert parser.parse_to_topology(top_file, topology), "Failed to parse propane topology file"
    
    # Find propane residues
    prpx_residues = set()
    for i in range(topology.get_num_residues()):
        res = topology.get_residue(i)
        if res.name == "PRPX":
            prpx_residues.add(i)
    
    assert len(prpx_residues) > 0, "No propane residues found"
    
    # Check each propane residue
    for res_id in prpx_residues:
        res = topology.get_residue(res_id)
        expected_atoms = {
            "C1", "C2", "C3",
            "H11", "H12", "H13",
            "H21", "H22",
            "H31", "H32", "H33",
            "LPA"
        }
        found_atoms = set()
        
        for atom_idx in res.atoms:
            atom = topology.get_atom(atom_idx)
            found_atoms.add(atom.name)
            
            # Check specific atoms
            if atom.name == "C1":
                assert atom.type == "CG331", "Wrong type for C1 atom"
                assert abs(atom.charge + 0.27) < 1e-6, "Wrong charge for C1 atom"
            elif atom.name == "C2":
                assert atom.type == "CG321", "Wrong type for C2 atom"
            elif atom.name == "C3":
                assert atom.type == "CG331", "Wrong type for C3 atom"
            elif atom.name.startswith("H"):
                if atom.name.startswith("H2"):
                    assert atom.type == "HGA2", "Wrong type for H2x atom"
                else:
                    assert atom.type == "HGA3", "Wrong type for H1x/H3x atom"
        
        assert found_atoms == expected_atoms, f"Missing or extra atoms in propane residue"

def test_parse_bonds(test_data_dir):
    """Test parsing bonds section from topology file."""
    top_file = os.path.join(test_data_dir, "test.top")
    parser = TopParser()
    topology = Topology()
    
    assert parser.parse_to_topology(top_file, topology), "Failed to parse topology file"
    
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
    """Test parsing angles section from topology file."""
    top_file = os.path.join(test_data_dir, "test.top")
    parser = TopParser()
    topology = Topology()
    
    assert parser.parse_to_topology(top_file, topology), "Failed to parse topology file"
    
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
    """Test parsing dihedrals section from topology file."""
    top_file = os.path.join(test_data_dir, "test.top")
    parser = TopParser()
    topology = Topology()
    
    assert parser.parse_to_topology(top_file, topology), "Failed to parse topology file"
    
    # Check total number of dihedrals
    assert topology.get_num_dihedrals() == 362, "Wrong number of dihedrals"
    
    # Check backbone dihedrals in VAL-8
    val_id = topology.find_residue("VAL", 8)
    val = topology.get_residue(val_id)
    
    # Find backbone atoms for phi/psi angles
    n_idx = None
    ca_idx = None
    c_idx = None
    
    # Find VAL-8 backbone atoms
    for atom_idx in val.atoms:
        atom = topology.get_atom(atom_idx)
        if atom.name == "N":
            n_idx = atom_idx
        elif atom.name == "CA":
            ca_idx = atom_idx
        elif atom.name == "C":
            c_idx = atom_idx
    
    # Get previous residue's C atom
    ala_id = topology.find_residue("ALA", 7)
    ala = topology.get_residue(ala_id)
    prev_c_idx = None
    for atom_idx in ala.atoms:
        atom = topology.get_atom(atom_idx)
        if atom.name == "C":
            prev_c_idx = atom_idx
            break
    
    # Check phi dihedral
    assert topology.has_dihedral(prev_c_idx, n_idx, ca_idx, c_idx), "Missing phi dihedral"

def test_parse_nonexistent_file(test_data_dir):
    """Test parsing a non-existent topology file."""
    top_file = os.path.join(test_data_dir, "nonexistent.top")
    parser = TopParser()
    topology = Topology()
    
    assert not parser.parse_to_topology(top_file, topology), "Should fail for non-existent file"

def test_parse_invalid_top(test_data_dir, tmp_path):
    """Test parsing an invalid topology file."""
    # Create an invalid topology file
    invalid_top = tmp_path / "invalid.top"
    with open(invalid_top, "w") as f:
        f.write("This is not a topology file\n")
    
    parser = TopParser()
    topology = Topology()
    assert not parser.parse_to_topology(str(invalid_top), topology), "Should fail for invalid topology file"


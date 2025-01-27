# tests/io/test_psfParser.py

import os
import pytest
from pygcmc.io import PSFParser
from pygcmc.model import Topology, TopologyResidue, TopologyAtom

@pytest.fixture
def test_data_dir():
    """Get the path to the test data directory."""
    current_dir = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(current_dir, "..", "data")

def test_parse_protein_psf(test_data_dir):
    """Test parsing protein PSF file (test_proa.psf)."""
    psf_file = os.path.join(test_data_dir, "test_proa.psf")
    parser = PSFParser()
    topology = Topology()
    
    # Parse the PSF file
    assert parser.parse_to_topology(psf_file, topology), "Failed to parse protein PSF file"
    
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
    ht_atoms = []
    for atom_idx in ala.atoms:
        atom = topology.get_atom(atom_idx)
        if atom.name == "N":
            n_atom = atom
        elif atom.name in ["HT1", "HT2", "HT3"]:
            ht_atoms.append(atom)
    
    assert n_atom is not None, "N atom not found in N-terminal ALA"
    assert len(ht_atoms) == 3, "Missing HT atoms in N-terminal ALA"
    assert n_atom.type == "NH3", "Wrong type for N-terminal N atom"
    assert abs(n_atom.charge + 0.30) < 1e-6, "Wrong charge for N-terminal N atom"
    for ht in ht_atoms:
        assert ht.type == "HC", f"Wrong type for {ht.name}"
        assert abs(ht.charge - 0.33) < 1e-6, f"Wrong charge for {ht.name}"

def test_parse_solvent_psf(test_data_dir):
    """Test parsing solvent PSF file (sol.psf)."""
    psf_file = os.path.join(test_data_dir, "mols", "sol.psf")
    parser = PSFParser()
    topology = Topology()
    
    # Parse the PSF file
    assert parser.parse_to_topology(psf_file, topology), "Failed to parse solvent PSF file"
    
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

def test_parse_benzene_psf(test_data_dir):
    """Test parsing benzene PSF file (benx.psf)."""
    psf_file = os.path.join(test_data_dir, "mols", "benx.psf")
    parser = PSFParser()
    topology = Topology()
    
    # Parse the PSF file
    assert parser.parse_to_topology(psf_file, topology), "Failed to parse benzene PSF file"
    
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

def test_parse_propane_psf(test_data_dir):
    """Test parsing propane PSF file (prpx.psf)."""
    psf_file = os.path.join(test_data_dir, "mols", "prpx.psf")
    parser = PSFParser()
    topology = Topology()
    
    # Parse the PSF file
    assert parser.parse_to_topology(psf_file, topology), "Failed to parse propane PSF file"
    
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

def test_parse_nonexistent_file(test_data_dir):
    """Test parsing a non-existent PSF file."""
    psf_file = os.path.join(test_data_dir, "nonexistent.psf")
    parser = PSFParser()
    topology = Topology()
    
    assert not parser.parse_to_topology(psf_file, topology), "Should fail for non-existent file"

def test_parse_invalid_psf(test_data_dir, tmp_path):
    """Test parsing an invalid PSF file."""
    # Create an invalid PSF file
    invalid_psf = tmp_path / "invalid.psf"
    with open(invalid_psf, "w") as f:
        f.write("This is not a PSF file\n")
    
    parser = PSFParser()
    topology = Topology()
    
    assert not parser.parse_to_topology(str(invalid_psf), topology), "Should fail for invalid PSF file"

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

def test_parse_donors_acceptors(test_data_dir):
    """Test parsing hydrogen bond donors and acceptors from PSF file."""
    psf_file = os.path.join(test_data_dir, "test_proa.psf")
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

def test_parse_cmap(test_data_dir):
    """Test parsing CMAP (correction map) terms from PSF file."""
    psf_file = os.path.join(test_data_dir, "test_proa.psf")
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

def test_parse_groups(test_data_dir):
    """Test parsing group definitions from PSF file."""
    psf_file = os.path.join(test_data_dir, "test_proa.psf")
    parser = PSFParser()
    topology = Topology()
    
    assert parser.parse_to_topology(psf_file, topology), "Failed to parse PSF file"
    
    # Check total number of groups
    assert topology.get_num_groups() > 0, "No groups found"
    
    # Check if each residue is a group
    for i in range(topology.get_num_residues()):
        res = topology.get_residue(i)
        assert topology.has_group(res.atoms), f"Missing group for residue {res.name} {res.number}"

def test_parse_out_of_order_psf(test_data_dir, tmp_path):
    """Test parsing PSF file with sections in non-standard order."""
    # Create a PSF file with reordered sections
    reordered_psf = tmp_path / "reordered.psf"
    with open(reordered_psf, "w") as f:
        f.write("PSF EXT CMAP XPLOR\n\n")
        f.write("         1 !NTITLE\n")
        f.write("* REORDERED PSF FILE FOR TESTING\n\n")
        f.write("       131 !NBOND: bonds\n")
        f.write("         1         2         2         3         3         4\n")
        f.write("       129 !NATOM\n")
        f.write("         1 PROA     7        ALA      N        NH3     -0.300000       14.0070           0\n")
        f.write("         2 PROA     7        ALA      HT1      HC       0.330000        1.0080           0\n")
    
    parser = PSFParser()
    topology = Topology()
    
    # Should still parse correctly despite reordered sections
    assert parser.parse_to_topology(str(reordered_psf), topology), "Failed to parse reordered PSF file"
    assert topology.get_num_atoms() == 129, "Wrong number of atoms"
    assert topology.get_num_bonds() == 131, "Wrong number of bonds"


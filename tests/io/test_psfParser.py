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
    assert topology.get_num_atoms() > 0, "No atoms found in topology"
    assert topology.get_num_residues() > 0, "No residues found in topology"
    assert topology.get_num_segments() > 0, "No segments found in topology"
    
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


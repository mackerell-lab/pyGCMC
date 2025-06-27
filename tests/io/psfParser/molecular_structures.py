# tests/io/psfParser/molecular_structures.py
"""PSF Parser molecular structure tests."""

import os
import pytest
from pygcmc.io import PSFParser
from pygcmc.model import Topology, TopologyResidue, TopologyAtom


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


def test_parse_step1_psf(test_data_dir):
    """Test parsing step1_pdbreader.psf file (which includes MG and TIP3 molecules) using PSFParser."""
    psf_file = os.path.join(test_data_dir, "step1_pdbreader.psf")
    parser = PSFParser()
    topology = Topology()
    
    # Parse the PSF file
    assert parser.parse_to_topology(psf_file, topology), "Failed to parse step1_pdbreader.psf"
    
    # Verify the total number of atoms according to the PSF header (e.g., 22068 atoms)
    expected_atoms = 22068
    actual_atoms = topology.get_num_atoms()
    assert actual_atoms == expected_atoms, f"Expected {expected_atoms} atoms, found {actual_atoms}"
    
    mg_count = 0
    tip3_count = 0
    rna_segments = 0
    
    # Print all segments and their residues for debugging
    print("\nSegments and their residues:")
    for i in range(topology.get_num_segments()):
        segment = topology.get_segment(i)
        print(f"\nSegment {segment.name}:")
        for res_idx in segment.residues:
            res = topology.get_residue(res_idx)
            print(f"  Residue {res.name} {res.number} (segment: {res.segment})")
            if res.name == "MG":
                mg_count += 1
            elif res.name == "TIP3":
                tip3_count += 1
    
    # Print total counts
    print(f"\nTotal counts:")
    print(f"MG residues: {mg_count}")
    print(f"TIP3 residues: {tip3_count}")
    print(f"Total residues: {topology.get_num_residues()}")
    print(f"Total segments: {topology.get_num_segments()}")
    
    # Count RNA segments
    current_segment = None
    for i in range(topology.get_num_residues()):
        res = topology.get_residue(i)
        if res.name not in {"MG", "TIP3"} and res.segment != current_segment:
            rna_segments += 1
            current_segment = res.segment
    
    print(f"RNA segments: {rna_segments}")
    
    # Based on the system setup (as in the TOP test):
    # Expected RNA segments: 12, MG residues: 60, TIP3 residues: 600
    assert mg_count == 60, f"Expected 60 MG residues, found {mg_count}"
    assert tip3_count == 600, f"Expected 600 TIP3 residues, found {tip3_count}"
    assert rna_segments == 12, f"Expected 12 RNA segments, found {rna_segments}"

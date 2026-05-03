# tests/io/psfParser/basic_parsing.py
"""PSF Parser basic parsing tests."""

import os
import pytest
from pygcmc.io import PSFParser
from pygcmc.model import Topology, TopologyResidue, TopologyAtom


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


def test_parse_psf_string_matches_file(test_data_dir):
    """Test parsing PSF content directly from memory."""
    psf_file = os.path.join(test_data_dir, "test_proa.psf")
    file_topology = PSFParser.parse_file(psf_file)

    with open(psf_file, "r") as handle:
        string_topology = PSFParser.parse_string(handle.read())

    assert string_topology.get_num_atoms() == file_topology.get_num_atoms()
    assert string_topology.get_num_bonds() == file_topology.get_num_bonds()
    assert string_topology.get_num_angles() == file_topology.get_num_angles()
    assert string_topology.get_num_dihedrals() == file_topology.get_num_dihedrals()
    assert string_topology.get_num_impropers() == file_topology.get_num_impropers()


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


def test_parse_out_of_order_psf(test_data_dir, tmp_path):
    """Test parsing PSF file with sections in non-standard order."""
    # Read the original PSF file
    original_psf = os.path.join(test_data_dir, "test_proa.psf")
    with open(original_psf, "r") as f:
        lines = f.readlines()

    # Create sections dictionary
    sections = {}
    current_section = None
    current_lines = []

    # Collect sections
    for line in lines:
        if "!NTITLE" in line:
            current_section = "NTITLE"
            current_lines = [line]
        elif "!NATOM" in line:
            if current_section:
                sections[current_section] = current_lines
            current_section = "NATOM"
            current_lines = [line]
        elif "!NBOND" in line:
            if current_section:
                sections[current_section] = current_lines
            current_section = "NBOND"
            current_lines = [line]
        elif "!NTHETA" in line:
            if current_section:
                sections[current_section] = current_lines
            current_section = "NTHETA"
            current_lines = [line]
        elif "!NPHI" in line:
            if current_section:
                sections[current_section] = current_lines
            current_section = "NPHI"
            current_lines = [line]
        elif "!NIMPHI" in line:
            if current_section:
                sections[current_section] = current_lines
            current_section = "NIMPHI"
            current_lines = [line]
        elif "!NDON" in line:
            if current_section:
                sections[current_section] = current_lines
            current_section = "NDON"
            current_lines = [line]
        elif "!NACC" in line:
            if current_section:
                sections[current_section] = current_lines
            current_section = "NACC"
            current_lines = [line]
        else:
            if current_section:
                current_lines.append(line)

    # Add the last section
    if current_section:
        sections[current_section] = current_lines

    # Create reordered PSF file
    reordered_psf = tmp_path / "reordered.psf"
    with open(reordered_psf, "w") as f:
        # Write header
        f.write("PSF EXT CMAP XPLOR\n\n")
        
        # Write title section
        f.writelines(sections["NTITLE"])
        
        # Write sections in non-standard order
        section_order = ["NBOND", "NTHETA", "NATOM", "NPHI", "NIMPHI", "NDON", "NACC", "CMAP"]
        for section in section_order:
            if section in sections:
                f.write("\n")  # Add spacing between sections
                f.writelines(sections[section])

    parser = PSFParser()
    topology = Topology()

    # Should still parse correctly despite reordered sections
    assert parser.parse_to_topology(str(reordered_psf), topology), "Failed to parse reordered PSF file"
    
    # Verify the topology has the correct content
    assert topology.get_num_atoms() == 129, "Wrong number of atoms"
    assert topology.get_num_bonds() == 131, "Wrong number of bonds"
    assert topology.get_num_angles() == 243, "Wrong number of angles"
    assert topology.get_num_dihedrals() == 362, "Wrong number of dihedrals"
    assert topology.get_num_impropers() == 29, "Wrong number of impropers"
    assert topology.get_num_cmaps() == 7, "Wrong number of CMAP terms"

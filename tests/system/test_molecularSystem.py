# tests/system/test_molecularSystem.py

import pytest
import os
import pygcmc
import math
from pygcmc.model import Topology
from pygcmc.io import TOPParser

# Get the directory containing test data files
TEST_DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")

@pytest.fixture
def test_structure():
    """Fixture to create a test Structure object."""
    pdb_path = os.path.join(TEST_DATA_DIR, "test.pdb")
    return pygcmc.PDBParser.parse_file(pdb_path)

@pytest.fixture
def test_topology():
    """Fixture to create a test Topology object."""
    top_path = os.path.join(TEST_DATA_DIR, "test.top")
    return TOPParser.parse_file(top_path)

def test_combine_structure_data(test_structure, test_topology):
    """Test if Structure data is correctly combined into Molecular object."""
    mol_system = pygcmc.MolecularSystem()
    molecular = mol_system.combine(test_structure, test_topology)
    
    # Check structure data
    assert len(molecular.atoms) == len(test_structure.atoms)
    assert len(molecular.residues) == len(test_structure.residues)
    
    # Check specific atom properties
    first_atom = molecular.atoms[0]
    first_struct_atom = test_structure.atoms[0]
    assert first_atom.get_type() == first_struct_atom.get_type()
    assert first_atom.get_resname() == first_struct_atom.get_resname()
    assert first_atom.get_chain() == first_struct_atom.get_chain()
    
    # Check box dimensions
    assert molecular.boxDimensions == pytest.approx(test_structure.box_dimensions)

def test_combine_topology_data(test_structure, test_topology):
    """Test if Topology data is correctly combined into Molecular object."""
    mol_system = pygcmc.MolecularSystem()
    molecular = mol_system.combine(test_structure, test_topology)
    
    # Check topology data
    assert len(molecular.topology_atoms) == test_topology.get_num_atoms()
    assert len(molecular.topology_residues) == test_topology.get_num_residues()
    assert len(molecular.segments) == test_topology.get_num_segments()
    
    # Check specific topology properties
    first_top_atom = molecular.topology_atoms[0]
    orig_top_atom = test_topology.get_atom(0)
    assert first_top_atom.name == orig_top_atom.name
    assert first_top_atom.type == orig_top_atom.type
    assert first_top_atom.charge == pytest.approx(orig_top_atom.charge)
    assert first_top_atom.mass == pytest.approx(orig_top_atom.mass)

def test_combine_bonds_and_angles(test_structure, test_topology):
    """Test if bonds and angles are correctly combined into Molecular object."""
    mol_system = pygcmc.MolecularSystem()
    molecular = mol_system.combine(test_structure, test_topology)
    
    # Check bonds
    assert molecular.get_num_bonds() == test_topology.get_num_bonds()
    if molecular.get_num_bonds() > 0:
        first_bond = molecular.bonds[0]
        assert len(first_bond) == 2  # Bond should connect 2 atoms
        assert all(isinstance(idx, int) for idx in first_bond)
    
    # Check angles
    assert molecular.get_num_angles() == test_topology.get_num_angles()
    if molecular.get_num_angles() > 0:
        first_angle = molecular.angles[0]
        assert len(first_angle) == 3  # Angle should involve 3 atoms
        assert all(isinstance(idx, int) for idx in first_angle)

def test_combine_mapping_data(test_structure, test_topology):
    """Test if mapping data is correctly combined into Molecular object."""
    mol_system = pygcmc.MolecularSystem()
    molecular = mol_system.combine(test_structure, test_topology)
    
    # Check segment mapping
    assert len(molecular.segment_map) == test_topology.get_num_segments()
    
    # Check residue mapping
    assert len(molecular.residue_map) == test_topology.get_num_residues()
    
    # Check atom mapping
    assert len(molecular.atom_map) == test_topology.get_num_atoms()
    
    # Check specific mapping
    if test_topology.get_num_atoms() > 0:
        first_atom = test_topology.get_atom(0)
        first_res = test_topology.get_residue(first_atom.residue_id)
        key = (first_res.name, first_res.number, first_atom.name)
        assert key in molecular.atom_map
        assert molecular.atom_map[key] == 0

def test_combine_null_inputs():
    """Test handling of null inputs."""
    mol_system = pygcmc.MolecularSystem()
    
    with pytest.raises(ValueError):
        mol_system.combine(None, None)
    
    with pytest.raises(ValueError):
        mol_system.combine(pygcmc.Structure(), None)
    
    with pytest.raises(ValueError):
        mol_system.combine(None, TOPParser.parse_string(""))

def test_combine_mismatched_data(test_structure):
    """Test handling of mismatched Structure and Topology data."""
    mol_system = pygcmc.MolecularSystem()
    
    # Create a topology with different number of atoms
    topology = TOPParser.parse_file(os.path.join(TEST_DATA_DIR, "mols", "sol.itp"))
    topology.add_atom("CA", "CT", 0.0, 12.01, "ALA", 1, "PROT")
    topology.add_atom("CB", "CT", 0.0, 12.01, "ALA", 1, "PROT")
    
    with pytest.raises(RuntimeError):
        mol_system.combine(test_structure, topology)

def test_combine_empty_data():
    """Test combining empty Structure and Topology objects."""
    mol_system = pygcmc.MolecularSystem()
    structure = pygcmc.Structure()
    topology = TOPParser.parse_string("")
    
    molecular = mol_system.combine(structure, topology)
    
    assert molecular.get_num_atoms() == 0
    assert molecular.get_num_residues() == 0
    assert molecular.get_num_segments() == 0
    assert molecular.get_num_bonds() == 0
    assert molecular.get_num_angles() == 0
    assert len(molecular.segment_map) == 0
    assert len(molecular.residue_map) == 0
    assert len(molecular.atom_map) == 0

def test_combine_incompatible_files():
    """Test error handling when combining incompatible structure and topology files."""
    # Load test.pdb which contains protein, BENX, PRPX, and SOL
    pdb_path = os.path.join(TEST_DATA_DIR, "test.pdb")
    structure = pygcmc.PDBParser.parse_file(pdb_path)
    
    # Load test_proa.psf which only contains protein chain
    psf_path = os.path.join(TEST_DATA_DIR, "test_proa.psf")
    topology = pygcmc.io.PSFParser.parse_file(psf_path)
    
    # Attempt to combine should raise an error
    mol_system = pygcmc.MolecularSystem()
    with pytest.raises(RuntimeError) as excinfo:
        molecular = mol_system.combine(structure, topology)
    
    # Check that the error message is descriptive
    error_msg = str(excinfo.value)
    assert "Inconsistent total number of" in error_msg
    assert "Structure has" in error_msg
    assert "but Topology has" in error_msg

def test_combine_multiple_topologies():
    """Test combining structure with multiple topology files."""
    # Load structure file
    pdb_path = os.path.join(TEST_DATA_DIR, "test.pdb")
    structure = pygcmc.PDBParser.parse_file(pdb_path)
    
    # Load topology files
    proa_path = os.path.join(TEST_DATA_DIR, "test_proa.psf")
    benx_path = os.path.join(TEST_DATA_DIR, "mols/benx.psf")
    prpx_path = os.path.join(TEST_DATA_DIR, "mols/prpx.itp")
    sol_path = os.path.join(TEST_DATA_DIR, "mols/sol.psf")
    
    topologies = [
        pygcmc.io.PSFParser.parse_file(proa_path),
        pygcmc.io.PSFParser.parse_file(benx_path),
        pygcmc.io.TOPParser.parse_file(prpx_path),
        pygcmc.io.PSFParser.parse_file(sol_path)
    ]
    
    # Combine structure with multiple topologies
    mol_system = pygcmc.MolecularSystem()
    molecular = mol_system.combine_multiple(structure, topologies)
    
    # Verify the combined result
    assert len(molecular.atoms) == len(structure.atoms)
    assert len(molecular.residues) == len(structure.residues)
    
    # Check specific residues
    residues = molecular.residues
    assert residues[0].get_resname() == "ALA"  # From test_proa.psf
    assert residues[1].get_resname() == "VAL"  # From test_proa.psf
    assert residues[2].get_resname() == "PRO"  # From test_proa.psf
    
    # Find BENX residue
    benx_found = False
    for res in residues:
        if res.get_resname() == "BENX":
            benx_found = True
            break
    assert benx_found
    
    # Find PRPX residues
    prpx_count = sum(1 for res in residues if res.get_resname() == "PRPX")
    assert prpx_count == 2
    
    # Find SOL residues
    sol_count = sum(1 for res in residues if res.get_resname() == "SOL")
    assert sol_count == 10


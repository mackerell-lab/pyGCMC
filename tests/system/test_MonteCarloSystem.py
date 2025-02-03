# tests/system/test_MonteCarloSystem.py

import pytest
import pygcmc
import os
import math

# Get the directory containing test data files
TEST_DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")

@pytest.fixture
def molecular_system():
    """Create a test molecular system from PDB and TOP files."""
    pdb_path = os.path.join(TEST_DATA_DIR, "test.pdb")
    top_path = os.path.join(TEST_DATA_DIR, "test.top")
    
    structure = pygcmc.PDBParser.parse_file(pdb_path)
    topology = pygcmc.TOPParser.parse_file(top_path)
    
    mol_system = pygcmc.MolecularSystem()
    return mol_system.combine(structure, topology)

def assert_arrays_almost_equal(arr1, arr2, tol=1e-6):
    """Compare two arrays for approximate equality."""
    if len(arr1) != len(arr2):
        return False
    return all(abs(a - b) < tol for a, b in zip(arr1, arr2))

def test_initialize_from_molecular(molecular_system):
    """Test conversion from MolecularSystem to MonteCarloSystem."""
    # Create Monte Carlo system
    mc_system = pygcmc.MonteCarloSystem()
    
    # Set reasonable max capacity
    info = pygcmc.GCMCInfo()
    info.max_residues = 1000
    info.max_atoms = 10000
    mc_system.initialize(info)
    
    # Initialize from molecular system
    mc_system.initialize_from_molecular(molecular_system)
    
    # Check basic properties
    assert mc_system.get_active_residue_count() == molecular_system.get_num_residues()
    assert mc_system.get_active_atom_count() == molecular_system.get_num_atoms()
    
    # Check box dimensions
    mc_box = mc_system.get_state().info.box
    mol_box = molecular_system.boxDimensions
    assert_arrays_almost_equal(mc_box, mol_box)
    
    # Check residues
    mc_state = mc_system.get_state()
    for i in range(mc_system.get_active_residue_count()):
        mc_res = mc_state.residues[i]
        mol_res = molecular_system.residues[i]
        
        assert mc_res.atom_count == len(mol_res.get_atoms())
        assert mc_res.active == True
        
        # Check atoms in residue
        for j in range(mc_res.atom_count):
            mc_atom = mc_state.atoms[mc_res.atom_start + j]
            mol_atom = mol_res.get_atoms()[j]
            
            # Check atom properties
            assert mc_system.get_type_maps().get_type_name(mc_atom.type) == mol_atom.get_type()
            assert abs(mc_atom.charge - mol_atom.get_charge()) < 1e-6
            assert_arrays_almost_equal(
                [mc_atom.x, mc_atom.y, mc_atom.z],
                [mol_atom.get_x(), mol_atom.get_y(), mol_atom.get_z()]
            )
        
        # Check center of mass calculation
        com = [0.0, 0.0, 0.0]
        for j in range(mc_res.atom_count):
            atom = mc_state.atoms[mc_res.atom_start + j]
            com[0] += atom.x
            com[1] += atom.y
            com[2] += atom.z
        
        if mc_res.atom_count > 0:
            com = [x / mc_res.atom_count for x in com]
        
        assert_arrays_almost_equal(mc_res.com, com)

def test_initialize_from_molecular_empty():
    """Test conversion with empty molecular system."""
    mol_system = pygcmc.MolecularSystem()
    mc_system = pygcmc.MonteCarloSystem()
    
    # Set reasonable max capacity
    info = pygcmc.GCMCInfo()
    info.max_residues = 1000
    info.max_atoms = 10000
    mc_system.initialize(info)
    
    # Should raise an exception
    with pytest.raises(RuntimeError, match="MolecularSystem has no molecular data"):
        mc_system.initialize_from_molecular(mol_system)

def test_initialize_from_molecular_large_system(molecular_system):
    """Test conversion with a system that exceeds max capacity."""
    mc_system = pygcmc.MonteCarloSystem()
    
    # Set small max capacity
    info = pygcmc.GCMCInfo()
    info.max_residues = 1
    info.max_atoms = 1
    mc_system.initialize(info)
    
    # Should raise an exception
    with pytest.raises(RuntimeError) as excinfo:
        mc_system.initialize_from_molecular(molecular_system)
    assert "Initial system exceeds max capacity" in str(excinfo.value)

def test_type_mapping():
    """Test the atom type mapping system."""
    mc_system = pygcmc.MonteCarloSystem()
    
    # Test adding new types
    type1_idx = mc_system.get_type_maps().get_or_add_type("C")
    type2_idx = mc_system.get_type_maps().get_or_add_type("O")
    type3_idx = mc_system.get_type_maps().get_or_add_type("N")
    
    # Test retrieving existing types
    assert mc_system.get_type_maps().get_or_add_type("C") == type1_idx
    assert mc_system.get_type_maps().get_or_add_type("O") == type2_idx
    assert mc_system.get_type_maps().get_or_add_type("N") == type3_idx
    
    # Test getting type names
    assert mc_system.get_type_maps().get_type_name(type1_idx) == "C"
    assert mc_system.get_type_maps().get_type_name(type2_idx) == "O"
    assert mc_system.get_type_maps().get_type_name(type3_idx) == "N"
    
    # Test invalid index
    assert mc_system.get_type_maps().get_type_name(-1) == ""
    assert mc_system.get_type_maps().get_type_name(1000) == ""


# tests/system/MonteCarloSystem/mapping_tests.py

import pytest
from .helpers import *

def test_residue_atom_properties():
    """Test detailed properties of residues and atoms in MonteCarloSystem."""
    # Create Monte Carlo system and initialize it
    mc_system = pygcmc.MonteCarloSystem()
    info = pygcmc.MCInfo()
    info.max_residues = 1000
    info.max_atoms = 10000
    mc_system.initialize(info)
    
    # Get test data paths
    pdb_path = os.path.join(TEST_DATA_DIR, "test.pdb")
    top_path = os.path.join(TEST_DATA_DIR, "test.top")
    
    # Parse PDB and TOP files
    structure = pygcmc.PDBParser.parse_file(pdb_path)
    topology = pygcmc.TOPParser.parse_file(top_path)
    
    # Create molecular system and combine structure with topology
    mol_system = pygcmc.MolecularSystem()
    molecular = mol_system.combine(structure, topology)
    
    # Initialize Monte Carlo system from molecular system
    mc_system.initialize_from_molecular(molecular)
    
    # Get initial state
    state = mc_system.get_state()
    
    # Verify all residues
    expected_atom_start = 0
    for res_idx in range(state.activeResidueCount):
        mc_res = state.residues[res_idx]
        mol_res = molecular.residues[res_idx]
        top_res = molecular.topology_residues[res_idx]
        mol_atoms = mol_res.get_atoms()
        
        # Verify residue properties
        assert mc_res.active == True, f"Residue {res_idx} should be active"
        assert mc_res.atom_count == len(mol_atoms), \
            f"Residue {res_idx} has incorrect atom count: {mc_res.atom_count} != {len(mol_atoms)}"
        assert mc_res.atom_start == expected_atom_start, \
            f"Residue {res_idx} has incorrect atom_start: {mc_res.atom_start} != {expected_atom_start}"
        
        # Store and verify atom properties
        for atom_idx in range(mc_res.atom_count):
            mc_atom = state.atoms[mc_res.atom_start + atom_idx]
            mol_atom = mol_atoms[atom_idx]
            top_atom = molecular.topology_atoms[top_res.atoms[atom_idx]]
            
            # Verify atom type
            atom_type_name = mc_system.get_type_maps().get_type_name(mc_atom.type)
            assert atom_type_name == top_atom.type, \
                f"Atom {atom_idx} in residue {res_idx} has incorrect type: {atom_type_name} != {top_atom.type}"
            
            # Verify atom charge
            assert abs(mc_atom.charge - top_atom.charge) < 1e-6, \
                f"Atom {atom_idx} in residue {res_idx} has incorrect charge"
            
            # Verify atom coordinates
            assert_arrays_almost_equal(
                [mc_atom.x, mc_atom.y, mc_atom.z],
                [mol_atom.get_x(), mol_atom.get_y(), mol_atom.get_z()]
            ), f"Atom {atom_idx} in residue {res_idx} has incorrect coordinates"
        
        # Calculate and verify center of mass
        com = [0.0, 0.0, 0.0]
        for atom_idx in range(mc_res.atom_count):
            atom = state.atoms[mc_res.atom_start + atom_idx]
            com[0] += atom.x
            com[1] += atom.y
            com[2] += atom.z
        
        if mc_res.atom_count > 0:
            com = [x / mc_res.atom_count for x in com]
        
        assert_arrays_almost_equal(mc_res.center, com)
        
        # Update expected_atom_start for next residue
        expected_atom_start += mc_res.atom_count
    
    # Verify total atom count
    assert expected_atom_start == state.activeAtomCount, \
        f"Total atom count mismatch: {expected_atom_start} != {state.activeAtomCount}"

def test_residue_atom_properties_empty():
    """Test residue and atom properties with empty molecular system."""
    mol_system = pygcmc.MolecularSystem()
    mc_system = pygcmc.MonteCarloSystem()
    
    # Set reasonable max capacity
    info = pygcmc.MCInfo()
    info.max_residues = 1000
    info.max_atoms = 10000
    mc_system.initialize(info)
    
    # Should raise an exception
    with pytest.raises(RuntimeError, match="MolecularSystem has no molecular data"):
        mc_system.initialize_from_molecular(mol_system)

def test_type_mapping_from_molecular():
    """Test that type mapping is correctly created when initializing from molecular system."""
    # Create Monte Carlo system and initialize it
    mc_system = pygcmc.MonteCarloSystem()
    info = pygcmc.MCInfo()
    info.max_residues = 1000
    info.max_atoms = 10000
    mc_system.initialize(info)
    
    # Get test data paths
    pdb_path = os.path.join(TEST_DATA_DIR, "test.pdb")
    top_path = os.path.join(TEST_DATA_DIR, "test.top")
    
    # Parse PDB and TOP files
    structure = pygcmc.PDBParser.parse_file(pdb_path)
    topology = pygcmc.TOPParser.parse_file(top_path)
    
    # Create molecular system and combine structure with topology
    mol_system = pygcmc.MolecularSystem()
    molecular = mol_system.combine(structure, topology)
    
    # Initialize Monte Carlo system from molecular system
    mc_system.initialize_from_molecular(molecular)
    
    # Get state and type maps
    state = mc_system.get_state()
    type_maps = mc_system.get_type_maps()
    
    # Create a set of all unique atom types in molecular system
    mol_types = set()
    for res in molecular.topology_residues:
        for atom_idx in res.atoms:
            mol_types.add(molecular.topology_atoms[atom_idx].type)
    
    # Create a set of all unique atom types in monte carlo system
    mc_types = set()
    for i in range(state.activeAtomCount):
        atom = state.atoms[i]
        type_name = type_maps.get_type_name(atom.type)
        mc_types.add(type_name)
    
    # Verify that both systems have the same atom types
    assert mol_types == mc_types, f"Type mismatch: molecular types {mol_types} != monte carlo types {mc_types}"
    
    # Verify that each atom's type is correctly mapped
    for res_idx in range(state.activeResidueCount):
        mc_res = state.residues[res_idx]
        top_res = molecular.topology_residues[res_idx]
        
        for atom_idx in range(mc_res.atom_count):
            mc_atom = state.atoms[mc_res.atom_start + atom_idx]
            top_atom = molecular.topology_atoms[top_res.atoms[atom_idx]]
            
            mc_type = type_maps.get_type_name(mc_atom.type)
            top_type = top_atom.type
            
            assert mc_type == top_type, \
                f"Type mismatch in residue {res_idx}, atom {atom_idx}: {mc_type} != {top_type}"
            
            # Verify that the type index is consistent
            assert mc_atom.type == type_maps.get_or_add_type(top_type), \
                f"Type index mismatch in residue {res_idx}, atom {atom_idx}"
    
    # Verify that type indices are continuous and start from 0
    type_indices = set()
    for atom_type in mol_types:
        idx = type_maps.get_or_add_type(atom_type)
        type_indices.add(idx)
    
    assert len(type_indices) == len(mol_types), "Number of type indices doesn't match number of types"
    assert min(type_indices) == 0, "Type indices should start from 0"
    assert max(type_indices) == len(mol_types) - 1, "Type indices should be continuous"

def test_residue_type_mapping_from_molecular():
    """Test that residue type mapping is correctly created when initializing from molecular system."""
    # Create Monte Carlo system and initialize it
    mc_system = pygcmc.MonteCarloSystem()
    info = pygcmc.MCInfo()
    info.max_residues = 1000
    info.max_atoms = 10000
    mc_system.initialize(info)
    
    # Get test data paths
    pdb_path = os.path.join(TEST_DATA_DIR, "test.pdb")
    top_path = os.path.join(TEST_DATA_DIR, "test.top")
    
    # Parse PDB and TOP files
    structure = pygcmc.PDBParser.parse_file(pdb_path)
    topology = pygcmc.TOPParser.parse_file(top_path)
    
    # Create molecular system and combine structure with topology
    mol_system = pygcmc.MolecularSystem()
    molecular = mol_system.combine(structure, topology)
    
    # Initialize Monte Carlo system from molecular system
    mc_system.initialize_from_molecular(molecular)
    
    # Get state and type maps
    state = mc_system.get_state()
    residue_type_maps = state.residueTypes
    
    # Create a set of all unique residue types in molecular system
    mol_res_types = set()
    for residue in molecular.residues:
        mol_res_types.add(residue.get_resname())
    
    # Create a set of all unique residue types in monte carlo system
    mc_res_types = set()
    for i in range(state.activeResidueCount):
        residue = state.residues[i]
        type_name = residue_type_maps.get_type_name(residue.type)
        mc_res_types.add(type_name)
    
    # Verify that both systems have the same residue types
    assert mol_res_types == mc_res_types, \
        f"Type mismatch: molecular residue types {mol_res_types} != monte carlo residue types {mc_res_types}"
    
    # Verify that each residue's type is correctly mapped
    for res_idx in range(state.activeResidueCount):
        mc_res = state.residues[res_idx]
        mol_res = molecular.residues[res_idx]
        
        mc_type = residue_type_maps.get_type_name(mc_res.type)
        mol_type = mol_res.get_resname()
        
        assert mc_type == mol_type, \
            f"Type mismatch for residue {res_idx}: {mc_type} != {mol_type}"
        
        # Verify that the type index is consistent
        assert mc_res.type == residue_type_maps.get_or_add_type(mol_type), \
            f"Type index mismatch for residue {res_idx}"
    
    # Verify that type indices are continuous and start from 0
    type_indices = set()
    for res_type in mol_res_types:
        idx = residue_type_maps.get_or_add_type(res_type)
        type_indices.add(idx)
    
    assert len(type_indices) == len(mol_res_types), "Number of residue type indices doesn't match number of types"
    assert min(type_indices) == 0, "Residue type indices should start from 0"
    assert max(type_indices) == len(mol_res_types) - 1, "Residue type indices should be continuous"


# tests/system/test_MonteCarloSystem.py

import pytest
import pygcmc
import os
import math
from pygcmc.model import Molecular
from pygcmc.io import PDBParser, TOPParser

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
    info = pygcmc.MCInfo()
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
        top_res = molecular_system.topology_residues[i]
        
        assert mc_res.atom_count == len(mol_res.get_atoms())
        assert mc_res.active == True
        
        # Check atoms in residue
        for j in range(mc_res.atom_count):
            mc_atom = mc_state.atoms[mc_res.atom_start + j]
            mol_atom = mol_res.get_atoms()[j]
            top_atom = molecular_system.topology_atoms[top_res.atoms[j]]
            
            # Check atom properties
            assert mc_system.get_type_maps().get_type_name(mc_atom.type) == top_atom.type
            assert abs(mc_atom.charge - top_atom.charge) < 1e-6
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
        
        assert_arrays_almost_equal(mc_res.center, com)

def test_initialize_from_molecular_empty():
    """Test conversion with empty molecular system."""
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

def test_initialize_from_molecular_large_system(molecular_system):
    """Test conversion with a system that exceeds max capacity."""
    mc_system = pygcmc.MonteCarloSystem()
    
    # Set small max capacity
    info = pygcmc.MCInfo()
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
        expected_com = [0.0, 0.0, 0.0]
        for atom_idx in range(mc_res.atom_count):
            atom = state.atoms[mc_res.atom_start + atom_idx]
            expected_com[0] += atom.x
            expected_com[1] += atom.y
            expected_com[2] += atom.z
        
        if mc_res.atom_count > 0:
            expected_com = [x / mc_res.atom_count for x in expected_com]
        
        assert_arrays_almost_equal(mc_res.center, expected_com), \
            f"Residue {res_idx} has incorrect geometric center"
        
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

def test_add_movement_molecules():
    """Test adding benx movement molecules to MonteCarloSystem."""
    # Create Monte Carlo system and initialize it
    mc_system = pygcmc.MonteCarloSystem()
    info = pygcmc.MCInfo()
    info.max_residues = 10000  # Large enough for base system + movement molecules
    info.max_atoms = 100000
    mc_system.initialize(info)
    
    # Load base system
    pdb_path = os.path.join(TEST_DATA_DIR, "test.pdb")
    top_path = os.path.join(TEST_DATA_DIR, "test.top")
    structure = pygcmc.PDBParser.parse_file(pdb_path)
    topology = pygcmc.TOPParser.parse_file(top_path)
    mol_system = pygcmc.MolecularSystem()
    base_molecular = mol_system.combine(structure, topology)
    
    # Initialize base system
    mc_system.initialize_from_molecular(base_molecular)
    
    # Store initial state
    initial_res_count = mc_system.get_active_residue_count()
    initial_atom_count = mc_system.get_active_atom_count()
    
    # Load benx molecule
    MOLS_DIR = os.path.join(TEST_DATA_DIR, "mols")
    benx_pdb = os.path.join(MOLS_DIR, "benx.pdb")
    benx_psf = os.path.join(MOLS_DIR, "benx.psf")
    
    # Verify test files exist
    for test_file in [benx_pdb, benx_psf]:
        if not os.path.exists(test_file):
            pytest.skip(f"Test file not found: {test_file}")
    
    # Create molecular system for benx
    benx_structure = pygcmc.PDBParser.parse_file(benx_pdb)
    benx_topology = pygcmc.PSFParser.parse_file(benx_psf)
    benx_molecular = pygcmc.MolecularSystem().combine(benx_structure, benx_topology)
    
    # Get expected atom types from BENX
    benx_types = set()
    for res in benx_molecular.topology_residues:
        for atom_idx in res.atoms:
            benx_types.add(benx_molecular.topology_atoms[atom_idx].type)
    print("\nBENX atom types:", sorted(list(benx_types)))
    
    # Create movement molecule info list with only benx
    movement_mols = [
        pygcmc.MovementMolecularInfo(benx_molecular, 20)
    ]
    
    # Add movement molecules
    mc_system.add_movement_molecules(movement_mols)
    
    # Get final state
    state = mc_system.get_state()
    
    # Test 1: Check movement residue info
    assert len(state.movementResidues) == 1, "Should have info for benx movement molecule"
    
    # Verify movement atom types
    print("\nMovement atom types:", sorted(list(state.movementAtomTypes)))
    assert state.numMovementAtomTypes == len(benx_types), \
           f"Wrong number of movement atom types: {state.numMovementAtomTypes} != {len(benx_types)}"
    assert set(state.atomTypes.get_type_name(t) for t in state.movementAtomTypes) == benx_types, \
           "Movement atom types do not match BENX types"
    
    # Verify benx movement info
    benx_info = state.movementResidues[0]
    benx_name = benx_molecular.residues[0].get_resname()
    
    # Print debug info
    print("\nBenx Movement Info:")
    print(f"Expected name: {benx_name}")
    print(f"Actual name: {benx_info.resName}")
    print(f"Start index: {benx_info.startIndex}")
    print(f"Active count: {benx_info.activeCount}")
    print(f"Total count: {benx_info.totalCount}")
    
    # Count actual benx residues in initial state
    initial_benx_count = sum(1 for i in range(initial_res_count)
                           if state.residueTypes.get_type_name(state.residues[i].type) == benx_name)
    
    print(f"\nInitial benx count: {initial_benx_count}")
    
    assert benx_info.resName == benx_name, f"Wrong residue name: {benx_info.resName} != {benx_name}"
    assert benx_info.activeCount == initial_benx_count, \
           f"Active benx count mismatch: {benx_info.activeCount} != {initial_benx_count}"
    assert benx_info.totalCount == initial_benx_count + 20, \
           f"Total benx count should be active + 20"
    
    # Test 2: Verify residue organization
    # Check that active benx residues are contiguous and at the correct position
    for i in range(benx_info.startIndex, benx_info.startIndex + benx_info.activeCount):
        res = state.residues[i]
        assert res.active == True, f"Residue at {i} should be active"
        assert res.fixed == False, f"Residue at {i} should not be fixed"
        assert state.residueTypes.get_type_name(res.type) == benx_name, \
               f"Residue at {i} has wrong type"
    
    # Check inactive benx residues
    for i in range(benx_info.startIndex + benx_info.activeCount, 
                  benx_info.startIndex + benx_info.totalCount):
        res = state.residues[i]
        assert res.active == False, f"Residue at {i} should be inactive"
        assert res.fixed == False, f"Residue at {i} should not be fixed"
        assert state.residueTypes.get_type_name(res.type) == benx_name, \
               f"Residue at {i} has wrong type"
    
    # Test 3: Verify atom organization
    # Check that atoms for each residue are contiguous
    for i in range(state.activeResidueCount):
        res = state.residues[i]
        # Verify atom indices are within bounds
        assert res.atomStart >= 0 and res.atomStart + res.atomCount <= len(state.atoms), \
               f"Residue {i} has invalid atom range"
        
        if i > 0:
            prev_res = state.residues[i-1]
            # Verify atoms are contiguous
            assert res.atomStart == prev_res.atomStart + prev_res.atomCount, \
                   f"Gap in atom indices between residues {i-1} and {i}"
    
    # Test 4: Verify total counts
    final_active_res = mc_system.get_active_residue_count()
    final_active_atoms = mc_system.get_active_atom_count()
    
    expected_new_res = 20  # We added 20 inactive benx residues
    expected_new_atoms = 20 * benx_molecular.residues[0].atom_count()
    
    print(f"\nFinal counts:")
    print(f"Active residues: {final_active_res} (initial: {initial_res_count})")
    print(f"Active atoms: {final_active_atoms} (initial: {initial_atom_count})")
    print(f"Expected new residues: {expected_new_res}")
    print(f"Expected new atoms: {expected_new_atoms}")
    
    assert final_active_res == initial_res_count + expected_new_res, \
           "Final residue count mismatch"
    assert final_active_atoms == initial_atom_count + expected_new_atoms, \
           "Final atom count mismatch"

def test_add_two_movement_molecules():
    """Test adding both benx and sol movement molecules to MonteCarloSystem."""
    # Create Monte Carlo system and initialize it
    mc_system = pygcmc.MonteCarloSystem()
    info = pygcmc.MCInfo()
    info.max_residues = 10000  # Large enough for base system + movement molecules
    info.max_atoms = 100000
    mc_system.initialize(info)
    
    # Load base system
    pdb_path = os.path.join(TEST_DATA_DIR, "test.pdb")
    top_path = os.path.join(TEST_DATA_DIR, "test.top")
    structure = pygcmc.PDBParser.parse_file(pdb_path)
    topology = pygcmc.TOPParser.parse_file(top_path)
    mol_system = pygcmc.MolecularSystem()
    base_molecular = mol_system.combine(structure, topology)
    
    # Initialize base system
    mc_system.initialize_from_molecular(base_molecular)
    
    # Store initial state
    initial_res_count = mc_system.get_active_residue_count()
    initial_atom_count = mc_system.get_active_atom_count()
    
    # Load molecules
    MOLS_DIR = os.path.join(TEST_DATA_DIR, "mols")
    benx_pdb = os.path.join(MOLS_DIR, "benx.pdb")
    benx_psf = os.path.join(MOLS_DIR, "benx.psf")
    sol_pdb = os.path.join(MOLS_DIR, "sol.pdb")
    sol_itp = os.path.join(MOLS_DIR, "sol.itp")
    
    # Verify test files exist
    for test_file in [benx_pdb, benx_psf, sol_pdb, sol_itp]:
        if not os.path.exists(test_file):
            pytest.skip(f"Test file not found: {test_file}")
    
    # Create molecular systems
    benx_structure = pygcmc.PDBParser.parse_file(benx_pdb)
    benx_topology = pygcmc.PSFParser.parse_file(benx_psf)
    benx_molecular = pygcmc.MolecularSystem().combine(benx_structure, benx_topology)
    
    sol_structure = pygcmc.PDBParser.parse_file(sol_pdb)
    sol_topology = pygcmc.TOPParser.parse_file(sol_itp)
    sol_molecular = pygcmc.MolecularSystem().combine(sol_structure, sol_topology)
    
    # Get expected atom types from each molecule
    benx_types = set()
    for res in benx_molecular.topology_residues:
        for atom_idx in res.atoms:
            benx_types.add(benx_molecular.topology_atoms[atom_idx].type)
    print("\nBENX atom types:", sorted(list(benx_types)))
    
    sol_types = set()
    for res in sol_molecular.topology_residues:
        for atom_idx in res.atoms:
            sol_types.add(sol_molecular.topology_atoms[atom_idx].type)
    print("\nSOL atom types:", sorted(list(sol_types)))
    
    # Create movement molecule info list
    movement_mols = [
        pygcmc.MovementMolecularInfo(benx_molecular, 20),
        pygcmc.MovementMolecularInfo(sol_molecular, 20)
    ]

    # Get molecule names
    benx_name = benx_molecular.residues[0].get_resname()
    sol_name = sol_molecular.residues[0].get_resname()

    state = mc_system.get_state()

    # Count initial residues of each type
    initial_benx_count = sum(1 for i in range(initial_res_count)
                           if state.residueTypes.get_type_name(state.residues[i].type) == benx_name)
    initial_sol_count = sum(1 for i in range(initial_res_count)
                          if state.residueTypes.get_type_name(state.residues[i].type) == sol_name)
    
    # Add movement molecules
    mc_system.add_movement_molecules(movement_mols)
    
    # Get final state
    state = mc_system.get_state()
    
    # Test 1: Check movement residue info
    assert len(state.movementResidues) == 2, "Should have info for both movement molecules"
    
    # Print debug info for both molecules
    print("\nMovement Info:")
    for info in state.movementResidues:
        print(f"\nMolecule: {info.resName}")
        print(f"Start index: {info.startIndex}")
        print(f"Active count: {info.activeCount}")
        print(f"Total count: {info.totalCount}")
    
    print(f"\nInitial counts:")
    print(f"Benx: {initial_benx_count}")
    print(f"Sol: {initial_sol_count}")
    
    # Verify movement info for each molecule
    benx_info = next(info for info in state.movementResidues if info.resName == benx_name)
    sol_info = next(info for info in state.movementResidues if info.resName == sol_name)
    
    # Verify benx info
    assert benx_info.activeCount == initial_benx_count, \
           f"Active benx count mismatch: {benx_info.activeCount} != {initial_benx_count}"
    assert benx_info.totalCount == initial_benx_count + 20, \
           f"Total benx count should be active + 20"
    
    # Verify sol info
    assert sol_info.activeCount == initial_sol_count, \
           f"Active sol count mismatch: {sol_info.activeCount} != {initial_sol_count}"
    assert sol_info.totalCount == initial_sol_count + 20, \
           f"Total sol count should be active + 20"
    
    # Test 2: Verify residue organization
    # Helper function to check residue range
    def verify_residue_range(info, name):
        # Check active residues
        for i in range(info.startIndex, info.startIndex + info.activeCount):
            res = state.residues[i]
            assert res.active == True, f"{name} residue at {i} should be active"
            assert res.fixed == False, f"{name} residue at {i} should not be fixed"
            assert state.residueTypes.get_type_name(res.type) == name, \
                   f"Residue at {i} should be {name}"
        
        # Check inactive residues
        for i in range(info.startIndex + info.activeCount, 
                      info.startIndex + info.totalCount):
            res = state.residues[i]
            assert res.active == False, f"{name} residue at {i} should be inactive"
            assert res.fixed == False, f"{name} residue at {i} should not be fixed"
            assert state.residueTypes.get_type_name(res.type) == name, \
                   f"Residue at {i} should be {name}"
    
    verify_residue_range(benx_info, benx_name)
    verify_residue_range(sol_info, sol_name)
    
    # Test 3: Verify atom organization
    # Check that atoms for each residue are contiguous
    for i in range(state.activeResidueCount):
        res = state.residues[i]
        # Verify atom indices are within bounds
        assert res.atomStart >= 0 and res.atomStart + res.atomCount <= len(state.atoms), \
               f"Residue {i} has invalid atom range"
        
        if i > 0:
            prev_res = state.residues[i-1]
            # Verify atoms are contiguous
            assert res.atomStart == prev_res.atomStart + prev_res.atomCount, \
                   f"Gap in atom indices between residues {i-1} and {i}"
    
    # Test 4: Verify total counts
    final_active_res = mc_system.get_active_residue_count()
    final_active_atoms = mc_system.get_active_atom_count()
    
    expected_new_res = 40  # 20 inactive benx + 20 inactive sol
    expected_new_atoms = (20 * benx_molecular.residues[0].atom_count() + 
                         20 * sol_molecular.residues[0].atom_count())
    
    print(f"\nFinal counts:")
    print(f"Active residues: {final_active_res} (initial: {initial_res_count})")
    print(f"Active atoms: {final_active_atoms} (initial: {initial_atom_count})")
    print(f"Expected new residues: {expected_new_res}")
    print(f"Expected new atoms: {expected_new_atoms}")
    
    assert final_active_res == initial_res_count + expected_new_res, \
           "Final residue count mismatch"
    assert final_active_atoms == initial_atom_count + expected_new_atoms, \
           "Final atom count mismatch"
    
    # Test 5: Verify movement residues are properly ordered
    # Movement residues should be in the order they were added
    assert benx_info.startIndex < sol_info.startIndex, \
           "Benx section should come before sol section"

def test_add_three_movement_molecules():
    """Test adding sol, imia, and benx movement molecules to MonteCarloSystem."""
    # Create Monte Carlo system and initialize it
    mc_system = pygcmc.MonteCarloSystem()
    info = pygcmc.MCInfo()
    info.max_residues = 10000  # Large enough for base system + movement molecules
    info.max_atoms = 100000
    mc_system.initialize(info)
    
    # Load base system
    pdb_path = os.path.join(TEST_DATA_DIR, "test.pdb")
    top_path = os.path.join(TEST_DATA_DIR, "test.top")
    structure = pygcmc.PDBParser.parse_file(pdb_path)
    topology = pygcmc.TOPParser.parse_file(top_path)
    mol_system = pygcmc.MolecularSystem()
    base_molecular = mol_system.combine(structure, topology)
    
    # Initialize base system
    mc_system.initialize_from_molecular(base_molecular)
    
    # Store initial state
    initial_res_count = mc_system.get_active_residue_count()
    initial_atom_count = mc_system.get_active_atom_count()
    
    # Load molecules
    MOLS_DIR = os.path.join(TEST_DATA_DIR, "mols")
    benx_pdb = os.path.join(MOLS_DIR, "benx.pdb")
    benx_psf = os.path.join(MOLS_DIR, "benx.psf")
    sol_pdb = os.path.join(MOLS_DIR, "sol.pdb")
    sol_itp = os.path.join(MOLS_DIR, "sol.itp")
    imia_pdb = os.path.join(MOLS_DIR, "imia.pdb")
    imia_psf = os.path.join(MOLS_DIR, "imia.psf")
    
    # Verify test files exist
    for test_file in [benx_pdb, benx_psf, sol_pdb, sol_itp, imia_pdb, imia_psf]:
        if not os.path.exists(test_file):
            pytest.skip(f"Test file not found: {test_file}")
    
    # Create molecular systems
    benx_structure = pygcmc.PDBParser.parse_file(benx_pdb)
    benx_topology = pygcmc.PSFParser.parse_file(benx_psf)
    benx_molecular = pygcmc.MolecularSystem().combine(benx_structure, benx_topology)
    
    sol_structure = pygcmc.PDBParser.parse_file(sol_pdb)
    sol_topology = pygcmc.TOPParser.parse_file(sol_itp)
    sol_molecular = pygcmc.MolecularSystem().combine(sol_structure, sol_topology)
    
    imia_structure = pygcmc.PDBParser.parse_file(imia_pdb)
    imia_topology = pygcmc.PSFParser.parse_file(imia_psf)
    imia_molecular = pygcmc.MolecularSystem().combine(imia_structure, imia_topology)
    
    # Get expected atom types from each molecule
    benx_types = set()
    for res in benx_molecular.topology_residues:
        for atom_idx in res.atoms:
            benx_types.add(benx_molecular.topology_atoms[atom_idx].type)
    print("\nBENX atom types:", sorted(list(benx_types)))
    
    sol_types = set()
    for res in sol_molecular.topology_residues:
        for atom_idx in res.atoms:
            sol_types.add(sol_molecular.topology_atoms[atom_idx].type)
    print("\nSOL atom types:", sorted(list(sol_types)))
    
    imia_types = set()
    for res in imia_molecular.topology_residues:
        for atom_idx in res.atoms:
            imia_types.add(imia_molecular.topology_atoms[atom_idx].type)
    print("\nIMIA atom types:", sorted(list(imia_types)))
    
    # Create movement molecule info list
    movement_mols = [
        pygcmc.MovementMolecularInfo(sol_molecular, 20),
        pygcmc.MovementMolecularInfo(imia_molecular, 20),
        pygcmc.MovementMolecularInfo(benx_molecular, 20)
    ]

    # Get molecule names
    sol_name = sol_molecular.residues[0].get_resname()
    imia_name = imia_molecular.residues[0].get_resname()
    benx_name = benx_molecular.residues[0].get_resname()

    state = mc_system.get_state()

    # Count initial residues of each type
    initial_sol_count = sum(1 for i in range(initial_res_count)
                          if state.residueTypes.get_type_name(state.residues[i].type) == sol_name)
    initial_imia_count = sum(1 for i in range(initial_res_count)
                           if state.residueTypes.get_type_name(state.residues[i].type) == imia_name)
    initial_benx_count = sum(1 for i in range(initial_res_count)
                           if state.residueTypes.get_type_name(state.residues[i].type) == benx_name)
    
    # Add movement molecules
    mc_system.add_movement_molecules(movement_mols)
    
    # Get final state
    state = mc_system.get_state()
    
    # Test 1: Check movement residue info
    assert len(state.movementResidues) == 3, "Should have info for all three movement molecules"
    
    # Print debug info for all molecules
    print("\nMovement Info:")
    for info in state.movementResidues:
        print(f"\nMolecule: {info.resName}")
        print(f"Start index: {info.startIndex}")
        print(f"Active count: {info.activeCount}")
        print(f"Total count: {info.totalCount}")
    
    print(f"\nInitial counts:")
    print(f"Sol: {initial_sol_count}")
    print(f"Imia: {initial_imia_count}")
    print(f"Benx: {initial_benx_count}")
    
    # Verify movement info for each molecule
    sol_info = next(info for info in state.movementResidues if info.resName == sol_name)
    imia_info = next(info for info in state.movementResidues if info.resName == imia_name)
    benx_info = next(info for info in state.movementResidues if info.resName == benx_name)
    
    # Verify sol info
    assert sol_info.activeCount == initial_sol_count, \
           f"Active sol count mismatch: {sol_info.activeCount} != {initial_sol_count}"
    assert sol_info.totalCount == initial_sol_count + 20, \
           f"Total sol count should be active + 20"
    
    # Verify imia info
    assert imia_info.activeCount == initial_imia_count, \
           f"Active imia count mismatch: {imia_info.activeCount} != {initial_imia_count}"
    assert imia_info.totalCount == initial_imia_count + 20, \
           f"Total imia count should be active + 20"
    
    # Verify benx info
    assert benx_info.activeCount == initial_benx_count, \
           f"Active benx count mismatch: {benx_info.activeCount} != {initial_benx_count}"
    assert benx_info.totalCount == initial_benx_count + 20, \
           f"Total benx count should be active + 20"
    
    # Test 2: Verify residue organization
    # Helper function to check residue range
    def verify_residue_range(info, name):
        # Check active residues
        for i in range(info.startIndex, info.startIndex + info.activeCount):
            res = state.residues[i]
            assert res.active == True, f"{name} residue at {i} should be active"
            assert res.fixed == False, f"{name} residue at {i} should not be fixed"
            assert state.residueTypes.get_type_name(res.type) == name, \
                   f"Residue at {i} should be {name}"
        
        # Check inactive residues
        for i in range(info.startIndex + info.activeCount, 
                      info.startIndex + info.totalCount):
            res = state.residues[i]
            assert res.active == False, f"{name} residue at {i} should be inactive"
            assert res.fixed == False, f"{name} residue at {i} should not be fixed"
            assert state.residueTypes.get_type_name(res.type) == name, \
                   f"Residue at {i} should be {name}"
    
    verify_residue_range(sol_info, sol_name)
    verify_residue_range(imia_info, imia_name)
    verify_residue_range(benx_info, benx_name)
    
    # Test 3: Verify atom organization
    # Check that atoms for each residue are contiguous
    for i in range(state.activeResidueCount):
        res = state.residues[i]
        # Verify atom indices are within bounds
        assert res.atomStart >= 0 and res.atomStart + res.atomCount <= len(state.atoms), \
               f"Residue {i} has invalid atom range"
        
        if i > 0:
            prev_res = state.residues[i-1]
            # Verify atoms are contiguous
            assert res.atomStart == prev_res.atomStart + prev_res.atomCount, \
                   f"Gap in atom indices between residues {i-1} and {i}"
    
    # Test 4: Verify total counts
    final_active_res = mc_system.get_active_residue_count()
    final_active_atoms = mc_system.get_active_atom_count()
    
    expected_new_res = 60  # 20 inactive each for sol, imia, and benx
    expected_new_atoms = (20 * sol_molecular.residues[0].atom_count() + 
                         20 * imia_molecular.residues[0].atom_count() +
                         20 * benx_molecular.residues[0].atom_count())
    
    print(f"\nFinal counts:")
    print(f"Active residues: {final_active_res} (initial: {initial_res_count})")
    print(f"Active atoms: {final_active_atoms} (initial: {initial_atom_count})")
    print(f"Expected new residues: {expected_new_res}")
    print(f"Expected new atoms: {expected_new_atoms}")
    
    assert final_active_res == initial_res_count + expected_new_res, \
           "Final residue count mismatch"
    assert final_active_atoms == initial_atom_count + expected_new_atoms, \
           "Final atom count mismatch"
    
    # Test 5: Verify movement residues are properly ordered
    # Movement residues should be in the order they were added
    assert sol_info.startIndex < imia_info.startIndex < benx_info.startIndex, \
           "Movement sections should be in order: sol, imia, benx"

def test_atom_type_maps_after_movement_molecules():
    """Test that atom type maps are correctly maintained after adding movement molecules."""
    # Create Monte Carlo system and initialize it
    mc_system = pygcmc.MonteCarloSystem()
    info = pygcmc.MCInfo()
    info.max_residues = 10000  # Large enough for base system + movement molecules
    info.max_atoms = 100000
    mc_system.initialize(info)
    
    # Load base system
    pdb_path = os.path.join(TEST_DATA_DIR, "test.pdb")
    top_path = os.path.join(TEST_DATA_DIR, "test.top")
    structure = pygcmc.PDBParser.parse_file(pdb_path)
    topology = pygcmc.TOPParser.parse_file(top_path)
    mol_system = pygcmc.MolecularSystem()
    base_molecular = mol_system.combine(structure, topology)
    
    # Initialize base system
    mc_system.initialize_from_molecular(base_molecular)
    
    # Store initial type maps
    initial_type_maps = mc_system.get_type_maps()
    initial_types = set(initial_type_maps.atomTypes)
    print("\nInitial atom types:", sorted(list(initial_types)))
    
    # Load molecules
    MOLS_DIR = os.path.join(TEST_DATA_DIR, "mols")
    benx_pdb = os.path.join(MOLS_DIR, "benx.pdb")
    benx_psf = os.path.join(MOLS_DIR, "benx.psf")
    sol_pdb = os.path.join(MOLS_DIR, "sol.pdb")
    sol_itp = os.path.join(MOLS_DIR, "sol.itp")
    imia_pdb = os.path.join(MOLS_DIR, "imia.pdb")
    imia_psf = os.path.join(MOLS_DIR, "imia.psf")
    
    # Verify test files exist
    for test_file in [benx_pdb, benx_psf, sol_pdb, sol_itp, imia_pdb, imia_psf]:
        if not os.path.exists(test_file):
            pytest.skip(f"Test file not found: {test_file}")
    
    # Create molecular systems
    benx_structure = pygcmc.PDBParser.parse_file(benx_pdb)
    benx_topology = pygcmc.PSFParser.parse_file(benx_psf)
    benx_molecular = pygcmc.MolecularSystem().combine(benx_structure, benx_topology)
    
    sol_structure = pygcmc.PDBParser.parse_file(sol_pdb)
    sol_topology = pygcmc.TOPParser.parse_file(sol_itp)
    sol_molecular = pygcmc.MolecularSystem().combine(sol_structure, sol_topology)
    
    imia_structure = pygcmc.PDBParser.parse_file(imia_pdb)
    imia_topology = pygcmc.PSFParser.parse_file(imia_psf)
    imia_molecular = pygcmc.MolecularSystem().combine(imia_structure, imia_topology)
    
    # Get expected atom types from each molecule
    benx_types = set()
    for res in benx_molecular.topology_residues:
        for atom_idx in res.atoms:
            benx_types.add(benx_molecular.topology_atoms[atom_idx].type)
    print("\nBENX atom types:", sorted(list(benx_types)))
    
    sol_types = set()
    for res in sol_molecular.topology_residues:
        for atom_idx in res.atoms:
            sol_types.add(sol_molecular.topology_atoms[atom_idx].type)
    print("\nSOL atom types:", sorted(list(sol_types)))
    
    imia_types = set()
    for res in imia_molecular.topology_residues:
        for atom_idx in res.atoms:
            imia_types.add(imia_molecular.topology_atoms[atom_idx].type)
    print("\nIMIA atom types:", sorted(list(imia_types)))
    
    # Create movement molecule info list
    movement_mols = [
        pygcmc.MovementMolecularInfo(sol_molecular, 20),
        pygcmc.MovementMolecularInfo(imia_molecular, 20),
        pygcmc.MovementMolecularInfo(benx_molecular, 20)
    ]
    
    # Add movement molecules
    mc_system.add_movement_molecules(movement_mols)
    
    # Get final type maps
    final_type_maps = mc_system.get_type_maps()
    final_types = set(final_type_maps.atomTypes)
    print("\nFinal atom types:", sorted(list(final_types)))
    
    # Verify all atom types are present
    expected_types = initial_types | benx_types | sol_types | imia_types
    assert final_types == expected_types, \
        f"Type mismatch after adding movement molecules:\nExpected: {sorted(list(expected_types))}\nGot: {sorted(list(final_types))}"
    
    # Verify type indices are preserved
    for atom_type in initial_types:
        initial_idx = initial_type_maps.get_or_add_type(atom_type)
        final_idx = final_type_maps.get_or_add_type(atom_type)
        assert initial_idx == final_idx, \
            f"Type index changed for {atom_type}: {initial_idx} -> {final_idx}"
    
    # Verify type indices are continuous
    type_indices = set()
    for atom_type in final_types:
        idx = final_type_maps.get_or_add_type(atom_type)
        type_indices.add(idx)
    
    assert len(type_indices) == len(final_types), "Number of type indices doesn't match number of types"
    assert min(type_indices) == 0, "Type indices should start from 0"
    assert max(type_indices) == len(final_types) - 1, "Type indices should be continuous"

def test_compare_add_molecules_together_vs_separate():
    """Test that adding molecules together vs separately produces equivalent results."""
    # Create two Monte Carlo systems for comparison
    mc_system1 = pygcmc.MonteCarloSystem()
    mc_system2 = pygcmc.MonteCarloSystem()
    info = pygcmc.MCInfo()
    info.max_residues = 10000
    info.max_atoms = 100000
    mc_system1.initialize(info)
    mc_system2.initialize(info)

    # Load base system into both
    pdb_path = os.path.join(TEST_DATA_DIR, "test.pdb")
    top_path = os.path.join(TEST_DATA_DIR, "test.top")
    structure = pygcmc.PDBParser.parse_file(pdb_path)
    topology = pygcmc.TOPParser.parse_file(top_path)
    mol_system = pygcmc.MolecularSystem()
    base_molecular = mol_system.combine(structure, topology)
    
    mc_system1.initialize_from_molecular(base_molecular)
    mc_system2.initialize_from_molecular(base_molecular)

    # Load molecules
    MOLS_DIR = os.path.join(TEST_DATA_DIR, "mols")
    benx_pdb = os.path.join(MOLS_DIR, "benx.pdb")
    benx_psf = os.path.join(MOLS_DIR, "benx.psf")
    imia_pdb = os.path.join(MOLS_DIR, "imia.pdb")
    imia_psf = os.path.join(MOLS_DIR, "imia.psf")

    # Verify test files exist
    for test_file in [benx_pdb, benx_psf, imia_pdb, imia_psf]:
        if not os.path.exists(test_file):
            pytest.skip(f"Test file not found: {test_file}")

    # Create molecular systems
    benx_structure = pygcmc.PDBParser.parse_file(benx_pdb)
    benx_topology = pygcmc.PSFParser.parse_file(benx_psf)
    benx_molecular = pygcmc.MolecularSystem().combine(benx_structure, benx_topology)

    imia_structure = pygcmc.PDBParser.parse_file(imia_pdb)
    imia_topology = pygcmc.PSFParser.parse_file(imia_psf)
    imia_molecular = pygcmc.MolecularSystem().combine(imia_structure, imia_topology)

    # Get expected atom types from each molecule
    benx_types = set()
    for res in benx_molecular.topology_residues:
        for atom_idx in res.atoms:
            benx_types.add(benx_molecular.topology_atoms[atom_idx].type)
    print("\nBENX atom types:", sorted(list(benx_types)))

    imia_types = set()
    for res in imia_molecular.topology_residues:
        for atom_idx in res.atoms:
            imia_types.add(imia_molecular.topology_atoms[atom_idx].type)
    print("\nIMIA atom types:", sorted(list(imia_types)))

    # System 1: Add both molecules together
    movement_mols = [
        pygcmc.MovementMolecularInfo(benx_molecular, 20),
        pygcmc.MovementMolecularInfo(imia_molecular, 20)
    ]
    mc_system1.add_movement_molecules(movement_mols)

    # System 2: Add molecules separately
    mc_system2.add_movement_molecules([pygcmc.MovementMolecularInfo(benx_molecular, 20)])
    mc_system2.add_movement_molecules([pygcmc.MovementMolecularInfo(imia_molecular, 20)])

    # Get final states
    state1 = mc_system1.get_state()
    state2 = mc_system2.get_state()

    # Compare movement atom types
    print("\nSystem 1 movement atom types:", sorted(list(state1.movementAtomTypes)))
    print("System 2 movement atom types:", sorted(list(state2.movementAtomTypes)))
    
    expected_types = benx_types | imia_types
    assert state1.numMovementAtomTypes == len(expected_types), \
           f"System 1: Wrong number of movement atom types: {state1.numMovementAtomTypes} != {len(expected_types)}"
    assert state2.numMovementAtomTypes == len(expected_types), \
           f"System 2: Wrong number of movement atom types: {state2.numMovementAtomTypes} != {len(expected_types)}"
    
    assert set(state1.atomTypes.get_type_name(t) for t in state1.movementAtomTypes) == expected_types, \
           "System 1: Movement atom types do not match expected types"
    assert set(state2.atomTypes.get_type_name(t) for t in state2.movementAtomTypes) == expected_types, \
           "System 2: Movement atom types do not match expected types"

    # Test 1: Compare basic counts and box info
    assert state1.activeResidueCount == state2.activeResidueCount, \
           "Active residue counts differ"
    assert state1.activeAtomCount == state2.activeAtomCount, \
           "Active atom counts differ"
    assert len(state1.movementResidues) == len(state2.movementResidues), \
           "Number of movement residues differ"
    assert_arrays_almost_equal(state1.info.box, state2.info.box), \
           "Box dimensions differ"

    # Test 2: Compare type maps
    # 2.1: Compare residue type maps
    residue_types1 = set(state1.residueTypes.atomTypes)
    residue_types2 = set(state2.residueTypes.atomTypes)
    print("\nResidue Types:")
    print("System 1:", sorted(list(residue_types1)))
    print("System 2:", sorted(list(residue_types2)))
    assert residue_types1 == residue_types2, \
           f"Residue type sets differ: {residue_types1} != {residue_types2}"

    # 2.2: Compare atom type maps
    atom_types1 = set(mc_system1.get_type_maps().atomTypes)
    atom_types2 = set(mc_system2.get_type_maps().atomTypes)
    print("\nAtom Types:")
    print("System 1:", sorted(list(atom_types1)))
    print("System 2:", sorted(list(atom_types2)))
    assert atom_types1 == atom_types2, \
           f"Atom type sets differ: {atom_types1} != {atom_types2}"

    # Test 3: Compare movement residues info
    # Sort by residue name to ensure consistent comparison
    move_res1 = sorted(state1.movementResidues, key=lambda x: x.resName)
    move_res2 = sorted(state2.movementResidues, key=lambda x: x.resName)
    
    print("\nMovement Residues Comparison:")
    for info1, info2 in zip(move_res1, move_res2):
        print(f"\nComparing {info1.resName}:")
        print(f"System 1: start={info1.startIndex}, active={info1.activeCount}, total={info1.totalCount}")
        print(f"System 2: start={info2.startIndex}, active={info2.activeCount}, total={info2.totalCount}")
        assert info1.resName == info2.resName, \
               f"Movement residue names differ: {info1.resName} != {info2.resName}"
        assert info1.activeCount == info2.activeCount, \
               f"Active counts differ for {info1.resName}: {info1.activeCount} != {info2.activeCount}"
        assert info1.totalCount == info2.totalCount, \
               f"Total counts differ for {info1.resName}: {info1.totalCount} != {info2.totalCount}"

    # Test 4: Compare residue organization and properties
    print("\nResidue Organization Comparison:")
    for i in range(state1.activeResidueCount):
        res1 = state1.residues[i]
        res2 = state2.residues[i]
        name1 = state1.residueTypes.get_type_name(res1.type)
        name2 = state2.residueTypes.get_type_name(res2.type)
        
        print(f"\nResidue {i}:")
        print(f"System 1: type={name1}, active={res1.active}, atoms={res1.atomCount}, start={res1.atomStart}")
        print(f"System 2: type={name2}, active={res2.active}, atoms={res2.atomCount}, start={res2.atomStart}")
        
        assert name1 == name2, f"Residue at {i} has different types: {name1} != {name2}"
        assert res1.active == res2.active, \
               f"Residue at {i} has different active states: {res1.active} != {res2.active}"
        assert res1.atomCount == res2.atomCount, \
               f"Residue at {i} has different atom counts: {res1.atomCount} != {res2.atomCount}"
        assert res1.atomStart == res2.atomStart, \
               f"Residue at {i} has different atom start positions: {res1.atomStart} != {res2.atomStart}"
        assert_arrays_almost_equal(res1.center, res2.center), \
               f"Residue at {i} has different centers"

    # Test 5: Compare atoms in detail
    print("\nAtom Comparison:")
    for i in range(state1.activeAtomCount):
        atom1 = state1.atoms[i]
        atom2 = state2.atoms[i]
        type1 = mc_system1.get_type_maps().get_type_name(atom1.type)
        type2 = mc_system2.get_type_maps().get_type_name(atom2.type)
        
        print(f"\nAtom {i}:")
        print(f"System 1: type={type1}, charge={atom1.charge:.6f}, pos=({atom1.x:.3f}, {atom1.y:.3f}, {atom1.z:.3f})")
        print(f"System 2: type={type2}, charge={atom2.charge:.6f}, pos=({atom2.x:.3f}, {atom2.y:.3f}, {atom2.z:.3f})")
        
        assert type1 == type2, f"Atom at {i} has different types: {type1} != {type2}"
        assert abs(atom1.charge - atom2.charge) < 1e-6, \
               f"Atom at {i} has different charges: {atom1.charge} != {atom2.charge}"
        assert abs(atom1.x - atom2.x) < 1e-6 and \
               abs(atom1.y - atom2.y) < 1e-6 and \
               abs(atom1.z - atom2.z) < 1e-6, \
               f"Atom at {i} has different coordinates"

    # Test 6: Compare type indices mapping
    print("\nType Indices Mapping Comparison:")
    # 6.1: Residue type indices
    for res_type in residue_types1:
        idx1 = state1.residueTypes.get_or_add_type(res_type)
        idx2 = state2.residueTypes.get_or_add_type(res_type)
        print(f"Residue type '{res_type}': System1 index={idx1}, System2 index={idx2}")
        assert state1.residueTypes.get_type_name(idx1) == state2.residueTypes.get_type_name(idx2), \
               f"Residue type name mismatch for indices {idx1} and {idx2}"

    # 6.2: Atom type indices
    for atom_type in atom_types1:
        idx1 = mc_system1.get_type_maps().get_or_add_type(atom_type)
        idx2 = mc_system2.get_type_maps().get_or_add_type(atom_type)
        print(f"Atom type '{atom_type}': System1 index={idx1}, System2 index={idx2}")
        assert mc_system1.get_type_maps().get_type_name(idx1) == mc_system2.get_type_maps().get_type_name(idx2), \
               f"Atom type name mismatch for indices {idx1} and {idx2}"

@pytest.fixture
def charmm_ff():
    # Get the test data directory
    test_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    data_dir = os.path.join(test_dir, "data")
    
    # Load CHARMM force field files
    ff = pygcmc.ForceField()
    pygcmc.PRMParser.parse_file_to_forcefield(os.path.join(data_dir, "par_all36_cgenff.prm"), ff)
    pygcmc.PRMParser.parse_file_to_forcefield(os.path.join(data_dir, "par_all36m_prot.prm"), ff)
    pygcmc.PRMParser.parse_file_to_forcefield(os.path.join(data_dir, "silcs.str"), ff)
    pygcmc.PRMParser.parse_file_to_forcefield(os.path.join(data_dir, "toppar_water_ions.str"), ff)
    return ff

def test_initialize_force_field(molecular_system, charmm_ff):
    print("\n=== Starting test_initialize_force_field ===")
    
    # Create Monte Carlo system
    mc_system = pygcmc.MonteCarloSystem()
    print("Created MonteCarloSystem")
    
    # Set reasonable max capacity
    info = pygcmc.MCInfo()
    info.max_residues = 1000
    info.max_atoms = 10000
    mc_system.initialize(info)
    print("Initialized MonteCarloSystem with max capacity")
    
    # Print molecular system info
    print(f"\nMolecular system info:")
    print(f"Number of residues: {molecular_system.get_num_residues()}")
    print(f"Number of atoms: {molecular_system.get_num_atoms()}")
    print(f"Box dimensions: {molecular_system.boxDimensions}")
    
    # Initialize with molecular system
    mc_system.initialize_from_molecular(molecular_system)
    print("\nInitialized from molecular system")
    
    # Add movement molecules
    test_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    data_dir = os.path.join(test_dir, "data")
    
    # Load small molecules
    movement_molecules = []
    
    # Load benzene
    benx_pdb = os.path.join(data_dir, "mols", "benx.pdb")
    benx_psf = os.path.join(data_dir, "mols", "benx.psf")
    benx_structure = pygcmc.PDBParser.parse_file(benx_pdb)
    benx_topology = pygcmc.PSFParser.parse_file(benx_psf)
    benx_mol = pygcmc.MolecularSystem()
    benx = benx_mol.combine(benx_structure, benx_topology)
    movement_molecules.append(pygcmc.MovementMolecularInfo(benx, 10))  # 10 copies
    print("\nLoaded benzene as movement molecule")
    
    # Load imidazole
    imia_pdb = os.path.join(data_dir, "mols", "imia.pdb")
    imia_psf = os.path.join(data_dir, "mols", "imia.psf")
    imia_structure = pygcmc.PDBParser.parse_file(imia_pdb)
    imia_topology = pygcmc.PSFParser.parse_file(imia_psf)
    imia_mol = pygcmc.MolecularSystem()
    imia = imia_mol.combine(imia_structure, imia_topology)
    movement_molecules.append(pygcmc.MovementMolecularInfo(imia, 10))  # 10 copies
    print("Loaded imidazole as movement molecule")
    
    # Load water
    sol_pdb = os.path.join(data_dir, "mols", "sol.pdb")
    sol_top = os.path.join(data_dir, "mols", "sol.itp")
    sol_structure = pygcmc.PDBParser.parse_file(sol_pdb)
    sol_topology = pygcmc.TOPParser.parse_file(sol_top)
    sol_mol = pygcmc.MolecularSystem()
    sol = sol_mol.combine(sol_structure, sol_topology)
    movement_molecules.append(pygcmc.MovementMolecularInfo(sol, 20))  # 20 copies for water
    print("Loaded water as movement molecule")
    
    # Add movement molecules to system
    mc_system.add_movement_molecules(movement_molecules)
    print("\nAdded all movement molecules to system")
    
    # Get atom types before force field initialization
    atom_types = mc_system.get_type_maps()
    print(f"\nNumber of atom types: {len(atom_types.atomTypes)}")
    print("Atom types in the system:", atom_types.atomTypes)
    
    # Get state and check movement atom types
    state = mc_system.get_state()
    print("\nMovement atom types info:")
    print(f"Number of movement atom types: {state.numMovementAtomTypes}")
    print("Movement atom types indices:", state.movementAtomTypes)
    print("Movement residues info:")
    for info in state.movementResidues:
        print(f"  {info.resName}: start={info.startIndex}, active={info.activeCount}, total={info.totalCount}")
    
    # Print available LJ parameters in the force field
    print("\nAvailable LJ parameters in force field:")
    lj_params = charmm_ff.lj_params
    print(f"Number of LJ parameters: {len(lj_params)}")
    for atom_type, params in lj_params.items():
        print(f"{atom_type}: epsilon={params.epsilon:.4f}, rmin={params.rmin:.4f}")
    
    # Check if all atom types have LJ parameters
    missing_types = []
    for type_name in atom_types.atomTypes:
        try:
            params = charmm_ff.get_lj_params(type_name)
            print(f"Found LJ params for {type_name}: epsilon={params.epsilon:.4f}, rmin={params.rmin:.4f}")
        except Exception as e:
            print(f"Failed to get LJ params for {type_name}: {str(e)}")
            missing_types.append(type_name)
    
    if missing_types:
        print("\nMissing LJ parameters for atom types:", missing_types)
        pytest.skip(f"Missing LJ parameters for atom types: {missing_types}")
    
    print("\nAll atom types have LJ parameters, proceeding with force field initialization")
    
    # Initialize force field
    try:
        mc_system.initialize_force_field(charmm_ff)
        print("Force field initialization successful")
    except Exception as e:
        print(f"Force field initialization failed: {str(e)}")
        raise
    
    # Get state after force field initialization
    state = mc_system.get_state()
    print("\nGot state after force field initialization")
    
    # Check force field parameters
    print("\nForce field parameters:")
    print(f"maxTypes: {state.forcefield.maxTypes}")
    print(f"numMovementTypes: {state.forcefield.numMovementTypes}")
    print(f"ljSigma size: {len(state.forcefield.ljSigma)}")
    print(f"ljEps size: {len(state.forcefield.ljEps)}")
    
    # Basic assertions
    assert state.forcefield.maxTypes == len(atom_types.atomTypes)
    assert state.forcefield.numMovementTypes == state.numMovementAtomTypes
    assert len(state.forcefield.ljSigma) == state.numMovementAtomTypes * state.forcefield.maxTypes
    assert len(state.forcefield.ljEps) == state.numMovementAtomTypes * state.forcefield.maxTypes
    
    print("\n=== test_initialize_force_field completed successfully ===")


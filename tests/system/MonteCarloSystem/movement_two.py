# tests/system/MonteCarloSystem/movement_two.py

import pytest
from .helpers import *

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


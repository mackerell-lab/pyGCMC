# tests/system/MonteCarloSystem/comparison_tests.py

import pytest
from .helpers import *

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


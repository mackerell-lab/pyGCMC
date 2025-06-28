# tests/system/MonteCarloSystem/molecule_addition_comparison.py

import pytest
from .helpers import *

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


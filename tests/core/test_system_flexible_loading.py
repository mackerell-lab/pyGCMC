# tests/core/test_system_flexible_loading.py

import os
import pytest
from pygcmc import System

@pytest.fixture
def test_data_dir():
    """Get the path to the test data directory."""
    current_dir = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(current_dir, "..", "data")

def test_flexible_loading(test_data_dir):
    """Test various combinations of file loading."""
    # Prepare file paths
    pdb_file = os.path.join(test_data_dir, "test.pdb")
    psf_files = [
        os.path.join(test_data_dir, "test_proa.psf"),
        os.path.join(test_data_dir, "mols", "benx.psf"),
        os.path.join(test_data_dir, "mols", "prpx.psf")
    ]
    itp_file = os.path.join(test_data_dir, "mols", "sol.itp")
    top_file = os.path.join(test_data_dir, "test.top")

    # Test different loading combinations
    systems = {
        "PSF only": System(pdb=pdb_file, psf=psf_files[0]),
        "Multiple PSF": System(pdb=pdb_file, psf=psf_files),
        "PSF + ITP": System(pdb=pdb_file, psf=psf_files, itp=itp_file),
        "TOP only": System(pdb=pdb_file, top=top_file)
    }

    # List of residues to check
    residues_to_check = ["ALA", "VAL", "PRO", "ASN", "GLN", "BENX", "PRPX", "TIP3"]

    # Compare all systems against the TOP system (reference)
    reference_system = systems["TOP only"]

    for sys_name, system in systems.items():
        if sys_name == "TOP only":
            continue

        print(f"\nComparing {sys_name} with TOP system:")
        for residue_name in residues_to_check:
            # Get atoms for both systems
            ref_atoms = reference_system.get_pdb_atoms_by_residue(residue_name)
            test_atoms = system.get_pdb_atoms_by_residue(residue_name)

            # Sort atoms by sequence and name
            ref_atoms = sorted(ref_atoms, key=lambda x: (x.sequence, x.name))
            test_atoms = sorted(test_atoms, key=lambda x: (x.sequence, x.name))

            # Print atom counts for debugging
            print(f"\n{residue_name}:")
            print(f"  Reference system has {len(ref_atoms)} atoms")
            print(f"  Test system has {len(test_atoms)} atoms")

            # Skip if residue not present in either system
            if not ref_atoms and not test_atoms:
                continue

            # Verify same number of atoms
            assert len(ref_atoms) == len(test_atoms), \
                f"Different number of atoms for {residue_name}: ref={len(ref_atoms)}, test={len(test_atoms)}"

            # Compare atoms
            for atom_ref, atom_test in zip(ref_atoms, test_atoms):
                print(f"\nComparing atom {atom_ref.name}:")
                print(f"  Reference: type={atom_ref.topo_type}, charge={atom_ref.topo_charge:.3f}, mass={atom_ref.topo_mass:.3f}")
                print(f"  Test: type={atom_test.topo_type}, charge={atom_test.topo_charge:.3f}, mass={atom_test.topo_mass:.3f}")

                # Check atom properties match
                assert atom_ref.name == atom_test.name, \
                    f"Atom name mismatch in {residue_name} {atom_ref.sequence}: {atom_ref.name} vs {atom_test.name}"
                
                assert atom_ref.sequence == atom_test.sequence, \
                    f"Sequence number mismatch for {residue_name} {atom_ref.name}: {atom_ref.sequence} vs {atom_test.sequence}"
                
                # Skip topology checks if either atom has no topology type
                if not atom_ref.topo_type or not atom_test.topo_type:
                    print(f"Warning: Missing topology type for {atom_ref.name}")
                    continue

                assert atom_ref.topo_type == atom_test.topo_type, \
                    f"Topology type mismatch for {residue_name} {atom_ref.sequence} {atom_ref.name}: {atom_ref.topo_type} vs {atom_test.topo_type}"
                
                assert pytest.approx(atom_ref.topo_charge) == atom_test.topo_charge, \
                    f"Charge mismatch for {residue_name} {atom_ref.sequence} {atom_ref.name}: {atom_ref.topo_charge} vs {atom_test.topo_charge}"
                
                assert pytest.approx(atom_ref.topo_mass) == atom_test.topo_mass, \
                    f"Mass mismatch for {residue_name} {atom_ref.sequence} {atom_ref.name}: {atom_ref.topo_mass} vs {atom_test.topo_mass}"

                # Also check coordinates
                assert pytest.approx(atom_ref.x) == atom_test.x, \
                    f"X coordinate mismatch for {residue_name} {atom_ref.sequence} {atom_ref.name}"
                assert pytest.approx(atom_ref.y) == atom_test.y, \
                    f"Y coordinate mismatch for {residue_name} {atom_ref.sequence} {atom_ref.name}"
                assert pytest.approx(atom_ref.z) == atom_test.z, \
                    f"Z coordinate mismatch for {residue_name} {atom_ref.sequence} {atom_ref.name}"

def test_load_structure_flexible(test_data_dir):
    """Test flexible loading using load_structure method."""
    # Prepare file paths
    pdb_file = os.path.join(test_data_dir, "test.pdb")
    psf_files = [
        os.path.join(test_data_dir, "test_proa.psf"),
        os.path.join(test_data_dir, "mols", "benx.psf"),
        os.path.join(test_data_dir, "mols", "prpx.psf"),
        os.path.join(test_data_dir, "mols", "sol.psf")
    ]

    # Test different loading methods
    print("\nTesting different loading methods:")
    
    # Method 1: Using constructor (known to work)
    print("\nMethod 1: Using constructor")
    system1 = System(pdb=pdb_file, psf=psf_files)
    atoms1 = system1.get_pdb_atoms_by_residue("ALA")
    print("\nConstructor method results:")
    for atom in atoms1:
        print(f"Constructor - Atom {atom.name} in ALA: topo_type={atom.topo_type}, charge={atom.topo_charge}, mass={atom.topo_mass}")
    
    # Method 2: Using load_structure with kwargs
    print("\nMethod 2: Using load_structure with kwargs")
    system2 = System()
    # Try both ways of kwargs
    print("\nTrying dict-style kwargs:")
    system2.load_structure(**{"pdb": pdb_file, "psf": psf_files})
    atoms2 = system2.get_pdb_atoms_by_residue("ALA")
    for atom in atoms2:
        print(f"load_structure - Atom {atom.name} in ALA: topo_type={atom.topo_type}, charge={atom.topo_charge}, mass={atom.topo_mass}")

    # Verify system has atoms
    assert system2.has_pdb_atoms(), "System should have PDB atoms after loading"

    # Get some residues to check
    for residue_name in ["ALA", "SOL"]:
        atoms = system2.get_pdb_atoms_by_residue(residue_name)
        assert len(atoms) > 0, f"Should have atoms for residue {residue_name}"
        for atom in atoms:
            assert atom.topo_type, f"Atom {atom.name} in {residue_name} should have topology type"
            assert not pytest.approx(atom.topo_charge) == 0, f"Atom {atom.name} in {residue_name} should have non-zero charge"
            assert not pytest.approx(atom.topo_mass) == 0, f"Atom {atom.name} in {residue_name} should have non-zero mass" 
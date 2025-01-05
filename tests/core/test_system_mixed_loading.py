import os
import pytest
from pygcmc import System, Structure

@pytest.fixture
def test_data_dir():
    """Get the path to the test data directory."""
    current_dir = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(current_dir, "..", "data")

def test_mixed_loading(test_data_dir):
    """Test loading structure using a combination of PSF and ITP files."""
    # Prepare file paths
    pdb_file = os.path.join(test_data_dir, "test.pdb")
    
    # Method 1: All PSF files
    psf_files_all = [
        os.path.join(test_data_dir, "test_proa.psf"),
        os.path.join(test_data_dir, "mols", "benx.psf"),
        os.path.join(test_data_dir, "mols", "prpx.psf"),
        os.path.join(test_data_dir, "mols", "sol.psf")
    ]
    
    # Method 2: Mixed PSF and ITP files
    psf_files_mixed = [
        os.path.join(test_data_dir, "test_proa.psf"),
        os.path.join(test_data_dir, "mols", "benx.psf"),
        os.path.join(test_data_dir, "mols", "prpx.psf")
    ]
    itp_file = os.path.join(test_data_dir, "mols", "sol.itp")

    # Load structures using different methods
    system_all_psf = System(pdb=pdb_file, psf=psf_files_all)
    system_mixed = System(pdb=pdb_file, psf=psf_files_mixed, itp=itp_file)

    # List of residues to check
    residues_to_check = ["ALA", "VAL", "PRO", "ASN", "GLN", "BENX", "PRPX", "TIP3"]

    # Compare systems
    systems = {
        "All PSF": system_all_psf,
        "Mixed PSF+ITP": system_mixed
    }

    # First verify all systems have the same residues
    for residue_name in residues_to_check:
        atoms_by_system = {}
        for sys_name, system in systems.items():
            # Get all atoms for this residue
            atoms = system.get_pdb_atoms_by_residue(residue_name)
            atoms_by_system[sys_name] = sorted(atoms, key=lambda x: (x.sequence, x.name))
            
            # Print atom counts for debugging
            print(f"\n{sys_name} has {len(atoms)} atoms for residue {residue_name}")
            if atoms:
                print(f"First atom: {atoms[0].name} in sequence {atoms[0].sequence}")
                print(f"Last atom: {atoms[-1].name} in sequence {atoms[-1].sequence}")

        # Skip if no system has this residue
        if all(len(atoms) == 0 for atoms in atoms_by_system.values()):
            continue

        # Verify all systems have the same number of atoms for this residue
        atom_counts = {name: len(atoms) for name, atoms in atoms_by_system.items()}
        assert len(set(atom_counts.values())) == 1, \
            f"Different number of atoms for {residue_name}: {atom_counts}"

        # Compare atoms between systems
        psf_atoms = atoms_by_system["All PSF"]
        mixed_atoms = atoms_by_system["Mixed PSF+ITP"]

        print(f"\nComparing All PSF vs Mixed PSF+ITP for {residue_name}:")
        for atom_psf, atom_mixed in zip(psf_atoms, mixed_atoms):
            print(f"\nComparing atom {atom_psf.name}:")
            print(f"All PSF: type={atom_psf.topo_type}, charge={atom_psf.topo_charge:.3f}, mass={atom_psf.topo_mass:.3f}")
            print(f"Mixed: type={atom_mixed.topo_type}, charge={atom_mixed.topo_charge:.3f}, mass={atom_mixed.topo_mass:.3f}")

            # Check atom properties match
            assert atom_psf.name == atom_mixed.name, \
                f"Atom name mismatch in {residue_name} {atom_psf.sequence}: {atom_psf.name} vs {atom_mixed.name}"
            
            assert atom_psf.sequence == atom_mixed.sequence, \
                f"Sequence number mismatch for {residue_name} {atom_psf.name}: {atom_psf.sequence} vs {atom_mixed.sequence}"
            
            # Skip topology checks if either atom has no topology type
            if not atom_psf.topo_type or not atom_mixed.topo_type:
                print(f"Warning: Missing topology type for {atom_psf.name}")
                continue

            assert atom_psf.topo_type == atom_mixed.topo_type, \
                f"Topology type mismatch for {residue_name} {atom_psf.sequence} {atom_psf.name}: {atom_psf.topo_type} vs {atom_mixed.topo_type}"
            
            assert pytest.approx(atom_psf.topo_charge) == atom_mixed.topo_charge, \
                f"Charge mismatch for {residue_name} {atom_psf.sequence} {atom_psf.name}: {atom_psf.topo_charge} vs {atom_mixed.topo_charge}"
            
            assert pytest.approx(atom_psf.topo_mass) == atom_mixed.topo_mass, \
                f"Mass mismatch for {residue_name} {atom_psf.sequence} {atom_psf.name}: {atom_psf.topo_mass} vs {atom_mixed.topo_mass}"

            # Also check coordinates
            assert pytest.approx(atom_psf.x) == atom_mixed.x, \
                f"X coordinate mismatch for {residue_name} {atom_psf.sequence} {atom_psf.name}"
            assert pytest.approx(atom_psf.y) == atom_mixed.y, \
                f"Y coordinate mismatch for {residue_name} {atom_psf.sequence} {atom_psf.name}"
            assert pytest.approx(atom_psf.z) == atom_mixed.z, \
                f"Z coordinate mismatch for {residue_name} {atom_psf.sequence} {atom_psf.name}" 
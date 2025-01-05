import os
import pytest
from pygcmc import System, Structure

@pytest.fixture
def test_data_dir():
    """Get the path to the test data directory."""
    current_dir = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(current_dir, "..", "data")

def test_compare_psf_itp_loading(test_data_dir):
    """Compare loading structure using multiple PSF files versus PSF+ITP files."""
    # Prepare file paths
    pdb_file = os.path.join(test_data_dir, "test.pdb")
    
    # PSF files method
    psf_files = [
        os.path.join(test_data_dir, "test_proa.psf"),
        os.path.join(test_data_dir, "mols", "benx.psf"),
        os.path.join(test_data_dir, "mols", "prpx.psf"),
        os.path.join(test_data_dir, "mols", "sol.psf")
    ]
    
    # PSF+ITP method
    main_psf = os.path.join(test_data_dir, "test_proa.psf")
    itp_files = [
        os.path.join(test_data_dir, "mols", "benx.itp"),
        os.path.join(test_data_dir, "mols", "prpx.itp"),
        os.path.join(test_data_dir, "mols", "sol.itp")
    ]

    # Method 1: Load with multiple PSF files
    system_psf = System(pdb=pdb_file, psf=psf_files)

    # Method 2: Load with PSF+ITP files
    system_itp = System(pdb=pdb_file, psf=main_psf, itp=itp_files)

    # List of residues to check
    residues_to_check = ["ALA", "VAL", "PRO", "ASN", "GLN", "BENX", "PRPX", "TIP3"]

    # Special case mappings for N-terminal residue atoms (PSF -> TOP)
    n_terminal_mappings = {
        "HT1": "H1",
        "HT2": "H2",
        "HT3": "H3"
    }

    # Compare systems
    systems = {
        "PSF": system_psf,
        "PSF+ITP": system_itp
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

        # Compare atoms between PSF and PSF+ITP systems
        psf_atoms = atoms_by_system["PSF"]
        itp_atoms = atoms_by_system["PSF+ITP"]

        print(f"\nComparing PSF vs PSF+ITP for {residue_name}:")
        for atom_psf, atom_itp in zip(psf_atoms, itp_atoms):
            # Get atom names with N-terminal mapping if needed
            name_psf = atom_psf.name
            name_itp = atom_itp.name
            if residue_name == "ALA" and atom_psf.sequence == 7:
                if name_psf in n_terminal_mappings:
                    name_psf = n_terminal_mappings[name_psf]
                if name_itp in n_terminal_mappings:
                    name_itp = n_terminal_mappings[name_itp]

            print(f"\nComparing atom {name_psf}:")
            print(f"PSF: type={atom_psf.topo_type}, charge={atom_psf.topo_charge:.3f}, mass={atom_psf.topo_mass:.3f}")
            print(f"ITP: type={atom_itp.topo_type}, charge={atom_itp.topo_charge:.3f}, mass={atom_itp.topo_mass:.3f}")

            # Check atom properties match
            assert name_psf == name_itp, \
                f"Atom name mismatch in {residue_name} {atom_psf.sequence}: {name_psf} vs {name_itp}"
            
            assert atom_psf.sequence == atom_itp.sequence, \
                f"Sequence number mismatch for {residue_name} {name_psf}: {atom_psf.sequence} vs {atom_itp.sequence}"
            
            # Skip topology checks if either atom has no topology type
            if not atom_psf.topo_type or not atom_itp.topo_type:
                print(f"Warning: Missing topology type for {name_psf}")
                continue

            assert atom_psf.topo_type == atom_itp.topo_type, \
                f"Topology type mismatch for {residue_name} {atom_psf.sequence} {name_psf}: {atom_psf.topo_type} vs {atom_itp.topo_type}"
            
            assert pytest.approx(atom_psf.topo_charge) == atom_itp.topo_charge, \
                f"Charge mismatch for {residue_name} {atom_psf.sequence} {name_psf}: {atom_psf.topo_charge} vs {atom_itp.topo_charge}"
            
            assert pytest.approx(atom_psf.topo_mass) == atom_itp.topo_mass, \
                f"Mass mismatch for {residue_name} {atom_psf.sequence} {name_psf}: {atom_psf.topo_mass} vs {atom_itp.topo_mass}"

            # Also check coordinates
            assert pytest.approx(atom_psf.x) == atom_itp.x, \
                f"X coordinate mismatch for {residue_name} {atom_psf.sequence} {name_psf}"
            assert pytest.approx(atom_psf.y) == atom_itp.y, \
                f"Y coordinate mismatch for {residue_name} {atom_psf.sequence} {name_psf}"
            assert pytest.approx(atom_psf.z) == atom_itp.z, \
                f"Z coordinate mismatch for {residue_name} {atom_psf.sequence} {name_psf}" 
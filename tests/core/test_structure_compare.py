import os
import pytest
from pygcmc import System, Structure

@pytest.fixture
def test_data_dir():
    """Get the path to the test data directory."""
    current_dir = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(current_dir, "..", "data")

def test_compare_structure_creation(test_data_dir):
    """Compare different methods of creating and loading Structure."""
    # Prepare file paths
    pdb_file = os.path.join(test_data_dir, "test.pdb")
    psf_files = [
        os.path.join(test_data_dir, "test_proa.psf"),
        os.path.join(test_data_dir, "mols", "benx.psf"),
        os.path.join(test_data_dir, "mols", "prpx.psf"),
        os.path.join(test_data_dir, "mols", "sol.psf")
    ]
    top_file = os.path.join(test_data_dir, "test.top")

    # Method 1: Create Structure directly with PSF
    structure_psf = Structure()
    structure_psf.read_pdb(pdb_file)
    for psf_file in psf_files:
        structure_psf.read_psf(psf_file)

    # Method 2: Create Structure directly with TOP
    structure_top = Structure()
    structure_top.read_pdb(pdb_file)
    structure_top.read_top(top_file)

    # List of residues to check
    residues_to_check = ["ALA", "VAL", "PRO", "ASN", "GLN", "BENX", "PRPX", "TIP3"]

    # Special case mappings for N-terminal residue atoms (PSF -> TOP)
    n_terminal_mappings = {
        "HT1": "H1",
        "HT2": "H2",
        "HT3": "H3"
    }

    # Compare structures
    structures = {
        "PSF": structure_psf,
        "TOP": structure_top
    }

    # First verify all structures have the same residues
    for residue_name in residues_to_check:
        residues_by_structure = {}
        for struct_name, structure in structures.items():
            # Get all residues with this name
            residues = [res for res in structure.residues if res.name == residue_name]
            residues_by_structure[struct_name] = residues
            
            # Print residue counts for debugging
            print(f"\n{struct_name} has {len(residues)} residues of type {residue_name}")
            for res in residues:
                print(f"  Residue {res.sequence_number} has {len(res.atoms)} atoms")

        # Skip if no structure has this residue
        if all(len(residues) == 0 for residues in residues_by_structure.values()):
            continue

        # Verify all structures have the same number of residues
        residue_counts = {name: len(residues) for name, residues in residues_by_structure.items()}
        assert len(set(residue_counts.values())) == 1, \
            f"Different number of residues for {residue_name}: {residue_counts}"

        # Compare residues between structures
        psf_residues = residues_by_structure["PSF"]
        top_residues = residues_by_structure["TOP"]

        print(f"\nComparing PSF vs TOP for {residue_name}:")
        for res_psf, res_top in zip(psf_residues, top_residues):
            assert res_psf.sequence_number == res_top.sequence_number, \
                f"Sequence number mismatch: {res_psf.sequence_number} vs {res_top.sequence_number}"
            assert len(res_psf.atoms) == len(res_top.atoms), \
                f"Different number of atoms in residue {res_psf.sequence_number}"

            # Sort atoms by name for comparison
            atoms_psf = sorted(res_psf.atoms, key=lambda x: x.name)
            atoms_top = sorted(res_top.atoms, key=lambda x: x.name)

            for atom_psf, atom_top in zip(atoms_psf, atoms_top):
                # Get atom names with N-terminal mapping if needed
                name_psf = atom_psf.name
                name_top = atom_top.name
                if residue_name == "ALA" and res_psf.sequence_number == 7:
                    if name_psf in n_terminal_mappings:
                        name_psf = n_terminal_mappings[name_psf]
                    if name_top in n_terminal_mappings:
                        name_top = n_terminal_mappings[name_top]

                print(f"\nComparing atom {name_psf}:")
                print(f"PSF: type={atom_psf.topo_type}, charge={atom_psf.topo_charge:.3f}, mass={atom_psf.topo_mass:.3f}")
                print(f"TOP: type={atom_top.topo_type}, charge={atom_top.topo_charge:.3f}, mass={atom_top.topo_mass:.3f}")

                # Check atom properties match
                assert name_psf == name_top, \
                    f"Atom name mismatch in {residue_name} {res_psf.sequence_number}: {name_psf} vs {name_top}"

                # Skip topology checks if either atom has no topology type
                if not atom_psf.topo_type or not atom_top.topo_type:
                    print(f"Warning: Missing topology type for {name_psf}")
                    continue

                assert atom_psf.topo_type == atom_top.topo_type, \
                    f"Topology type mismatch for {residue_name} {res_psf.sequence_number} {name_psf}: {atom_psf.topo_type} vs {atom_top.topo_type}"
                
                assert pytest.approx(atom_psf.topo_charge) == atom_top.topo_charge, \
                    f"Charge mismatch for {residue_name} {res_psf.sequence_number} {name_psf}: {atom_psf.topo_charge} vs {atom_top.topo_charge}"
                
                assert pytest.approx(atom_psf.topo_mass) == atom_top.topo_mass, \
                    f"Mass mismatch for {residue_name} {res_psf.sequence_number} {name_psf}: {atom_psf.topo_mass} vs {atom_top.topo_mass}" 
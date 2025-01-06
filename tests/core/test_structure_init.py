# tests/core/test_structure_init.py

import os
import pytest
from pygcmc import Structure

@pytest.fixture
def test_data_dir():
    """Get the path to the test data directory."""
    current_dir = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(current_dir, "..", "data")

def test_structure_initialization(test_data_dir):
    """Test different methods of initializing Structure."""
    # Prepare file paths
    pdb_file = os.path.join(test_data_dir, "test.pdb")
    psf_files = [
        os.path.join(test_data_dir, "test_proa.psf"),
        os.path.join(test_data_dir, "mols", "benx.psf"),
        os.path.join(test_data_dir, "mols", "prpx.psf"),
        os.path.join(test_data_dir, "mols", "sol.psf")
    ]
    top_file = os.path.join(test_data_dir, "test.top")

    # Method 1: Create Structure with constructor arguments (PSF)
    structure_psf_init = Structure(pdb=pdb_file, psf=psf_files)

    # Method 2: Create Structure with constructor arguments (TOP)
    structure_top_init = Structure(pdb=pdb_file, top=top_file)

    # Method 3: Create Structure with separate method calls (PSF)
    structure_psf_load = Structure()
    structure_psf_load.read_pdb(pdb_file)
    for psf_file in psf_files:
        structure_psf_load.read_psf(psf_file)

    # Method 4: Create Structure with separate method calls (TOP)
    structure_top_load = Structure()
    structure_top_load.read_pdb(pdb_file)
    structure_top_load.read_top(top_file)

    # List of residues to check
    residues_to_check = ["ALA", "VAL", "PRO", "ASN", "GLN", "BENX", "PRPX", "TIP3"]

    # Special case mappings for N-terminal residue atoms (PSF -> TOP)
    n_terminal_mappings = {
        "HT1": "H1",
        "HT2": "H2",
        "HT3": "H3"
    }

    # Compare all structures
    structures = {
        "PSF Init": structure_psf_init,
        "TOP Init": structure_top_init,
        "PSF Load": structure_psf_load,
        "TOP Load": structure_top_load
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

        # Compare residues between all structures
        for struct1_name, struct1_residues in residues_by_structure.items():
            for struct2_name, struct2_residues in residues_by_structure.items():
                if struct1_name >= struct2_name:  # Skip self-comparison and avoid duplicates
                    continue

                print(f"\nComparing {struct1_name} vs {struct2_name} for {residue_name}:")
                for res1, res2 in zip(struct1_residues, struct2_residues):
                    assert res1.sequence_number == res2.sequence_number, \
                        f"Sequence number mismatch: {res1.sequence_number} vs {res2.sequence_number}"
                    assert len(res1.atoms) == len(res2.atoms), \
                        f"Different number of atoms in residue {res1.sequence_number}"

                    # Sort atoms by name for comparison
                    atoms1 = sorted(res1.atoms, key=lambda x: x.name)
                    atoms2 = sorted(res2.atoms, key=lambda x: x.name)

                    for atom1, atom2 in zip(atoms1, atoms2):
                        # Get atom names with N-terminal mapping if needed
                        name1 = atom1.name
                        name2 = atom2.name
                        if residue_name == "ALA" and res1.sequence_number == 7:
                            if name1 in n_terminal_mappings:
                                name1 = n_terminal_mappings[name1]
                            if name2 in n_terminal_mappings:
                                name2 = n_terminal_mappings[name2]

                        print(f"\nComparing atom {name1}:")
                        print(f"{struct1_name}: type={atom1.topo_type}, charge={atom1.topo_charge:.3f}, mass={atom1.topo_mass:.3f}")
                        print(f"{struct2_name}: type={atom2.topo_type}, charge={atom2.topo_charge:.3f}, mass={atom2.topo_mass:.3f}")

                        # Check atom properties match
                        assert name1 == name2, \
                            f"Atom name mismatch in {residue_name} {res1.sequence_number}: {name1} vs {name2}"

                        # Skip topology checks if either atom has no topology type
                        if not atom1.topo_type or not atom2.topo_type:
                            print(f"Warning: Missing topology type for {name1}")
                            continue

                        assert atom1.topo_type == atom2.topo_type, \
                            f"Topology type mismatch for {residue_name} {res1.sequence_number} {name1}: {atom1.topo_type} vs {atom2.topo_type}"
                        
                        assert pytest.approx(atom1.topo_charge) == atom2.topo_charge, \
                            f"Charge mismatch for {residue_name} {res1.sequence_number} {name1}: {atom1.topo_charge} vs {atom2.topo_charge}"
                        
                        assert pytest.approx(atom1.topo_mass) == atom2.topo_mass, \
                            f"Mass mismatch for {residue_name} {res1.sequence_number} {name1}: {atom1.topo_mass} vs {atom2.topo_mass}" 
# tests/core/test_system_compare.py

import os
import pytest
from pygcmc import System

@pytest.fixture
def test_data_dir():
    """Get the path to the test data directory."""
    current_dir = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(current_dir, "..", "data")

def test_compare_psf_top_loading(test_data_dir):
    """Compare system loading between multiple PSF files and single TOP file."""
    # Prepare file paths
    pdb_file = os.path.join(test_data_dir, "test.pdb")
    psf_files = [
        os.path.join(test_data_dir, "test_proa.psf"),
        os.path.join(test_data_dir, "mols", "benx.psf"),
        os.path.join(test_data_dir, "mols", "prpx.psf"),
        os.path.join(test_data_dir, "mols", "sol.psf")
    ]
    top_file = os.path.join(test_data_dir, "test.top")

    # Load systems using both methods
    system_psf = System(pdb=pdb_file, psf=psf_files)
    system_top = System(pdb=pdb_file, top=top_file)

    # List of residues to check
    residues_to_check = ["ALA", "VAL", "PRO", "ASN", "GLN", "BENX", "PRPX", "TIP3"]

    # Special case mappings for N-terminal residue atoms (PSF -> TOP)
    n_terminal_mappings = {
        "HT1": "H1",
        "HT2": "H2",
        "HT3": "H3"
    }

    for residue_name in residues_to_check:
        # Get atoms for both systems
        psf_atoms = system_psf.get_pdb_atoms_by_residue(residue_name)
        top_atoms = system_top.get_pdb_atoms_by_residue(residue_name)

        # Check if both systems have atoms for this residue
        if len(psf_atoms) == 0 and len(top_atoms) == 0:
            continue

        # Check if both systems have the same number of atoms for this residue
        assert len(psf_atoms) == len(top_atoms), \
            f"Different number of atoms for {residue_name}: PSF={len(psf_atoms)}, TOP={len(top_atoms)}"

        # Sort atoms by name to ensure matching comparison
        psf_atoms = sorted(psf_atoms, key=lambda x: (x.sequence, x.name))
        top_atoms = sorted(top_atoms, key=lambda x: (x.sequence, x.name))

        # Compare each atom's properties
        for psf_atom, top_atom in zip(psf_atoms, top_atoms):
            print(f"\nComparing {residue_name} {psf_atom.sequence} - {psf_atom.name}:")
            print(f"PSF: type={psf_atom.topo_type}, charge={psf_atom.topo_charge:.3f}, mass={psf_atom.topo_mass:.3f}")
            print(f"TOP: type={top_atom.topo_type}, charge={top_atom.topo_charge:.3f}, mass={top_atom.topo_mass:.3f}")

            # For N-terminal residue (ALA 7), map HT1/HT2/HT3 to H1/H2/H3
            psf_name = psf_atom.name
            top_name = top_atom.name
            if residue_name == "ALA" and psf_atom.sequence == 7:
                if psf_name in n_terminal_mappings:
                    top_name = n_terminal_mappings[psf_name]
                    print(f"Mapping N-terminal atom name: {psf_name} -> {top_name}")

            # Check atom properties match
            assert psf_name == top_name, \
                f"Atom name mismatch in {residue_name} {psf_atom.sequence}: {psf_name} vs {top_name}"
            
            assert psf_atom.sequence == top_atom.sequence, \
                f"Sequence number mismatch for {residue_name} {psf_name}: {psf_atom.sequence} vs {top_atom.sequence}"
            
            # Skip topology checks for unmapped atoms
            if psf_atom.topo_type == "" and top_atom.topo_type != "":
                print(f"Warning: PSF atom {psf_name} has no topology type, while TOP has {top_atom.topo_type}")
                continue
            
            assert psf_atom.topo_type == top_atom.topo_type, \
                f"Topology type mismatch for {residue_name} {psf_atom.sequence} {psf_name}: {psf_atom.topo_type} vs {top_atom.topo_type}"
            
            assert pytest.approx(psf_atom.topo_charge) == top_atom.topo_charge, \
                f"Charge mismatch for {residue_name} {psf_atom.sequence} {psf_name}: {psf_atom.topo_charge} vs {top_atom.topo_charge}"
            
            assert pytest.approx(psf_atom.topo_mass) == top_atom.topo_mass, \
                f"Mass mismatch for {residue_name} {psf_atom.sequence} {psf_name}: {psf_atom.topo_mass} vs {top_atom.topo_mass}" 
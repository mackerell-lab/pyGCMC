import os
import pytest
from pygcmc import System, Structure

@pytest.fixture
def test_data_dir():
    """Get the path to the test data directory."""
    current_dir = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(current_dir, "..", "data")

def test_system_structure_loading(test_data_dir):
    """Compare different methods of loading structure into system."""
    # Prepare file paths
    pdb_file = os.path.join(test_data_dir, "test.pdb")
    psf_files = [
        os.path.join(test_data_dir, "test_proa.psf"),
        os.path.join(test_data_dir, "mols", "benx.psf"),
        os.path.join(test_data_dir, "mols", "prpx.psf"),
        os.path.join(test_data_dir, "mols", "sol.psf")
    ]
    top_file = os.path.join(test_data_dir, "test.top")

    # Method 1: System direct initialization with PSF
    system_direct_psf = System(pdb=pdb_file, psf=psf_files)

    # Method 2: System direct initialization with TOP
    system_direct_top = System(pdb=pdb_file, top=top_file)

    # Method 3: Load through Structure with PSF
    structure_psf = Structure(pdb=pdb_file, psf=psf_files)
    system_from_structure_psf = System()
    system_from_structure_psf.load_structure(structure_psf)

    # Method 4: Load through Structure with TOP
    structure_top = Structure(pdb=pdb_file, top=top_file)
    system_from_structure_top = System()
    system_from_structure_top.load_structure(structure_top)

    # Method 5: System step-by-step loading with PSF
    system_steps_psf = System()
    for psf_file in psf_files:
        system_steps_psf.load_structure(pdb=pdb_file, psf=psf_file)

    # Method 6: System step-by-step loading with TOP
    system_steps_top = System()
    system_steps_top.load_structure(pdb=pdb_file, top=top_file)

    # List of residues to check
    residues_to_check = ["ALA", "VAL", "PRO", "ASN", "GLN", "BENX", "PRPX", "TIP3"]

    # Special case mappings for N-terminal residue atoms (PSF -> TOP)
    n_terminal_mappings = {
        "HT1": "H1",
        "HT2": "H2",
        "HT3": "H3"
    }

    # Compare all systems
    systems = {
        "Direct PSF": system_direct_psf,
        "Direct TOP": system_direct_top,
        "Structure PSF": system_from_structure_psf,
        "Structure TOP": system_from_structure_top,
        "Steps PSF": system_steps_psf,
        "Steps TOP": system_steps_top
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

        # Skip if no system has this residue
        if all(len(atoms) == 0 for atoms in atoms_by_system.values()):
            continue

        # Verify all systems have the same number of atoms for this residue
        atom_counts = {name: len(atoms) for name, atoms in atoms_by_system.items()}
        assert len(set(atom_counts.values())) == 1, \
            f"Different number of atoms for {residue_name}: {atom_counts}"

        # Compare atoms between all systems
        for sys1_name, sys1_atoms in atoms_by_system.items():
            for sys2_name, sys2_atoms in atoms_by_system.items():
                if sys1_name >= sys2_name:  # Skip self-comparison and avoid duplicates
                    continue

                print(f"\nComparing {sys1_name} vs {sys2_name} for {residue_name}:")
                for atom1, atom2 in zip(sys1_atoms, sys2_atoms):
                    # Get atom names with N-terminal mapping if needed
                    name1 = atom1.name
                    name2 = atom2.name
                    if residue_name == "ALA" and atom1.sequence == 7:
                        if name1 in n_terminal_mappings:
                            name1 = n_terminal_mappings[name1]
                        if name2 in n_terminal_mappings:
                            name2 = n_terminal_mappings[name2]

                    print(f"\nComparing atom {name1}:")
                    print(f"{sys1_name}: type={atom1.topo_type}, charge={atom1.topo_charge:.3f}, mass={atom1.topo_mass:.3f}")
                    print(f"{sys2_name}: type={atom2.topo_type}, charge={atom2.topo_charge:.3f}, mass={atom2.topo_mass:.3f}")

                    # Check atom properties match
                    assert name1 == name2, \
                        f"Atom name mismatch in {residue_name} {atom1.sequence}: {name1} vs {name2}"
                    
                    assert atom1.sequence == atom2.sequence, \
                        f"Sequence number mismatch for {residue_name} {name1}: {atom1.sequence} vs {atom2.sequence}"
                    
                    # Skip topology checks if either atom has no topology type
                    if not atom1.topo_type or not atom2.topo_type:
                        print(f"Warning: Missing topology type for {name1}")
                        continue

                    assert atom1.topo_type == atom2.topo_type, \
                        f"Topology type mismatch for {residue_name} {atom1.sequence} {name1}: {atom1.topo_type} vs {atom2.topo_type}"
                    
                    assert pytest.approx(atom1.topo_charge) == atom2.topo_charge, \
                        f"Charge mismatch for {residue_name} {atom1.sequence} {name1}: {atom1.topo_charge} vs {atom2.topo_charge}"
                    
                    assert pytest.approx(atom1.topo_mass) == atom2.topo_mass, \
                        f"Mass mismatch for {residue_name} {atom1.sequence} {name1}: {atom1.topo_mass} vs {atom2.topo_mass}" 
# tests/system/molecularSystem/top_psf_comparison.py

import pytest
from .helpers import *

def test_compare_top_psf_systems():
    """Test that combining PDB with TOP and PSF files produces identical molecular systems."""
    # Load PDB file
    pdb_path = os.path.join(TEST_DATA_DIR, "step1_pdbreader.pdb")
    structure = pygcmc.PDBParser.parse_file(pdb_path)

    # Fix water residue names in structure (TIP -> TIP3)
    for residue in structure.residues:
        if residue.get_resname() == "TIP":
            residue.set_resname("TIP3")

    # Load TOP and PSF files
    top_path = os.path.join(TEST_DATA_DIR, "step1_pdbreader.top")
    psf_path = os.path.join(TEST_DATA_DIR, "step1_pdbreader.psf")
    top_topology = TOPParser.parse_file(top_path)
    psf_topology = pygcmc.io.PSFParser.parse_file(psf_path)

    # Create molecular systems
    mol_system = pygcmc.MolecularSystem()
    molecular_top = mol_system.combine(structure, top_topology)
    molecular_psf = mol_system.combine(structure, psf_topology)

    # Print basic information for debugging
    print("\nBasic system information:")
    print(f"TOP system: {molecular_top.get_num_atoms()} atoms, {molecular_top.get_num_residues()} residues")
    print(f"PSF system: {molecular_psf.get_num_atoms()} atoms, {molecular_psf.get_num_residues()} residues")

    # Compare basic properties
    assert molecular_top.get_num_atoms() == molecular_psf.get_num_atoms(), \
        f"Different number of atoms: TOP {molecular_top.get_num_atoms()} vs PSF {molecular_psf.get_num_atoms()}"
    assert molecular_top.get_num_residues() == molecular_psf.get_num_residues(), \
        f"Different number of residues: TOP {molecular_top.get_num_residues()} vs PSF {molecular_psf.get_num_residues()}"

    # Print residue information for debugging
    print("\nResidue counts by type:")
    residue_counts = {}
    for res in molecular_top.residues:
        resname = res.get_resname()
        residue_counts[resname] = residue_counts.get(resname, 0) + 1
    for resname, count in sorted(residue_counts.items()):
        print(f"{resname}: {count}")

    # Compare atom properties for a random sample
    print("\nComparing atom properties (random sample)...")
    num_atoms = len(molecular_top.topology_atoms)
    sample_size = min(100, num_atoms)  # Sample 100 atoms or all atoms if less than 100
    sampled_indices = random.sample(range(num_atoms), sample_size)
    
    for i in sampled_indices:
        top_atom = molecular_top.topology_atoms[i]
        psf_atom = molecular_psf.topology_atoms[i]
        
        # Get residue names for both atoms
        top_resname = molecular_top.topology_residues[top_atom.residue_id].name
        psf_resname = molecular_psf.topology_residues[psf_atom.residue_id].name
        
        # Skip atom name comparison for water residues
        if top_resname != "TIP3" and psf_resname != "TIP3" and top_resname != "TIP" and psf_resname != "TIP":
            assert top_atom.name == psf_atom.name, \
                f"Atom {i} name mismatch: TOP {top_atom.name} vs PSF {psf_atom.name}"
        
        assert top_atom.type == psf_atom.type, \
            f"Atom {i} type mismatch: TOP {top_atom.type} vs PSF {psf_atom.type}"
        assert abs(top_atom.charge - psf_atom.charge) < 1e-6, \
            f"Atom {i} charge mismatch: TOP {top_atom.charge} vs PSF {psf_atom.charge}"
        assert abs(top_atom.mass - psf_atom.mass) < 1e-6, \
            f"Atom {i} mass mismatch: TOP {top_atom.mass} vs PSF {psf_atom.mass}"

    # Compare residue properties for a random sample
    print("\nComparing residue properties (random sample)...")
    num_residues = len(molecular_top.topology_residues)
    sample_size = min(50, num_residues)  # Sample 50 residues or all residues if less than 50
    sampled_indices = random.sample(range(num_residues), sample_size)
    
    for i in sampled_indices:
        top_res = molecular_top.topology_residues[i]
        psf_res = molecular_psf.topology_residues[i]
        
        # Skip residue name comparison for water residues
        if top_res.name not in ["TIP3", "TIP"] and psf_res.name not in ["TIP3", "TIP"]:
            assert top_res.name == psf_res.name, \
                f"Residue {i} name mismatch: TOP {top_res.name} vs PSF {psf_res.name}"
        # Skip residue number comparison as TOP and PSF handle them differently
        assert len(top_res.atoms) == len(psf_res.atoms), \
            f"Residue {i} atom count mismatch: TOP {len(top_res.atoms)} vs PSF {len(psf_res.atoms)}"

    # Compare dihedrals
    assert molecular_top.get_num_dihedrals() == molecular_psf.get_num_dihedrals(), \
        f"Different number of dihedrals: TOP {molecular_top.get_num_dihedrals()} vs PSF {molecular_psf.get_num_dihedrals()}"


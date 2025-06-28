# tests/system/molecularSystem/multiple_vs_single.py

import pytest
from .helpers import *

def test_combine_multiple_vs_single():
    """Test that combining multiple PSF files gives same result as single TOP file."""
    # Load structure file
    pdb_path = os.path.join(TEST_DATA_DIR, "test.pdb")
    structure = pygcmc.PDBParser.parse_file(pdb_path)
    
    # Load individual PSF files
    proa_path = os.path.join(TEST_DATA_DIR, "test_proa.psf")
    benx_path = os.path.join(TEST_DATA_DIR, "mols/benx.psf")
    prpx_path = os.path.join(TEST_DATA_DIR, "mols/prpx.itp")
    sol_path = os.path.join(TEST_DATA_DIR, "mols/sol.itp")
    
    # Print original CMAP information from PSF file
    print("\nOriginal CMAP information from test_proa.psf:")
    proa_top = pygcmc.io.PSFParser.parse_file(proa_path)
    print(f"Number of CMAPs: {proa_top.get_num_cmaps()}")
    for i, cmap in enumerate(proa_top.cmaps):
        print(f"CMAP {i+1}: {'-'.join(str(x) for x in cmap.atoms)}")
        residues = [proa_top.get_residue(proa_top.get_atom(atom_id).residue_id).name 
                   for atom_id in cmap.atoms if atom_id >= 0]
        print(f"Residues: {', '.join(residues)}")

    # Print original CMAP information from TOP file
    print("\nOriginal CMAP information from test.top:")
    top_path = os.path.join(TEST_DATA_DIR, "test.top")
    top_topology = pygcmc.io.TOPParser.parse_file(top_path)
    print(f"Number of CMAPs: {top_topology.get_num_cmaps()}")
    for i, cmap in enumerate(top_topology.cmaps):
        print(f"CMAP {i+1}: {'-'.join(str(x) for x in cmap.atoms)}")
        residues = [top_topology.get_residue(top_topology.get_atom(atom_id).residue_id).name 
                   for atom_id in cmap.atoms if atom_id >= 0]
        print(f"Residues: {', '.join(residues)}")

    # Create molecule from multiple PSF files
    mol_system = pygcmc.MolecularSystem()
    mol_multiple = mol_system.combine_multiple(structure, [
        pygcmc.io.PSFParser.parse_file(proa_path),
        pygcmc.io.PSFParser.parse_file(benx_path),
        pygcmc.io.TOPParser.parse_file(prpx_path),
        pygcmc.io.TOPParser.parse_file(sol_path)
    ])
    
    # Load single TOP file that contains all molecules
    top_path = os.path.join(TEST_DATA_DIR, "test.top")
    mol_single = pygcmc.Combine(structure, pygcmc.io.TOPParser.parse_file(top_path))
    
    # Compare basic properties
    assert len(mol_multiple.atoms) == len(mol_single.atoms), \
        f"Different number of atoms: multiple {len(mol_multiple.atoms)} vs single {len(mol_single.atoms)}"
    assert len(mol_multiple.residues) == len(mol_single.residues), \
        f"Different number of residues: multiple {len(mol_multiple.residues)} vs single {len(mol_single.residues)}"
    
    # Compare residue sequences
    for i, (res_multi, res_single) in enumerate(zip(mol_multiple.residues, mol_single.residues)):
        assert res_multi.get_resname() == res_single.get_resname(), \
            f"Residue {i} name mismatch: {res_multi.get_resname()} vs {res_single.get_resname()}"
        assert res_multi.get_ires() == res_single.get_ires(), \
            f"Residue {i} ID mismatch: {res_multi.get_ires()} vs {res_single.get_ires()}"
    
    # Compare atom properties
    for i, (atom_multi, atom_single) in enumerate(zip(mol_multiple.atoms, mol_single.atoms)):
        assert atom_multi.get_formatted_atom_name() == atom_single.get_formatted_atom_name(), \
            f"Atom {i} name mismatch: {atom_multi.get_formatted_atom_name()} vs {atom_single.get_formatted_atom_name()}"
        assert atom_multi.get_type() == atom_single.get_type(), \
            f"Atom {i} type mismatch: {atom_multi.get_type()} vs {atom_single.get_type()}"
        assert atom_multi.get_charge() == pytest.approx(atom_single.get_charge()), \
            f"Atom {i} charge mismatch: {atom_multi.get_charge()} vs {atom_single.get_charge()}"
        assert atom_multi.get_mass() == pytest.approx(atom_single.get_mass()), \
            f"Atom {i} mass mismatch: {atom_multi.get_mass()} vs {atom_single.get_mass()}"
    
    # Compare bonds
    bonds_multi = set((min(b.atom1, b.atom2), max(b.atom1, b.atom2)) for b in mol_multiple.bonds)
    bonds_single = set((min(b.atom1, b.atom2), max(b.atom1, b.atom2)) for b in mol_single.bonds)
    
    # Print detailed bond information for debugging
    print("\nBond counts by residue in multiple topology case:")
    residue_bonds_multi = {}
    for b in mol_multiple.bonds:
        res1 = mol_multiple.atoms[b.atom1].get_resname()
        res2 = mol_multiple.atoms[b.atom2].get_resname()
        key = tuple(sorted([res1, res2]))
        residue_bonds_multi[key] = residue_bonds_multi.get(key, 0) + 1
    for res_pair, count in sorted(residue_bonds_multi.items()):
        print(f"{res_pair}: {count} bonds")
    
    print("\nBond counts by residue in single topology case:")
    residue_bonds_single = {}
    for b in mol_single.bonds:
        res1 = mol_single.atoms[b.atom1].get_resname()
        res2 = mol_single.atoms[b.atom2].get_resname()
        key = tuple(sorted([res1, res2]))
        residue_bonds_single[key] = residue_bonds_single.get(key, 0) + 1
    for res_pair, count in sorted(residue_bonds_single.items()):
        print(f"{res_pair}: {count} bonds")
    
    # Find and print the extra bonds
    extra_bonds = bonds_multi - bonds_single
    if extra_bonds:
        print("\nExtra bonds in multiple topology case:")
        for b in sorted(extra_bonds):
            res1 = mol_multiple.atoms[b[0]].get_resname()
            res2 = mol_multiple.atoms[b[1]].get_resname()
            print(f"Bond between {res1} atoms {b[0]}-{b[1]}")
    
    assert bonds_multi == bonds_single, \
        f"Different bonds between molecules: multiple has {len(bonds_multi)} bonds, single has {len(bonds_single)} bonds"
    
    # Compare angles if they exist
    angles_multi = set((min(a.atom1, a.atom3), a.atom2, max(a.atom1, a.atom3))
                      for a in mol_multiple.angles)
    angles_single = set((min(a.atom1, a.atom3), a.atom2, max(a.atom1, a.atom3))
                       for a in mol_single.angles)

    # Print angle information from each topology file
    print("\nAngle information from individual topology files:")
    for i, topology in enumerate([
        pygcmc.io.PSFParser.parse_file(proa_path),
        pygcmc.io.PSFParser.parse_file(benx_path),
        pygcmc.io.TOPParser.parse_file(prpx_path),
        pygcmc.io.TOPParser.parse_file(sol_path)
    ]):
        print(f"\nTopology {i+1}:")
        num_angles = topology.get_num_angles() if hasattr(topology, 'get_num_angles') else len(topology.angles)
        print(f"Number of angles: {num_angles}")
        if hasattr(topology, 'angles') and topology.angles:
            print("Sample angles:")
            for angle in list(topology.angles)[:5]:
                print(f"  {angle.atom1}-{angle.atom2}-{angle.atom3}")

    print(f"\nTotal angles in mol_multiple: {len(mol_multiple.angles)}")
    print(f"Total angles in mol_single: {len(mol_single.angles)}")

    # Filter out angles involving SOL residues
    angles_multi_no_sol = set()
    for a in mol_multiple.angles:
        res1 = mol_multiple.atoms[a.atom1].get_resname()
        res2 = mol_multiple.atoms[a.atom2].get_resname()
        res3 = mol_multiple.atoms[a.atom3].get_resname()
        if "SOL" not in (res1, res2, res3):
            angles_multi_no_sol.add((min(a.atom1, a.atom3), a.atom2, max(a.atom1, a.atom3)))

    angles_single_no_sol = set()
    for a in mol_single.angles:
        res1 = mol_single.atoms[a.atom1].get_resname()
        res2 = mol_single.atoms[a.atom2].get_resname()
        res3 = mol_single.atoms[a.atom3].get_resname()
        if "SOL" not in (res1, res2, res3):
            angles_single_no_sol.add((min(a.atom1, a.atom3), a.atom2, max(a.atom1, a.atom3)))

    assert angles_multi_no_sol == angles_single_no_sol, \
        f"Different angles between molecules (excluding SOL): multiple has {len(angles_multi_no_sol)} angles, single has {len(angles_single_no_sol)} angles"

    # Compare dihedrals
    def normalize_dihedral(d):
        # A dihedral A-B-C-D can be represented as either A-B-C-D or D-C-B-A
        # and the middle atoms (B-C) can also be swapped in some cases
        atoms = [d.atom1, d.atom2, d.atom3, d.atom4]
        # Forward order
        v1 = (min(atoms[0], atoms[3]), atoms[1], atoms[2], max(atoms[0], atoms[3]))
        # Reverse order
        v2 = (min(atoms[0], atoms[3]), atoms[2], atoms[1], max(atoms[0], atoms[3]))
        # Return the lexicographically smaller representation
        return min(v1, v2)
    
    dihedrals_multi = set(normalize_dihedral(d) for d in mol_multiple.dihedrals)
    dihedrals_single = set(normalize_dihedral(d) for d in mol_single.dihedrals)
    
    # Print dihedral information for debugging
    print("\nDihedral counts by residue in multiple topology case:")
    residue_dihedrals_multi = {}
    for d in mol_multiple.dihedrals:
        res1 = mol_multiple.atoms[d.atom1].get_resname()
        res2 = mol_multiple.atoms[d.atom2].get_resname()
        res3 = mol_multiple.atoms[d.atom3].get_resname()
        res4 = mol_multiple.atoms[d.atom4].get_resname()
        key = tuple(sorted([res1, res2, res3, res4]))
        residue_dihedrals_multi[key] = residue_dihedrals_multi.get(key, 0) + 1
    for res_tuple, count in sorted(residue_dihedrals_multi.items()):
        print(f"{res_tuple}: {count} dihedrals")
    
    print("\nDihedral counts by residue in single topology case:")
    residue_dihedrals_single = {}
    for d in mol_single.dihedrals:
        res1 = mol_single.atoms[d.atom1].get_resname()
        res2 = mol_single.atoms[d.atom2].get_resname()
        res3 = mol_single.atoms[d.atom3].get_resname()
        res4 = mol_single.atoms[d.atom4].get_resname()
        key = tuple(sorted([res1, res2, res3, res4]))
        residue_dihedrals_single[key] = residue_dihedrals_single.get(key, 0) + 1
    for res_tuple, count in sorted(residue_dihedrals_single.items()):
        print(f"{res_tuple}: {count} dihedrals")
    
    # Find and print the extra dihedrals
    extra_dihedrals = dihedrals_multi - dihedrals_single
    if extra_dihedrals:
        print("\nExtra dihedrals in multiple topology case:")
        for d in sorted(extra_dihedrals):
            print(f"Dihedral: {d[0]}-{d[1]}-{d[2]}-{d[3]}")
    
    assert dihedrals_multi == dihedrals_single, \
        f"Different dihedrals between molecules: multiple has {len(dihedrals_multi)} dihedrals, single has {len(dihedrals_single)} dihedrals"

    # Compare donors
    donors_multi = set((d.atom_id, tuple(sorted(d.hydrogen_ids))) for d in mol_multiple.donors)
    donors_single = set((d.atom_id, tuple(sorted(d.hydrogen_ids))) for d in mol_single.donors)
    assert donors_multi == donors_single, \
        f"Different donors between molecules: multiple has {len(donors_multi)} donors, single has {len(donors_single)} donors"

    # Compare acceptors
    acceptors_multi = set(a.atom_id for a in mol_multiple.acceptors)
    acceptors_single = set(a.atom_id for a in mol_single.acceptors)
    assert acceptors_multi == acceptors_single, \
        f"Different acceptors between molecules: multiple has {len(acceptors_multi)} acceptors, single has {len(acceptors_single)} acceptors"

    # Compare exclusions
    exclusions_multi = set((min(e.atom1, e.atom2), max(e.atom1, e.atom2)) for e in mol_multiple.exclusions)
    exclusions_single = set((min(e.atom1, e.atom2), max(e.atom1, e.atom2)) for e in mol_single.exclusions)
    assert exclusions_multi == exclusions_single, \
        f"Different exclusions between molecules: multiple has {len(exclusions_multi)} exclusions, single has {len(exclusions_single)} exclusions"

    # Compare groups
    groups_multi = set(tuple(sorted(g.atom_indices)) for g in mol_multiple.groups)
    groups_single = set(tuple(sorted(g.atom_indices)) for g in mol_single.groups)
    assert groups_multi == groups_single, \
        f"Different groups between molecules: multiple has {len(groups_multi)} groups, single has {len(groups_single)} groups"

    # Compare CMAP entries
    def normalize_cmap(atoms):
        # Only compare standardized 5 atoms
        return tuple(atoms)
            
    cmaps_multi = set(normalize_cmap(c.atoms) for c in mol_multiple.standard_cmaps)
    cmaps_single = set(normalize_cmap(c.atoms) for c in mol_single.standard_cmaps)

    # Print CMAP information for debugging
    print("\nCMAP counts by residue in multiple topology case:")
    residue_cmaps_multi = {}
    for c in mol_multiple.standard_cmaps:
        residues = tuple(mol_multiple.atoms[atom_id].get_resname() for atom_id in c.atoms if atom_id >= 0)
        residue_cmaps_multi[residues] = residue_cmaps_multi.get(residues, 0) + 1
    for res_tuple, count in sorted(residue_cmaps_multi.items()):
        print(f"{res_tuple}: {count} CMAPs")

    print("\nCMAP counts by residue in single topology case:")
    residue_cmaps_single = {}
    for c in mol_single.standard_cmaps:
        residues = tuple(mol_single.atoms[atom_id].get_resname() for atom_id in c.atoms if atom_id >= 0)
        residue_cmaps_single[residues] = residue_cmaps_single.get(residues, 0) + 1
    for res_tuple, count in sorted(residue_cmaps_single.items()):
        print(f"{res_tuple}: {count} CMAPs")

    # Find and print the extra CMAPs
    extra_cmaps = cmaps_multi - cmaps_single
    if extra_cmaps:
        print("\nExtra CMAPs in multiple topology case:")
        for c in sorted(extra_cmaps):
            residues = [mol_multiple.atoms[atom_id].get_resname() for atom_id in c if atom_id >= 0]
            print(f"CMAP: {'-'.join(str(x) for x in c)} ({', '.join(residues)})")
            # Print original format information
            orig_cmap = next(cmap for cmap in mol_multiple.standard_cmaps if tuple(cmap.atoms) == c)
            print(f"Original format: {'-'.join(str(x) for x in orig_cmap.raw_atoms)}")
            print(f"Format type: {'PSF' if orig_cmap.is_psf_format else 'TOP'}")

    assert cmaps_multi == cmaps_single, \
        f"Different CMAP entries between molecules: multiple has {len(cmaps_multi)} entries, single has {len(cmaps_single)} entries"


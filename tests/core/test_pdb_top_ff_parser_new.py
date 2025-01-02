# tests/core/test_pdb_top_ff_parser_new.py
# tests/core/test_pdb_top_ff_parser_new.py

import os
import math
import pytest
import pygcmc as mc

@pytest.fixture(scope="module")
def test_data_dir():
    """
    Fixture to provide the path to the test data directory.
    """
    return os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "data"))


@pytest.fixture(scope="module")
def initialize_parsers(test_data_dir):
    """
    Fixture to initialize and parse all necessary files.
    Returns:
        atoms (list): List of atom objects with updated topology and force field information.
        ff_parsers (dict): Dictionary of force field parsers.
        merged_nbfix (dict): Merged NBFIX parameters.
    """
    # Initialize parsers
    pdb_parser = mc.PDBParser()
    top_parser = mc.TopParser()
    ff_parser_prot = mc.FFParser()
    ff_parser_cgenff = mc.FFParser()
    ff_parser_water = mc.FFParser()
    ff_parser_silcs = mc.FFParser()

    # Parse PDB file
    pdb_file = os.path.join(test_data_dir, "test.pdb")
    coords, residues = pdb_parser.parse(pdb_file)
    assert residues, "PDB file contains no residues"

    # Extract atoms from residues
    atoms = []
    for residue in residues:
        atoms.extend(residue.atoms)

    # Parse topology file with includes
    top_file = os.path.join(test_data_dir, "test.top")
    parse_success = top_parser.parse_with_includes(top_file)
    assert parse_success, "Failed to parse topology file"

    # Parse force field files
    prot_prm_file = os.path.join(test_data_dir, "par_all36m_prot.prm")
    parse_success = ff_parser_prot.parse(prot_prm_file)
    assert parse_success, "Failed to parse protein force field file"

    cgenff_prm_file = os.path.join(test_data_dir, "par_all36_cgenff.prm")
    parse_success = ff_parser_cgenff.parse(cgenff_prm_file)
    assert parse_success, "Failed to parse CGenFF force field file"

    water_ions_str = os.path.join(test_data_dir, "toppar_water_ions.str")
    parse_success = ff_parser_water.parse(water_ions_str)
    assert parse_success, "Failed to parse water and ions force field file"

    silcs_str = os.path.join(test_data_dir, "silcs.str")
    parse_success = ff_parser_silcs.parse(silcs_str)
    assert parse_success, "Failed to parse SILCS force field file"

    # Update atoms with topology information
    topo_updated = top_parser.update_pdb_atoms(atoms)
    assert topo_updated == len(atoms), f"Expected to update {len(atoms)} atoms, but updated {topo_updated}"

    # Update atoms with force field parameters
    ff_parser_prot.update_pdb_atoms(atoms)
    ff_parser_cgenff.update_pdb_atoms(atoms)
    ff_parser_water.update_pdb_atoms(atoms)
    ff_parser_silcs.update_pdb_atoms(atoms)

    # Collect force field parsers for summary
    ff_parsers = {
        "protein": ff_parser_prot,
        "cgenff": ff_parser_cgenff,
        "water_ions": ff_parser_water,
        "silcs": ff_parser_silcs
    }

    # Merge NBFIX parameters
    merged_nbfix = mc.FFParser.merge_nbfix_params(list(ff_parsers.values()))

    # Debug: Print merged NBFIX keys
    print(f"Merged NBFIX entries: {merged_nbfix.keys()}")

    return atoms, ff_parsers, merged_nbfix

class TestPDBTopFFParser:
    def test_number_of_residues_and_atoms(self, initialize_parsers):
        atoms, _, _ = initialize_parsers
        # From the original output: Found 22 residues and 196 atoms
        expected_residues = 22
        expected_atoms = 196
        # Count unique residues
        unique_residues = set((atom.residue, atom.sequence) for atom in atoms)
        assert len(unique_residues) == expected_residues, f"Expected {expected_residues} residues, found {len(unique_residues)}"
        assert len(atoms) == expected_atoms, f"Expected {expected_atoms} atoms, found {len(atoms)}"

    def test_specific_atom_topology(self, initialize_parsers):
        atoms, _, _ = initialize_parsers
        # Example: Atom N in residue ALA 7
        target_atom = next((atom for atom in atoms if atom.residue == "ALA" and atom.sequence == 7 and atom.name == "N"), None)
        assert target_atom is not None, "Atom N in residue ALA 7 not found"
        assert target_atom.type == "N", f"Expected type 'N', found '{target_atom.type}'"
        assert target_atom.topo_type == "NH3", f"Expected topo_type 'NH3', found '{target_atom.topo_type}'"
        assert math.isclose(target_atom.topo_charge, -0.3, abs_tol=1e-4), f"Expected topo_charge -0.3, found {target_atom.topo_charge}"
        assert math.isclose(target_atom.topo_mass, 14.007, abs_tol=1e-3), f"Expected topo_mass 14.007, found {target_atom.topo_mass}"

    def test_specific_atom_forcefield(self, initialize_parsers):
        atoms, _, _ = initialize_parsers
        # Example: Atom CA in residue ALA 7
        target_atom = next((atom for atom in atoms if atom.residue == "ALA" and atom.sequence == 7 and atom.name == "CA"), None)
        assert target_atom is not None, "Atom CA in residue ALA 7 not found"
        assert math.isclose(target_atom.forcefield_epsilon, -0.032, abs_tol=1e-4), f"Expected epsilon -0.032, found {target_atom.forcefield_epsilon}"
        assert math.isclose(target_atom.forcefield_rmin, 4.0, abs_tol=1e-4), f"Expected rmin 4.0, found {target_atom.forcefield_rmin}"

    def test_forcefield_parameters(self, initialize_parsers):
        _, ff_parsers, _ = initialize_parsers
        # Example: Protein FF parameter for atom type 'C'
        prot_ff = ff_parsers["protein"]
        c_params = prot_ff.get_nonbonded_params().get("C", None)
        assert c_params is not None, "Force field parameter for atom type 'C' not found in Protein FF"
        assert math.isclose(c_params.epsilon, -0.11, abs_tol=1e-4), f"Expected epsilon -0.11 for 'C', found {c_params.epsilon}"
        assert math.isclose(c_params.rmin, 4.0, abs_tol=1e-4), f"Expected rmin 4.0 for 'C', found {c_params.rmin}"

    def test_nbfix_parameters(self, initialize_parsers):
        _, _, merged_nbfix = initialize_parsers
        # Example: NBFIX parameter for ('SOD', 'CLA')
        key = ('SOD', 'CLA')
        sorted_key = tuple(sorted(key))
        nbfix = merged_nbfix.get(sorted_key, None)
        assert nbfix is not None, f"NBFIX parameter for {sorted_key} not found"
        assert math.isclose(nbfix.epsilon, -0.083875, abs_tol=1e-4), f"Expected epsilon -0.083875 for {sorted_key}, found {nbfix.epsilon}"
        assert math.isclose(nbfix.rmin, 3.7310, abs_tol=1e-4), f"Expected rmin 3.7310 for {sorted_key}, found {nbfix.rmin}"

    def test_multiple_nbfix_entries(self, initialize_parsers):
        _, _, merged_nbfix = initialize_parsers
        # Example: Check multiple NBFIX entries
        nbfix_entries = {
            ('BRGR1', 'HGP3'): (-0.24, 2.97),
            ('BRGR1', 'NC2'): (-1.10, 3.66),
            ('SOD', 'CLA'): (-0.083875, 3.7310),
            ('NC2', 'OC'): (-0.154919, 3.637),
            ('CLA', 'POT'): (-0.114236, 4.0810),
            # Add more entries as per your output
        }

        for key, (expected_epsilon, expected_rmin) in nbfix_entries.items():
            sorted_key = tuple(sorted(key))
            nbfix = merged_nbfix.get(sorted_key, None)
            if nbfix is None:
                print(f"NBFIX parameter for {sorted_key} not found")
            else:
                print(f"NBFIX {sorted_key}: epsilon={nbfix.epsilon}, rmin={nbfix.rmin}")
            assert nbfix is not None, f"NBFIX parameter for {sorted_key} not found"
            assert math.isclose(nbfix.epsilon, expected_epsilon, abs_tol=1e-4), f"Expected epsilon {expected_epsilon} for {sorted_key}, found {nbfix.epsilon}"
            assert math.isclose(nbfix.rmin, expected_rmin, abs_tol=1e-4), f"Expected rmin {expected_rmin} for {sorted_key}, found {nbfix.rmin}"

    def test_nonbonded_parameters_summary(self, initialize_parsers):
        _, ff_parsers, _ = initialize_parsers
        # Example: Check count of nonbonded parameters in each force field
        expected_counts = {
            "protein": 54,
            "cgenff": 161,
            "water_ions": 17,
            "silcs": 2
        }

        for ff_name, expected_count in expected_counts.items():
            ff_parser = ff_parsers.get(ff_name, None)
            assert ff_parser is not None, f"Force field parser for {ff_name} not found"
            params = ff_parser.get_nonbonded_params()
            assert len(params) == expected_count, f"Expected {expected_count} nonbonded parameters in force field '{ff_name}', found {len(params)}"

    def test_forcefield_global_parameters(self, initialize_parsers):
        _, ff_parsers, _ = initialize_parsers
        # Example: Check global nonbonded parameters
        cgenff_parser = ff_parsers.get("cgenff", None)
        assert cgenff_parser is not None, "CGenFF parser not found"

        global_params = {
            "cutnb": 14.0,
            "ctofnb": 12.0,
            "ctonnb": 10.0,
            "eps": 1.0,
            "e14fac": 1.0,
            "wmin": 1.5
        }

        for param_name, expected_value in global_params.items():
            actual_value = getattr(cgenff_parser, f"get_{param_name}")()
            assert math.isclose(actual_value, expected_value, abs_tol=1e-4), f"Expected {param_name} = {expected_value}, found {actual_value}"

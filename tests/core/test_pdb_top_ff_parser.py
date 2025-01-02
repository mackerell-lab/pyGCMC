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

    def test_water_molecule_parameters(self, initialize_parsers):
        """Test complete water molecule parameters including topology and force field."""
        atoms, _, _ = initialize_parsers
        
        # Find a complete water molecule
        water_atoms = {}
        for atom in atoms:
            if atom.residue == "SOL" and atom.sequence == 2901:  # Using a specific water molecule
                water_atoms[atom.name] = atom
        
        # Test oxygen parameters
        assert "OW" in water_atoms, "Water oxygen atom not found"
        ow = water_atoms["OW"]
        assert ow.type == "O", f"Expected type 'O', found '{ow.type}'"
        assert ow.topo_type == "OT", f"Expected topo_type 'OT', found '{ow.topo_type}'"
        assert math.isclose(ow.topo_charge, -0.834, abs_tol=1e-4)
        assert math.isclose(ow.topo_mass, 15.999, abs_tol=1e-3)
        assert math.isclose(ow.forcefield_epsilon, -0.1521, abs_tol=1e-4)
        assert math.isclose(ow.forcefield_rmin, 3.5364, abs_tol=1e-4)

        # Test hydrogen parameters
        for hw_name in ["HW1", "HW2"]:
            assert hw_name in water_atoms, f"Water hydrogen atom {hw_name} not found"
            hw = water_atoms[hw_name]
            assert hw.type == "H", f"Expected type 'H', found '{hw.type}'"
            assert hw.topo_type == "HT", f"Expected topo_type 'HT', found '{hw.topo_type}'"
            assert math.isclose(hw.topo_charge, 0.417, abs_tol=1e-4)
            assert math.isclose(hw.topo_mass, 1.008, abs_tol=1e-3)
            assert math.isclose(hw.forcefield_epsilon, -0.046, abs_tol=1e-4)
            assert math.isclose(hw.forcefield_rmin, 0.449, abs_tol=1e-3)

    def test_protein_backbone_parameters(self, initialize_parsers):
        """Test protein backbone parameters for multiple residues."""
        atoms, _, _ = initialize_parsers
        
        # Test backbone atoms for ALA-7, VAL-8, and PRO-9
        backbone_tests = [
            # residue, atom, type, topo_type, charge, mass, epsilon, rmin
            ("ALA", 7, "N", "N", "NH3", -0.3, 14.007, -0.2, 3.7),
            ("ALA", 7, "CA", "C", "CT1", 0.21, 12.011, -0.032, 4.0),
            ("ALA", 7, "C", "C", "C", 0.51, 12.011, -0.11, 4.0),
            ("ALA", 7, "O", "O", "O", -0.51, 15.999, -0.12, 3.4),
            ("VAL", 8, "N", "N", "NH1", -0.47, 14.007, -0.2, 3.7),
            ("VAL", 8, "CA", "C", "CT1", 0.07, 12.011, -0.032, 4.0),
            ("PRO", 9, "N", "N", "N", -0.29, 14.007, -0.2, 3.7),
            ("PRO", 9, "CA", "C", "CP1", 0.02, 12.011, -0.02, 4.55),
        ]
        
        for res, seq, name, type_, topo_type, charge, mass, epsilon, rmin in backbone_tests:
            atom = next((a for a in atoms 
                       if a.residue == res and a.sequence == seq and a.name == name), None)
            assert atom is not None, f"Atom {name} in {res}-{seq} not found"
            assert atom.type == type_, f"Expected type '{type_}', found '{atom.type}'"
            assert atom.topo_type == topo_type, f"Expected topo_type '{topo_type}', found '{atom.topo_type}'"
            assert math.isclose(atom.topo_charge, charge, abs_tol=1e-4)
            assert math.isclose(atom.topo_mass, mass, abs_tol=1e-3)
            assert math.isclose(atom.forcefield_epsilon, epsilon, abs_tol=1e-4)
            assert math.isclose(atom.forcefield_rmin, rmin, abs_tol=1e-4)

    def test_ligand_parameters(self, initialize_parsers):
        """Test parameters for ligand atoms (BENX and PRPX)."""
        atoms, _, _ = initialize_parsers
        
        # Test BENX (benzene) parameters
        benx_tests = [
            # atom, type, topo_type, charge, mass, epsilon, rmin
            ("CG", "C", "CG2R61", -0.115, 12.011, -0.07, 3.9848),
            ("HG", "H", "HGR61", 0.115, 1.008, -0.03, 2.7164),
            ("CD1", "C", "CG2R61", -0.115, 12.011, -0.07, 3.9848),
            ("HD1", "H", "HGR61", 0.115, 1.008, -0.03, 2.7164),
        ]
        
        for name, type_, topo_type, charge, mass, epsilon, rmin in benx_tests:
            atom = next((a for a in atoms 
                       if a.residue == "BENX" and a.sequence == 651 and a.name == name), None)
            assert atom is not None, f"BENX atom {name} not found"
            assert atom.type == type_, f"Expected type '{type_}', found '{atom.type}'"
            assert atom.topo_type == topo_type, f"Expected topo_type '{topo_type}', found '{atom.topo_type}'"
            assert math.isclose(atom.topo_charge, charge, abs_tol=1e-4)
            assert math.isclose(atom.topo_mass, mass, abs_tol=1e-3)
            assert math.isclose(atom.forcefield_epsilon, epsilon, abs_tol=1e-4)
            assert math.isclose(atom.forcefield_rmin, rmin, abs_tol=1e-4)

        # Test PRPX (propane) parameters
        prpx_tests = [
            # atom, type, topo_type, charge, mass, epsilon, rmin
            ("H11", "H", "HGA3", 0.09, 1.008, -0.024, 2.68),
            ("C1", "C", "CG331", -0.27, 12.011, -0.078, 4.1),
            ("C2", "C", "CG321", -0.18, 12.011, -0.056, 4.02),
            ("H21", "H", "HGA2", 0.09, 1.008, -0.035, 2.68),
        ]
        
        for name, type_, topo_type, charge, mass, epsilon, rmin in prpx_tests:
            atom = next((a for a in atoms 
                       if a.residue == "PRPX" and a.sequence == 792 and a.name == name), None)
            assert atom is not None, f"PRPX atom {name} not found"
            assert atom.type == type_, f"Expected type '{type_}', found '{atom.type}'"
            assert atom.topo_type == topo_type, f"Expected topo_type '{topo_type}', found '{atom.topo_type}'"
            assert math.isclose(atom.topo_charge, charge, abs_tol=1e-4)
            assert math.isclose(atom.topo_mass, mass, abs_tol=1e-3)
            assert math.isclose(atom.forcefield_epsilon, epsilon, abs_tol=1e-4)
            assert math.isclose(atom.forcefield_rmin, rmin, abs_tol=1e-4)

    def test_charge_conservation(self, initialize_parsers):
        """Test charge conservation within residues."""
        atoms, _, _ = initialize_parsers
        
        # Group atoms by residue
        residues = {}
        for atom in atoms:
            key = (atom.residue, atom.sequence)
            if key not in residues:
                residues[key] = []
            residues[key].append(atom)
        
        # Test charge conservation for each residue
        for (res_name, res_seq), res_atoms in residues.items():
            total_charge = sum(atom.topo_charge for atom in res_atoms)
            
            # Expected total charges
            if res_name == "SOL":
                assert math.isclose(total_charge, 0.0, abs_tol=1e-4), \
                    f"Water molecule {res_seq} should have neutral total charge"
            elif res_name == "BENX":
                assert math.isclose(total_charge, 0.0, abs_tol=1e-4), \
                    f"BENX molecule {res_seq} should have neutral total charge"
            elif res_name == "PRPX":
                assert math.isclose(total_charge, 0.0, abs_tol=1e-4), \
                    f"PRPX molecule {res_seq} should have neutral total charge"

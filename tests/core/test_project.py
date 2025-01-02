# tests/core/test_project_new.py

import os
import pytest
from pygcmc import Project
import sys
from io import StringIO
import time

@pytest.fixture(scope="module")
def test_data_dir():
    """Fixture to provide the path to the test data directory."""
    current_dir = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(current_dir, "..", "data")

@pytest.fixture(scope="module")
def test_files(test_data_dir):
    """Fixture to provide paths to test files."""
    files = {
        'test.pdb': os.path.join(test_data_dir, 'test.pdb'),
        'test.top': os.path.join(test_data_dir, 'test.top'),
        'silcs.str': os.path.join(test_data_dir, 'silcs.str'),
        'par_all36m_prot.prm': os.path.join(test_data_dir, 'par_all36m_prot.prm'),
        'par_all36_cgenff.prm': os.path.join(test_data_dir, 'par_all36_cgenff.prm'),
        'toppar_water_ions.str': os.path.join(test_data_dir, 'toppar_water_ions.str'),
        'water.pdb': os.path.join(test_data_dir, 'water.pdb'),
        'benx.pdb': os.path.join(test_data_dir, 'mols', 'benx.pdb')
    }
    
    # Verify all files exist
    for name, path in files.items():
        assert os.path.exists(path), f"Required test file not found: {path}"
    return files

@pytest.fixture(scope="module")
def project():
    """Fixture to provide a Project instance."""
    return Project("test_project")

class TestProject:
    def test_project_creation(self, project):
        """Test project creation and name."""
        assert project.get_name() == "test_project"

    def test_crystal_parameters_test_pdb(self, project, test_files):
        """Test crystal parameters for test.pdb."""
        structure = project.load_structure(test_files['test.pdb'])
        box = structure.get_box()
        
        assert box is not None, "test.pdb should have crystal information"
        assert len(box) == 6, "Box should have 6 parameters (a, b, c, alpha, beta, gamma)"
        assert abs(box[0] - 127.022) < 1e-3, "Incorrect a parameter in test.pdb"
        assert abs(box[1] - 133.419) < 1e-3, "Incorrect b parameter in test.pdb"
        assert abs(box[2] - 132.854) < 1e-3, "Incorrect c parameter in test.pdb"
        assert abs(box[3] - 90.0) < 1e-3, "Incorrect alpha angle in test.pdb"
        assert abs(box[4] - 90.0) < 1e-3, "Incorrect beta angle in test.pdb"
        assert abs(box[5] - 90.0) < 1e-3, "Incorrect gamma angle in test.pdb"

    def test_crystal_parameters_water_pdb(self, project, test_files):
        """Test crystal parameters for water.pdb."""
        structure = project.load_structure(test_files['water.pdb'])
        box = structure.get_box()
        
        assert box is not None, "water.pdb should have crystal information"
        assert len(box) == 6, "Water box should have 6 parameters"
        assert abs(box[0] - 10.0) < 1e-3, "Incorrect a parameter in water.pdb"
        assert abs(box[1] - 10.0) < 1e-3, "Incorrect b parameter in water.pdb"
        assert abs(box[2] - 10.0) < 1e-3, "Incorrect c parameter in water.pdb"
        assert abs(box[3] - 90.0) < 1e-3, "Incorrect alpha angle in water.pdb"
        assert abs(box[4] - 90.0) < 1e-3, "Incorrect beta angle in water.pdb"
        assert abs(box[5] - 90.0) < 1e-3, "Incorrect gamma angle in water.pdb"

    def test_crystal_parameters_benx_pdb(self, project, test_files):
        """Test crystal parameters for benx.pdb (should not have crystal info)."""
        structure = project.load_structure(test_files['benx.pdb'])
        box = structure.get_box()
        assert box is None, "benx.pdb should not have crystal information"

    def test_structure_with_topology(self, project, test_files):
        """Test loading structure with topology."""
        structure = project.load_structure(test_files['test.pdb'], test_files['test.top'])
        assert structure is not None, "Structure should be loaded successfully"
        assert len(structure.atoms) > 0, "Structure should contain atoms"

    def test_forcefield_loading_and_application(self, project, test_files):
        """Test loading and applying force field."""
        # Load structure first
        structure = project.load_structure(test_files['test.pdb'], test_files['test.top'])
        
        # Load force field
        param_files = [
            test_files['par_all36m_prot.prm'],
            test_files['par_all36_cgenff.prm'],
            test_files['toppar_water_ions.str'],
            test_files['silcs.str']
        ]
        forcefield = project.load_forcefield(param_files)
        assert forcefield is not None, "Force field should be loaded successfully"
        
        # Apply force field to structure
        structure.apply_forcefield(forcefield)
        
        # Verify force field parameters were applied
        # This could be expanded with more specific checks based on expected values
        for atom in structure.atoms:
            assert hasattr(atom, 'forcefield_epsilon'), "Atoms should have forcefield_epsilon after applying force field"
            assert hasattr(atom, 'forcefield_rmin'), "Atoms should have forcefield_rmin after applying force field" 

    def test_specific_atom_parameters(self, project, test_files):
        """Test specific atom parameters after loading structure and applying force field."""
        # Load structure with topology
        structure = project.load_structure(test_files['test.pdb'], test_files['test.top'])
        
        # Load and apply force field
        param_files = [
            test_files['par_all36m_prot.prm'],
            test_files['par_all36_cgenff.prm'],
            test_files['toppar_water_ions.str'],
            test_files['silcs.str']
        ]
        forcefield = project.load_forcefield(param_files)
        structure.apply_forcefield(forcefield)
        
        # Test ALA-7 N atom parameters
        ala_n = next((atom for atom in structure.atoms 
                     if atom.residue == "ALA" and atom.sequence == 7 and atom.name == "N"), None)
        assert ala_n is not None, "ALA-7 N atom not found"
        assert ala_n.type == "N", "Incorrect atom type"
        assert ala_n.topo_type == "NH3", "Incorrect topology type"
        assert abs(ala_n.topo_charge - (-0.3)) < 1e-6, "Incorrect topology charge"
        assert abs(ala_n.topo_mass - 14.007) < 1e-6, "Incorrect topology mass"
        assert abs(ala_n.forcefield_epsilon - (-0.2)) < 1e-6, "Incorrect epsilon"
        assert abs(ala_n.forcefield_rmin - 3.7) < 1e-6, "Incorrect rmin"

        # Test VAL-8 CA atom parameters
        val_ca = next((atom for atom in structure.atoms 
                      if atom.residue == "VAL" and atom.sequence == 8 and atom.name == "CA"), None)
        assert val_ca is not None, "VAL-8 CA atom not found"
        assert val_ca.type == "C", "Incorrect atom type"
        assert val_ca.topo_type == "CT1", "Incorrect topology type"
        assert abs(val_ca.topo_charge - 0.07) < 1e-6, "Incorrect topology charge"
        assert abs(val_ca.topo_mass - 12.011) < 1e-6, "Incorrect topology mass"
        assert abs(val_ca.forcefield_epsilon - (-0.032)) < 1e-6, "Incorrect epsilon"
        assert abs(val_ca.forcefield_rmin - 4.0) < 1e-6, "Incorrect rmin"

        # Test water molecule parameters
        water_ow = next((atom for atom in structure.atoms 
                        if atom.residue == "SOL" and atom.sequence == 2901 and atom.name == "OW"), None)
        assert water_ow is not None, "Water OW atom not found"
        assert water_ow.type == "O", "Incorrect atom type"
        assert water_ow.topo_type == "OT", "Incorrect topology type"
        assert abs(water_ow.topo_charge - (-0.834)) < 1e-6, "Incorrect topology charge"
        assert abs(water_ow.topo_mass - 15.999) < 1e-3, "Incorrect topology mass"
        assert abs(water_ow.forcefield_epsilon - (-0.1521)) < 1e-6, "Incorrect epsilon"
        assert abs(water_ow.forcefield_rmin - 3.5364) < 1e-6, "Incorrect rmin"

    def test_charge_conservation(self, project, test_files):
        """Test charge conservation within residues."""
        # Load structure with topology and force field
        structure = project.load_structure(test_files['test.pdb'], test_files['test.top'])
        param_files = [
            test_files['par_all36m_prot.prm'],
            test_files['par_all36_cgenff.prm'],
            test_files['toppar_water_ions.str'],
            test_files['silcs.str']
        ]
        forcefield = project.load_forcefield(param_files)
        structure.apply_forcefield(forcefield)
        
        # Group atoms by residue
        residues = {}
        for atom in structure.atoms:
            key = (atom.residue, atom.sequence)
            if key not in residues:
                residues[key] = []
            residues[key].append(atom)
        
        # Test charge conservation for specific residues
        test_cases = [
            ("SOL", 2901, 0.0),  # Water should have net charge 0
            ("ALA", 7, 1.0),     # ALA should have net charge +1 (N-terminal residue)
            ("VAL", 8, 0.0),     # VAL should have net charge 0
            ("BENX", 651, 0.0),  # BENX should have net charge 0
            ("PRPX", 792, 0.0)   # PRPX should have net charge 0
        ]
        
        for res_name, res_seq, expected_charge in test_cases:
            key = (res_name, res_seq)
            assert key in residues, f"Residue {res_name} {res_seq} not found"
            total_charge = sum(atom.topo_charge for atom in residues[key])
            assert abs(total_charge - expected_charge) < 1e-6, \
                f"Non-zero total charge {total_charge} for {res_name} {res_seq}"

    def test_forcefield_global_parameters(self, project, test_files):
        """Test global force field parameters."""
        # Load force field
        param_files = [
            test_files['par_all36m_prot.prm'],
            test_files['par_all36_cgenff.prm'],
            test_files['toppar_water_ions.str'],
            test_files['silcs.str']
        ]
        forcefield = project.load_forcefield(param_files)
        
        # Get global parameters
        global_params = forcefield.get_global_parameters()
        assert len(global_params) > 0, "No global parameters found"
        
        # Test expected values
        expected_params = {
            "cutoff": 14.0,
            "switching": 12.0,
            "pairlist_distance": 16.0
        }
        
        for param_name, expected_value in expected_params.items():
            assert param_name in global_params, f"Parameter {param_name} not found"
            assert abs(global_params[param_name] - expected_value) < 1e-6, \
                f"Incorrect {param_name}: expected {expected_value}, got {global_params[param_name]}"

    def test_nbfix_parameters(self, project, test_files):
        """Test NBFIX parameters in force field."""
        # Load force field
        param_files = [
            test_files['par_all36m_prot.prm'],
            test_files['par_all36_cgenff.prm'],
            test_files['toppar_water_ions.str'],
            test_files['silcs.str']
        ]
        forcefield = project.load_forcefield(param_files)
        
        # Get NBFIX parameters
        nbfix_params = forcefield.get_nbfix_parameters()
        assert len(nbfix_params) > 0, "No NBFIX parameters found"
        
        # Test specific NBFIX parameters
        test_cases = [
            # (type1, type2, epsilon, rmin)
            ("SOD", "CLA", -0.083875, 3.7310),  # Sodium-Chloride interaction
            ("BRGR1", "HGP3", -0.24, 2.97),     # Bromine-Hydrogen interaction
            ("BRGR1", "NC2", -1.10, 3.66),      # Bromine-Nitrogen interaction
            ("NC2", "OC", -0.154919, 3.637),    # Nitrogen-Oxygen interaction
            ("CLA", "POT", -0.114236, 4.0810)   # Chloride-Potassium interaction
        ]
        
        for type1, type2, expected_epsilon, expected_rmin in test_cases:
            # Find matching parameter in the list
            param = None
            for p in nbfix_params:
                if (p['type1'] == type1 and p['type2'] == type2) or \
                   (p['type1'] == type2 and p['type2'] == type1):
                    param = p
                    break
            
            assert param is not None, f"NBFIX parameters not found for {type1}-{type2}"
            assert abs(param['epsilon'] - expected_epsilon) < 1e-4, \
                f"Incorrect NBFIX epsilon for {type1}-{type2}"
            assert abs(param['rmin'] - expected_rmin) < 1e-4, \
                f"Incorrect NBFIX rmin for {type1}-{type2}" 

    def test_print_nonbonded_params(self, project, test_files, capsys):
        """Test printing of nonbonded parameters."""
        # Load structure with topology
        structure = project.load_structure(test_files['test.pdb'], test_files['test.top'])

        # Load force field
        param_files = [
            test_files['par_all36m_prot.prm'],
            test_files['par_all36_cgenff.prm'],
            test_files['toppar_water_ions.str'],
            test_files['silcs.str']
        ]
        forcefield = project.load_forcefield(param_files)
        structure.apply_forcefield(forcefield)

        # Just print nonbonded parameters
        print("\n=== Nonbonded Parameters ===")
        forcefield.print_nonbonded_params()

    def test_print_atoms(self, project, test_files, capsys):
        """Test printing of atom information."""
        # Load structure with topology
        structure = project.load_structure(test_files['test.pdb'], test_files['test.top'])

        # Load force field
        param_files = [
            test_files['par_all36m_prot.prm'],
            test_files['par_all36_cgenff.prm'],
            test_files['toppar_water_ions.str'],
            test_files['silcs.str']
        ]
        forcefield = project.load_forcefield(param_files)
        structure.apply_forcefield(forcefield)

        # Just print atom information
        print("\n=== Atom Information ===")
        project.print_all_atoms()

    def test_print_forcefield_info(self, project, test_files, capsys):
        """Test printing of force field information."""
        # Load structure with topology
        structure = project.load_structure(test_files['test.pdb'], test_files['test.top'])

        # Load force field
        param_files = [
            test_files['par_all36m_prot.prm'],
            test_files['par_all36_cgenff.prm'],
            test_files['toppar_water_ions.str'],
            test_files['silcs.str']
        ]
        forcefield = project.load_forcefield(param_files)
        structure.apply_forcefield(forcefield)

        # Just print force field information
        print("\n=== Force Field Information ===")
        project.print_forcefield_info() 
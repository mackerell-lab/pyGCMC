import os
import pytest
import math
from pygcmc import Project

@pytest.fixture
def test_data_dir():
    """Get the path to the test data directory."""
    current_dir = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(current_dir, "..", "data")

@pytest.fixture
def project():
    """Create a test project."""
    return Project("test_project")

@pytest.fixture
def param_files(test_data_dir):
    """Get list of parameter files."""
    return [
        os.path.join(test_data_dir, "par_all36m_prot.prm"),
        os.path.join(test_data_dir, "par_all36_cgenff.prm"),
        os.path.join(test_data_dir, "toppar_water_ions.str"),
        os.path.join(test_data_dir, "silcs.str")
    ]

@pytest.fixture
def structure_files(test_data_dir):
    """Get structure files."""
    return {
        'pdb': os.path.join(test_data_dir, "test.pdb"),
        'psf': os.path.join(test_data_dir, "test_proa.psf"),
        'sol_pdb': os.path.join(test_data_dir, "mols", "sol.pdb"),
        'sol_itp': os.path.join(test_data_dir, "mols", "sol.itp"),
        'benx_pdb': os.path.join(test_data_dir, "mols", "benx.pdb"),
        'benx_itp': os.path.join(test_data_dir, "mols", "benx.itp")
    }

def verify_lj_parameters(atoms_data, nb_params_dict):
    """Helper function to verify LJ parameters for atoms."""
    missing_params = []
    for atom in atoms_data:
        # Skip atoms without topology type
        if not atom['topo_type']:
            continue
        
        # Skip virtual atoms and lone pairs
        if atom['topo_type'] == 'LP' or atom.get('is_virtual', False):
            continue
        
        if atom['topo_type'] not in nb_params_dict:
            missing_params.append(atom['topo_type'])
            continue
        
        ff_params = nb_params_dict[atom['topo_type']]
        assert 'epsilon' in ff_params, f"Missing epsilon for {atom['topo_type']}"
        assert 'rmin' in ff_params, f"Missing rmin for {atom['topo_type']}"
        assert ff_params['epsilon'] != 0, f"Zero epsilon for {atom['topo_type']}"
        assert ff_params['rmin'] != 0, f"Zero rmin for {atom['topo_type']}"
    
    if missing_params:
        pytest.fail(f"Missing LJ parameters for atom types: {set(missing_params)}")

def test_multiple_structures_forcefield(project, structure_files, param_files):
    """Test applying forcefield to multiple structures."""
    # Load forcefield first
    forcefield = project.load_forcefield(param_files)
    nb_params = forcefield.get_nonbonded_parameters()
    nb_params_dict = {param['atom_type']: param for param in nb_params}
    
    # Create and load first structure (protein)
    protein_structure = project.create_structure()
    protein_structure.read_pdb(structure_files['pdb'])
    protein_structure.read_psf(structure_files['psf'])
    protein_structure.apply_forcefield(forcefield)
    
    # Create and load second structure (water)
    water_structure = project.create_structure()
    water_structure.read_pdb(structure_files['sol_pdb'])
    water_structure.read_itp(structure_files['sol_itp'])
    water_structure.apply_forcefield(forcefield)
    
    # Create and load third structure (benzene)
    benx_structure = project.create_structure()
    benx_structure.read_pdb(structure_files['benx_pdb'])
    benx_structure.read_itp(structure_files['benx_itp'])
    benx_structure.apply_forcefield(forcefield)
    
    # Verify protein structure parameters
    protein_atoms = protein_structure.get_atoms_data()
    assert len(protein_atoms) > 0, "No atoms found in protein structure"
    verify_lj_parameters(protein_atoms, nb_params_dict)
    
    # Verify specific protein parameters
    expected_protein_params = {
        'NH3': {'epsilon': -0.2, 'rmin': 3.7},    # N-terminal nitrogen
        'HC': {'epsilon': -0.046, 'rmin': 0.449}, # Hydrogen
        'CT1': {'epsilon': -0.032, 'rmin': 4.0},  # Alpha carbon
        'O': {'epsilon': -0.12, 'rmin': 3.4},     # Carbonyl oxygen
    }
    
    for atom in protein_atoms:
        if atom['topo_type'] in expected_protein_params:
            expected = expected_protein_params[atom['topo_type']]
            ff_params = nb_params_dict[atom['topo_type']]
            assert abs(ff_params['epsilon'] - expected['epsilon']) < 1e-3, \
                f"Incorrect epsilon for {atom['topo_type']}"
            assert abs(ff_params['rmin'] - expected['rmin']) < 1e-3, \
                f"Incorrect rmin for {atom['topo_type']}"
    
    # Verify water structure parameters
    water_atoms = water_structure.get_atoms_data()
    assert len(water_atoms) > 0, "No atoms found in water structure"
    verify_lj_parameters(water_atoms, nb_params_dict)
    
    # Verify specific water parameters (TIP3P model)
    expected_water_params = {
        'OT': {'epsilon': -0.1521, 'rmin': 3.5364},  # TIP3P oxygen (rmin is doubled from parameter file)
        'HT': {'epsilon': -0.0460, 'rmin': 0.4490},  # TIP3P hydrogen (rmin is doubled from parameter file)
    }
    
    for atom in water_atoms:
        if atom['topo_type'] in expected_water_params:
            expected = expected_water_params[atom['topo_type']]
            ff_params = nb_params_dict[atom['topo_type']]
            assert abs(ff_params['epsilon'] - expected['epsilon']) < 1e-3, \
                f"Incorrect epsilon for {atom['topo_type']}"
            assert abs(ff_params['rmin'] - expected['rmin']) < 1e-3, \
                f"Incorrect rmin for {atom['topo_type']}"
    
    # Verify benzene structure parameters
    benx_atoms = benx_structure.get_atoms_data()
    assert len(benx_atoms) > 0, "No atoms found in benzene structure"
    verify_lj_parameters(benx_atoms, nb_params_dict)
    
    # Verify specific benzene parameters
    expected_benx_params = {
        'CG2R61': {'epsilon': -0.070, 'rmin': 3.9848},  # Aromatic carbon (rmin is doubled from parameter file)
        'HGR61': {'epsilon': -0.030, 'rmin': 2.7164},   # Aromatic hydrogen (rmin is doubled from parameter file)
    }
    
    for atom in benx_atoms:
        if atom['topo_type'] in expected_benx_params:
            expected = expected_benx_params[atom['topo_type']]
            ff_params = nb_params_dict[atom['topo_type']]
            assert abs(ff_params['epsilon'] - expected['epsilon']) < 1e-3, \
                f"Incorrect epsilon for {atom['topo_type']}"
            assert abs(ff_params['rmin'] - expected['rmin']) < 1e-3, \
                f"Incorrect rmin for {atom['topo_type']}"

def test_forcefield_consistency(project, structure_files, param_files):
    """Test that forcefield parameters are consistent across multiple structures."""
    # Load forcefield
    forcefield = project.load_forcefield(param_files)
    
    # Create and load structures
    structures = []
    for name, files in [
        ('protein', {'pdb': structure_files['pdb'], 'top': structure_files['psf']}),
        ('water', {'pdb': structure_files['sol_pdb'], 'top': structure_files['sol_itp']}),
        ('benzene', {'pdb': structure_files['benx_pdb'], 'top': structure_files['benx_itp']})
    ]:
        structure = project.create_structure()
        structure.read_pdb(files['pdb'])
        if files['top'].endswith('.psf'):
            structure.read_psf(files['top'])
        else:
            structure.read_itp(files['top'])
        structure.apply_forcefield(forcefield)
        structures.append((name, structure))
    
    # Get all unique atom types across structures
    atom_types = set()
    for name, structure in structures:
        atoms_data = structure.get_atoms_data()
        for atom in atoms_data:
            if atom['topo_type']:
                atom_types.add(atom['topo_type'])
    
    # Verify parameters are consistent for each atom type across structures
    nb_params = forcefield.get_nonbonded_parameters()
    nb_params_dict = {param['atom_type']: param for param in nb_params}
    
    for atom_type in atom_types:
        # Find this atom type in each structure
        params_across_structures = []
        for name, structure in structures:
            atoms_data = structure.get_atoms_data()
            matching_atoms = [atom for atom in atoms_data if atom['topo_type'] == atom_type]
            if matching_atoms:
                ff_params = nb_params_dict[atom_type]
                params_across_structures.append((name, ff_params))
        
        # If atom type appears in multiple structures, verify parameters match
        if len(params_across_structures) > 1:
            base_name, base_params = params_across_structures[0]
            for other_name, other_params in params_across_structures[1:]:
                assert abs(base_params['epsilon'] - other_params['epsilon']) < 1e-6, \
                    f"Epsilon mismatch for {atom_type} between {base_name} and {other_name}"
                assert abs(base_params['rmin'] - other_params['rmin']) < 1e-6, \
                    f"Rmin mismatch for {atom_type} between {base_name} and {other_name}" 
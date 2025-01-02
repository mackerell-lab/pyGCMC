# tests/core/test_project_data.py

import os
import pytest
from pygcmc import Project
import math

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
    """Get PDB and topology files."""
    return {
        'pdb': os.path.join(test_data_dir, "test.pdb"),
        'top': os.path.join(test_data_dir, "test.top"),
        'water': os.path.join(test_data_dir, "water.pdb"),
        'benx': os.path.join(test_data_dir, "mols", "benx.pdb")
    }

@pytest.fixture
def loaded_structure(project, structure_files):
    """Load structure into project."""
    return project.load_structure(structure_files['pdb'], structure_files['top'])

@pytest.fixture
def loaded_forcefield(project, param_files):
    """Load forcefield into project."""
    return project.load_forcefield(param_files)

@pytest.fixture
def prepared_system(loaded_structure, loaded_forcefield):
    """Prepare system with structure and forcefield."""
    loaded_structure.apply_forcefield(loaded_forcefield)
    return loaded_structure, loaded_forcefield

def test_structure_data_access(prepared_system):
    """Test structure data access methods."""
    structure, ff = prepared_system
    
    # Get atoms data
    atoms_data = structure.get_atoms_data()
    assert atoms_data, "Atoms data should not be empty"
    
    # Check first atom's data
    first_atom = atoms_data[0]
    assert isinstance(first_atom, dict)
    assert all(key in first_atom for key in [
        'residue', 'sequence', 'name', 'type', 'topo_type',
        'topo_charge', 'topo_mass', 'x', 'y', 'z'
    ])
    
    # Get nonbonded parameters for verification
    nb_params = {param['atom_type']: param for param in ff.get_nonbonded_parameters()}
    
    # Verify each atom has correct data structure and valid values
    for atom in atoms_data:
        assert isinstance(atom['sequence'], int)
        assert isinstance(atom['topo_charge'], float)
        assert isinstance(atom['topo_mass'], float)
        assert isinstance(atom['x'], float)
        assert isinstance(atom['y'], float)
        assert isinstance(atom['z'], float)
        
        # Verify atom has corresponding forcefield parameters if it has a topo_type
        if atom['topo_type'] in nb_params:
            ff_params = nb_params[atom['topo_type']]
            assert 'epsilon' in ff_params
            assert 'rmin' in ff_params
            assert isinstance(ff_params['epsilon'], float)
            assert isinstance(ff_params['rmin'], float)
    
    # Test specific atom properties for known residues
    expected_atoms = {
        ('ALA', 7, 'N'): {
            'topo_type': 'NH3',
            'topo_charge': -0.3,
            'topo_mass': 14.007
        },
        ('ALA', 7, 'CA'): {
            'topo_type': 'CT1',
            'topo_charge': 0.21,
            'topo_mass': 12.011
        },
        ('ALA', 7, 'CB'): {
            'topo_type': 'CT3',
            'topo_charge': -0.27,
            'topo_mass': 12.011
        }
    }
    
    for atom in atoms_data:
        key = (atom['residue'], atom['sequence'], atom['name'])
        if key in expected_atoms:
            expected = expected_atoms[key]
            assert atom['topo_type'] == expected['topo_type'], \
                f"Incorrect topo_type for {key}"
            assert abs(atom['topo_charge'] - expected['topo_charge']) < 1e-6, \
                f"Incorrect topo_charge for {key}"
            assert abs(atom['topo_mass'] - expected['topo_mass']) < 1e-6, \
                f"Incorrect topo_mass for {key}"
            
            # Verify forcefield parameters are properly applied
            if atom['topo_type'] in nb_params:
                ff_params = nb_params[atom['topo_type']]
                assert not math.isnan(ff_params['epsilon']), \
                    f"Epsilon is NaN for {key}"
                assert not math.isnan(ff_params['rmin']), \
                    f"Rmin is NaN for {key}"

def test_forcefield_data_access(prepared_system):
    """Test forcefield data access methods."""
    _, ff = prepared_system
    
    # Test nonbonded parameters
    nb_params = ff.get_nonbonded_parameters()
    assert nb_params, "Nonbonded parameters should not be empty"
    
    # Verify structure of nonbonded parameters
    for param in nb_params:
        assert 'atom_type' in param
        assert 'epsilon' in param
        assert 'rmin' in param
        assert isinstance(param['epsilon'], float)
        assert isinstance(param['rmin'], float)
    
    # Test NBFIX parameters
    nbfix_params = ff.get_nbfix_parameters()
    if nbfix_params:  # NBFIX parameters might be empty for some forcefields
        for param in nbfix_params:
            assert 'type1' in param
            assert 'type2' in param
            assert 'epsilon' in param
            assert 'rmin' in param
            assert isinstance(param['epsilon'], float)
            assert isinstance(param['rmin'], float)
    
    # Test global parameters
    global_params = ff.get_global_parameters()
    assert global_params, "Global parameters should not be empty"
    
    # Test specific global parameter values
    expected_params = {
        'cutoff': 14.0,
        'switching': 12.0,
        'pairlist_distance': 16.0
    }
    for param_name, expected_value in expected_params.items():
        assert param_name in global_params, f"Parameter {param_name} not found"
        assert abs(global_params[param_name] - expected_value) < 1e-6, \
            f"Incorrect {param_name}: expected {expected_value}, got {global_params[param_name]}"
    
    # Test specific nonbonded parameters for common atom types
    expected_nb_params = {
        'NH3': {'epsilon': -0.2, 'rmin': 3.7},  # N-terminal nitrogen
        'HC': {'epsilon': -0.046, 'rmin': 0.449},  # Hydrogen
        'CT1': {'epsilon': -0.032, 'rmin': 4.0},   # Alpha carbon
        'O': {'epsilon': -0.12, 'rmin': 3.4},      # Carbonyl oxygen
    }
    
    nb_params_dict = {param['atom_type']: param for param in nb_params}
    for atom_type, expected in expected_nb_params.items():
        if atom_type in nb_params_dict:
            param = nb_params_dict[atom_type]
            assert abs(param['epsilon'] - expected['epsilon']) < 1e-3, \
                f"Incorrect epsilon for {atom_type}"
            assert abs(param['rmin'] - expected['rmin']) < 1e-3, \
                f"Incorrect rmin for {atom_type}"

def test_energy_data_access(prepared_system):
    """Test energy calculation data access methods."""
    structure, _ = prepared_system
    
    # Get energy components
    energy_components = structure.get_energy_components()
    assert energy_components, "Energy components should not be empty"
    
    # Get atom energy contributions
    atom_energies = structure.get_atom_energy_contributions()
    assert atom_energies, "Atom energy contributions should not be empty"
    
    # Verify structure of atom energy data
    for energy_data in atom_energies:
        assert 'atom_index' in energy_data
        assert 'vdw_energy' in energy_data
        assert 'electrostatic_energy' in energy_data
        assert 'total_energy' in energy_data
        assert isinstance(energy_data['atom_index'], int)
        assert isinstance(energy_data['vdw_energy'], float)
        assert isinstance(energy_data['electrostatic_energy'], float)
        assert isinstance(energy_data['total_energy'], float)
        # Verify total energy is sum of components
        assert abs(energy_data['total_energy'] - 
                  (energy_data['vdw_energy'] + energy_data['electrostatic_energy'])) < 1e-6

def test_residue_data_access(prepared_system):
    """Test residue data access methods."""
    structure, _ = prepared_system
    
    # Get residues data
    residues_data = structure.get_residues_data()
    assert residues_data, "Residues data should not be empty"
    
    # Verify structure of residue data
    for residue in residues_data:
        assert 'name' in residue
        assert 'sequence_number' in residue
        assert 'chain_id' in residue
        assert 'n_atoms' in residue
        assert 'center_of_mass' in residue
        
        assert isinstance(residue['sequence_number'], int)
        assert isinstance(residue['n_atoms'], int)
        assert isinstance(residue['center_of_mass'], list)
        assert len(residue['center_of_mass']) == 3
        assert all(isinstance(coord, float) for coord in residue['center_of_mass'])
        assert residue['n_atoms'] > 0 

def test_crystal_parameters(project, structure_files):
    """Test crystal parameters for different PDB files."""
    # Test test.pdb (should have crystal info)
    test_structure = project.load_structure(structure_files['pdb'])
    test_box = test_structure.get_box()
    assert test_box is not None, "test.pdb should have crystal information"
    assert len(test_box) == 6, "Box should have 6 parameters (a, b, c, alpha, beta, gamma)"
    assert abs(test_box[0] - 127.022) < 1e-3, "Incorrect a parameter in test.pdb"
    assert abs(test_box[1] - 133.419) < 1e-3, "Incorrect b parameter in test.pdb"
    assert abs(test_box[2] - 132.854) < 1e-3, "Incorrect c parameter in test.pdb"
    assert abs(test_box[3] - 90.0) < 1e-3, "Incorrect alpha angle in test.pdb"
    assert abs(test_box[4] - 90.0) < 1e-3, "Incorrect beta angle in test.pdb"
    assert abs(test_box[5] - 90.0) < 1e-3, "Incorrect gamma angle in test.pdb"

    # Test water.pdb (should have crystal info)
    water_structure = project.load_structure(structure_files['water'])
    water_box = water_structure.get_box()
    assert water_box is not None, "water.pdb should have crystal information"
    assert len(water_box) == 6, "Water box should have 6 parameters"
    assert abs(water_box[0] - 10.0) < 1e-3, "Incorrect a parameter in water.pdb"
    assert abs(water_box[1] - 10.0) < 1e-3, "Incorrect b parameter in water.pdb"
    assert abs(water_box[2] - 10.0) < 1e-3, "Incorrect c parameter in water.pdb"
    assert abs(water_box[3] - 90.0) < 1e-3, "Incorrect alpha angle in water.pdb"
    assert abs(water_box[4] - 90.0) < 1e-3, "Incorrect beta angle in water.pdb"
    assert abs(water_box[5] - 90.0) < 1e-3, "Incorrect gamma angle in water.pdb"

    # Test benx.pdb (should not have crystal info)
    benx_structure = project.load_structure(structure_files['benx'])
    benx_box = benx_structure.get_box()
    assert benx_box is None, "benx.pdb should not have crystal information" 
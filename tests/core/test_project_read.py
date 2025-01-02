# tests/core/test_project_read.py

import os
import pytest
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
    """Get PDB and topology files."""
    return {
        'pdb': os.path.join(test_data_dir, "test.pdb"),
        'top': os.path.join(test_data_dir, "test.top"),
        'water': os.path.join(test_data_dir, "water.pdb"),
        'benx': os.path.join(test_data_dir, "mols", "benx.pdb")
    }

def test_read_pdb(project, structure_files):
    """Test reading PDB file."""
    # Create a new structure
    structure = project.create_structure()
    
    # Read PDB file
    structure.read_pdb(structure_files['pdb'])
    
    # Verify structure is loaded correctly
    atoms_data = structure.get_atoms_data()
    assert len(atoms_data) > 0
    
    # Check first atom's data
    first_atom = atoms_data[0]
    assert isinstance(first_atom, dict)
    assert all(key in first_atom for key in [
        'residue', 'sequence', 'name', 'x', 'y', 'z'
    ])
    
    # Check box parameters
    box = structure.get_box()
    assert box is not None
    assert len(box) == 6
    assert abs(box[0] - 127.022) < 1e-3, "Incorrect a parameter"
    assert abs(box[1] - 133.419) < 1e-3, "Incorrect b parameter"
    assert abs(box[2] - 132.854) < 1e-3, "Incorrect c parameter"
    assert abs(box[3] - 90.0) < 1e-3, "Incorrect alpha angle"
    assert abs(box[4] - 90.0) < 1e-3, "Incorrect beta angle"
    assert abs(box[5] - 90.0) < 1e-3, "Incorrect gamma angle"

def test_read_top(project, structure_files):
    """Test reading topology file with includes."""
    # Create a new structure
    structure = project.create_structure()
    
    # Read topology file first
    structure.read_top(structure_files['top'])
    
    # Read PDB file
    structure.read_pdb(structure_files['pdb'])
    
    # Verify topology information is applied correctly
    atoms_data = structure.get_atoms_data()
    assert any(atom['topo_charge'] != 0 for atom in atoms_data), "No charges found"
    assert any(atom['topo_type'] != "" for atom in atoms_data), "No atom types found"
    
    # Check specific atoms from included topology files
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
        }
    }
    
    for atom in atoms_data:
        key = (atom['residue'], atom['sequence'], atom['name'])
        if key in expected_atoms:
            expected = expected_atoms[key]
            assert atom['topo_type'] == expected['topo_type']
            assert abs(atom['topo_charge'] - expected['topo_charge']) < 1e-6
            assert abs(atom['topo_mass'] - expected['topo_mass']) < 1e-6

def test_read_top_without_includes(project, structure_files):
    """Test reading topology file without includes."""
    # Create a new structure
    structure = project.create_structure()
    
    # Read topology file without includes
    structure.read_top_without_includes(structure_files['top'])
    
    # Read PDB file
    structure.read_pdb(structure_files['pdb'])
    
    # Verify topology information is applied correctly
    atoms_data = structure.get_atoms_data()
    
    # Check specific atoms from main topology file (ALA)
    expected_ala_atoms = {
        ('ALA', 7, 'N'): {
            'topo_type': 'NH3',
            'topo_charge': -0.3,
            'topo_mass': 14.007
        },
        ('ALA', 7, 'CA'): {
            'topo_type': 'CT1',
            'topo_charge': 0.21,
            'topo_mass': 12.011
        }
    }
    
    for atom in atoms_data:
        key = (atom['residue'], atom['sequence'], atom['name'])
        # ALA residue should have topology info as it's in the main file
        if key in expected_ala_atoms:
            expected = expected_ala_atoms[key]
            assert atom['topo_type'] == expected['topo_type'], f"Wrong type for {key}"
            assert abs(atom['topo_charge'] - expected['topo_charge']) < 1e-6, f"Wrong charge for {key}"
            assert abs(atom['topo_mass'] - expected['topo_mass']) < 1e-6, f"Wrong mass for {key}"
        # SOL/BENX/PRPX should not have topology info as they're in included files
        elif atom['residue'] in ['SOL', 'BENX', 'PRPX']:
            assert atom['topo_type'] == "", f"{atom['residue']} atoms should not have topology type without includes"
            assert 'topo_charge' not in atom or atom['topo_charge'] != atom['topo_charge'], f"{atom['residue']} atoms should have nan charge without includes"
            assert 'topo_mass' not in atom or atom['topo_mass'] != atom['topo_mass'], f"{atom['residue']} atoms should have nan mass without includes"

def test_forcefield_application(project, structure_files, param_files):
    """Test applying forcefield to structure."""
    # Create a new structure
    structure = project.create_structure()
    
    # Read PDB and topology files
    structure.read_pdb(structure_files['pdb'])
    structure.read_top(structure_files['top'])
    
    # Load and apply forcefield
    forcefield = project.load_forcefield(param_files)
    structure.apply_forcefield(forcefield)
    
    # Verify forcefield parameters are applied correctly
    atoms_data = structure.get_atoms_data()
    nb_params = forcefield.get_nonbonded_parameters()
    nb_params_dict = {param['atom_type']: param for param in nb_params}
    
    # Check specific atoms
    expected_params = {
        'NH3': {'epsilon': -0.2, 'rmin': 3.7},  # N-terminal nitrogen
        'HC': {'epsilon': -0.046, 'rmin': 0.449},  # Hydrogen
        'CT1': {'epsilon': -0.032, 'rmin': 4.0},   # Alpha carbon
        'O': {'epsilon': -0.12, 'rmin': 3.4},      # Carbonyl oxygen
    }
    
    for atom in atoms_data:
        if atom['topo_type'] in expected_params:
            expected = expected_params[atom['topo_type']]
            ff_params = nb_params_dict[atom['topo_type']]
            assert abs(ff_params['epsilon'] - expected['epsilon']) < 1e-3
            assert abs(ff_params['rmin'] - expected['rmin']) < 1e-3

def test_error_handling(project, structure_files):
    """Test error handling in structure reading functions."""
    structure = project.create_structure()
    
    # Test reading non-existent PDB file
    with pytest.raises(Exception) as excinfo:
        structure.read_pdb("non_existent.pdb")
    assert "无法打开文件" in str(excinfo.value)
    
    # Test reading non-existent topology file
    with pytest.raises(Exception) as excinfo:
        structure.read_top("non_existent.top")
    assert "Failed to parse topology file" in str(excinfo.value)
    
    # Test reading topology without includes for non-existent file
    with pytest.raises(Exception) as excinfo:
        structure.read_top_without_includes("non_existent.top")
    assert "Failed to parse topology file" in str(excinfo.value)
    
    # Test reading mismatched topology and PDB
    structure = project.create_structure()
    structure.read_top(structure_files['top'])
    with pytest.raises(Exception) as excinfo:
        structure.read_pdb(structure_files['water'])  # water.pdb doesn't match test.top
    assert "No atoms were updated with topology information" in str(excinfo.value)

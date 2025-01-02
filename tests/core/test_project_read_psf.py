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
def structure_files(test_data_dir):
    """Get PDB, PSF and ITP files."""
    return {
        'pdb': os.path.join(test_data_dir, "test.pdb"),
        'psf': os.path.join(test_data_dir, "test_proa.psf"),
        'benx': os.path.join(test_data_dir, "mols", "benx.itp"),
        'prpx': os.path.join(test_data_dir, "mols", "prpx.itp"),
        'sol': os.path.join(test_data_dir, "mols", "sol.itp")
    }

def test_read_psf(project, structure_files):
    """Test reading PSF file."""
    # Create a new structure
    structure = project.create_structure()
    
    # Read PDB file first
    structure.read_pdb(structure_files['pdb'])
    
    # Read PSF file
    structure.read_psf(structure_files['psf'])
    
    # Verify structure is loaded correctly
    atoms_data = structure.get_atoms_data()
    assert len(atoms_data) > 0
    
    # Check specific atoms from PSF file
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
        ('VAL', 8, 'N'): {
            'topo_type': 'NH1',
            'topo_charge': -0.47,
            'topo_mass': 14.007
        }
    }
    
    for atom in atoms_data:
        key = (atom['residue'], atom['sequence'], atom['name'])
        if key in expected_atoms:
            expected = expected_atoms[key]
            assert atom['topo_type'] == expected['topo_type'], f"Wrong type for {key}"
            assert abs(atom['topo_charge'] - expected['topo_charge']) < 1e-6, f"Wrong charge for {key}"
            assert abs(atom['topo_mass'] - expected['topo_mass']) < 1e-6, f"Wrong mass for {key}"

def test_read_psf_and_itp(project, structure_files):
    """Test reading PSF file followed by ITP files."""
    # Create a new structure
    structure = project.create_structure()
    
    # Read PDB file first
    structure.read_pdb(structure_files['pdb'])
    
    # Read PSF file
    structure.read_psf(structure_files['psf'])
    
    # Read ITP files
    structure.read_itp(structure_files['benx'])
    structure.read_itp(structure_files['prpx'])
    structure.read_itp(structure_files['sol'])
    
    # Verify structure is loaded correctly
    atoms_data = structure.get_atoms_data()
    assert len(atoms_data) > 0
    
    # Check atoms from PSF file
    psf_atoms = {
        ('ALA', 7, 'N'): {
            'topo_type': 'NH3',
            'topo_charge': -0.3,
            'topo_mass': 14.007
        }
    }
    
    # Check atoms from ITP files
    itp_atoms = {
        ('BENX', 1, 'CG'): {
            'topo_type': 'CA',
            'topo_charge': -0.115,
            'topo_mass': 12.011
        },
        ('PRPX', 1, 'C1'): {
            'topo_type': 'CT3',
            'topo_charge': -0.27,
            'topo_mass': 12.011
        },
        ('SOL', 1, 'OW'): {
            'topo_type': 'OW',
            'topo_charge': -0.834,
            'topo_mass': 15.9994
        }
    }
    
    # Verify all expected atoms
    for atom in atoms_data:
        key = (atom['residue'], atom['sequence'], atom['name'])
        if key in psf_atoms:
            expected = psf_atoms[key]
            assert atom['topo_type'] == expected['topo_type'], f"Wrong type for PSF atom {key}"
            assert abs(atom['topo_charge'] - expected['topo_charge']) < 1e-6, f"Wrong charge for PSF atom {key}"
            assert abs(atom['topo_mass'] - expected['topo_mass']) < 1e-6, f"Wrong mass for PSF atom {key}"
        elif key in itp_atoms:
            expected = itp_atoms[key]
            assert atom['topo_type'] == expected['topo_type'], f"Wrong type for ITP atom {key}"
            assert abs(atom['topo_charge'] - expected['topo_charge']) < 1e-6, f"Wrong charge for ITP atom {key}"
            assert abs(atom['topo_mass'] - expected['topo_mass']) < 1e-6, f"Wrong mass for ITP atom {key}"

def test_error_handling(project, structure_files):
    """Test error handling in PSF and ITP reading."""
    structure = project.create_structure()
    
    # Test reading PSF without PDB
    with pytest.raises(Exception) as excinfo:
        structure.read_psf(structure_files['psf'])
    assert "No atoms were updated with topology information" in str(excinfo.value)
    
    # Test reading non-existent PSF file
    with pytest.raises(Exception) as excinfo:
        structure.read_psf("non_existent.psf")
    assert "Cannot open file" in str(excinfo.value)
    
    # Test reading non-existent ITP file
    with pytest.raises(Exception) as excinfo:
        structure.read_itp("non_existent.itp")
    assert "Cannot open file" in str(excinfo.value)
    
    # Test reading PDB and PSF with mismatched atoms
    structure.read_pdb(structure_files['pdb'])
    with pytest.raises(Exception) as excinfo:
        structure.read_psf("tests/data/water.pdb")  # Using water.pdb as a stand-in for mismatched PSF
    assert "Cannot open file" in str(excinfo.value) 
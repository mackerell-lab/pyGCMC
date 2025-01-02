# tests/core/test_project_multi_structure_psf.py

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
    """Get PDB and PSF files."""
    return {
        'pdb': os.path.join(test_data_dir, "test.pdb"),
        'psf': os.path.join(test_data_dir, "test_proa.psf"),
        'top': os.path.join(test_data_dir, "test.top"),
        'sol_pdb': os.path.join(test_data_dir, "mols", "sol.pdb"),
        'sol_itp': os.path.join(test_data_dir, "mols", "sol.itp"),
        'benx_itp': os.path.join(test_data_dir, "mols", "benx.itp"),
        'prpx_itp': os.path.join(test_data_dir, "mols", "prpx.itp"),
        'benx': os.path.join(test_data_dir, "mols", "benx.pdb")
    }

def verify_test_structure(structure):
    """Helper function to verify test.pdb structure data."""
    # Verify structure is loaded correctly
    atoms_data = structure.get_atoms_data()
    assert len(atoms_data) > 0
    
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
    
    # Check specific atoms from PSF
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

def verify_sol_structure(structure):
    """Helper function to verify sol.pdb structure data."""
    atoms_data = structure.get_atoms_data()
    assert len(atoms_data) > 0
    
    # Check water molecules
    water_atoms = [atom for atom in atoms_data if atom['residue'] == 'SOL']
    assert len(water_atoms) > 0
    
    # Check topology information for water (TIP3P model)
    for atom in water_atoms:
        if atom['name'] == 'OW':
            assert atom['topo_type'] == 'OT'  # TIP3P oxygen type
            assert abs(atom['topo_charge'] - (-0.834)) < 1e-6
            assert abs(atom['topo_mass'] - 15.9994) < 1e-6
        elif atom['name'] == 'HW1' or atom['name'] == 'HW2':
            assert atom['topo_type'] == 'HT'  # TIP3P hydrogen type
            assert abs(atom['topo_charge'] - 0.417) < 1e-6
            assert abs(atom['topo_mass'] - 1.008) < 1e-6

def test_multiple_structures(project, structure_files, param_files):
    """Test handling of multiple structures in a project using PSF files."""
    # Create and load first structure (test.pdb)
    test_structure = project.create_structure()
    test_structure.read_pdb(structure_files['pdb'])
    test_structure.read_psf(structure_files['psf'])
    forcefield = project.load_forcefield(param_files)
    test_structure.apply_forcefield(forcefield)
    
    # Verify test structure data before loading water
    verify_test_structure(test_structure)
    
    # Create and load second structure (sol.pdb)
    sol_structure = project.create_structure()
    sol_structure.read_pdb(structure_files['sol_pdb'])
    sol_structure.read_itp(structure_files['sol_itp'])
    sol_forcefield = project.load_forcefield(param_files)
    sol_structure.apply_forcefield(sol_forcefield)
    
    # Verify sol structure data
    verify_sol_structure(sol_structure)
    
    # Verify test structure data again to ensure it wasn't affected
    verify_test_structure(test_structure)

def test_error_handling_multi_structure(project, structure_files):
    """Test error handling with multiple structures."""
    # Create first structure
    test_structure = project.create_structure()
    
    # Create second structure
    sol_structure = project.create_structure()
    
    # Test reading non-existent PDB file for first structure
    with pytest.raises(Exception) as excinfo:
        test_structure.read_pdb("non_existent.pdb")
    assert "无法打开文件" in str(excinfo.value)
    
    # Test reading non-existent PSF file for second structure
    with pytest.raises(Exception) as excinfo:
        sol_structure.read_psf("non_existent.psf")
    assert "Cannot open file" in str(excinfo.value)
    
    # Test reading mismatched PSF and PDB
    test_structure.read_psf(structure_files['psf'])
    with pytest.raises(Exception) as excinfo:
        test_structure.read_pdb(structure_files['sol_pdb'])
    assert "No atoms were updated with PSF information" in str(excinfo.value)

def test_structure_independence(project, structure_files, param_files):
    """Test that structures maintain independence in their data."""
    # Create and load first structure (test.pdb)
    test_structure = project.create_structure()
    test_structure.read_pdb(structure_files['pdb'])
    test_structure.read_psf(structure_files['psf'])
    test_forcefield = project.load_forcefield(param_files)
    test_structure.apply_forcefield(test_forcefield)
    
    # Get initial data from test structure
    test_atoms_initial = test_structure.get_atoms_data()
    test_box_initial = test_structure.get_box()
    
    # Create and load second structure (sol.pdb)
    sol_structure = project.create_structure()
    sol_structure.read_pdb(structure_files['sol_pdb'])
    sol_structure.read_itp(structure_files['sol_itp'])
    sol_forcefield = project.load_forcefield(param_files)
    sol_structure.apply_forcefield(sol_forcefield)
    
    # Get data from test structure after loading water
    test_atoms_after = test_structure.get_atoms_data()
    test_box_after = test_structure.get_box()
    
    # Verify that test structure data hasn't changed
    assert len(test_atoms_initial) == len(test_atoms_after)
    for i in range(len(test_atoms_initial)):
        for key in test_atoms_initial[i]:
            if isinstance(test_atoms_initial[i][key], (int, float, str)):
                val1 = test_atoms_initial[i][key]
                val2 = test_atoms_after[i][key]
                if isinstance(val1, float) and math.isnan(val1):
                    assert math.isnan(val2)
                else:
                    assert val1 == val2
    
    assert len(test_box_initial) == len(test_box_after)
    for i in range(len(test_box_initial)):
        assert abs(test_box_initial[i] - test_box_after[i]) < 1e-6

def test_structure_data_access(project, structure_files):
    """Test data access methods for multiple structures."""
    # Create and load structures
    test_structure = project.create_structure()
    test_structure.read_pdb(structure_files['pdb'])
    test_structure.read_psf(structure_files['psf'])
    
    sol_structure = project.create_structure()
    sol_structure.read_pdb(structure_files['sol_pdb'])
    sol_structure.read_itp(structure_files['sol_itp'])
    
    # Test coordinates access
    test_coords = test_structure.get_coordinates()
    sol_coords = sol_structure.get_coordinates()
    
    assert len(test_coords) > 0
    assert len(sol_coords) > 0
    assert len(test_coords) != len(sol_coords)  # Different number of atoms
    
    # Test box vectors access
    test_box = test_structure.get_box_vectors()
    sol_box = sol_structure.get_box_vectors()
    
    assert len(test_box) == 3
    assert len(sol_box) == 3
    
    # Verify each structure has its own box dimensions
    assert any(abs(test_box[i][j] - sol_box[i][j]) > 1e-6 
              for i in range(3) for j in range(3))

def test_topology_loading_order(project, structure_files):
    """Test different orders of loading topology and PDB files."""
    # Test structure: Load PSF first
    test_structure1 = project.create_structure()
    test_structure1.read_psf(structure_files['psf'])
    test_structure1.read_pdb(structure_files['pdb'])
    verify_test_structure(test_structure1)
    
    # Sol structure: Load PDB first
    sol_structure1 = project.create_structure()
    sol_structure1.read_pdb(structure_files['sol_pdb'])
    sol_structure1.read_itp(structure_files['sol_itp'])
    verify_sol_structure(sol_structure1)
    
    # Verify test structure wasn't affected
    verify_test_structure(test_structure1)
    
    # Create new structures with opposite loading order
    test_structure2 = project.create_structure()
    test_structure2.read_pdb(structure_files['pdb'])
    test_structure2.read_psf(structure_files['psf'])
    verify_test_structure(test_structure2)
    
    sol_structure2 = project.create_structure()
    sol_structure2.read_itp(structure_files['sol_itp'])
    sol_structure2.read_pdb(structure_files['sol_pdb'])
    verify_sol_structure(sol_structure2) 

# def test_topology_file_equivalence(project, structure_files, param_files):
#     """Test that loading the same structure with different topology files gives identical results."""
#     # Create and load the first structure using TOP file
#     top_structure = project.create_structure()
#     top_structure.read_pdb(structure_files['pdb'])
#     top_structure.read_top(structure_files['top'])

#     # Create and load the second structure using PSF file
#     psf_structure = project.create_structure()
#     psf_structure.read_pdb(structure_files['pdb'])
#     # Read PSF file first to establish base topology
#     psf_structure.read_psf(structure_files['psf'])
#     # Then read ITP files to add any additional information
#     psf_structure.read_itp(structure_files['sol_itp'])
#     psf_structure.read_itp(structure_files['benx_itp'])
#     psf_structure.read_itp(structure_files['prpx_itp'])

#     # Get atoms data from both structures
#     top_atoms = top_structure.get_atoms_data()
#     psf_atoms = psf_structure.get_atoms_data()

#     # Compare number of atoms
#     assert len(top_atoms) == len(psf_atoms), "Number of atoms mismatch"

#     # Compare atom properties
#     for i, (top_atom, psf_atom) in enumerate(zip(top_atoms, psf_atoms)):
#         # Print debug information for the first few atoms
#         if i < 5:
#             print(f"Atom {i}:")
#             print(f"  TOP: {top_atom}")
#             print(f"  PSF: {psf_atom}")

#         # Basic properties
#         assert top_atom['residue'] == psf_atom['residue'], f"Residue mismatch at atom {i}"
#         assert top_atom['sequence'] == psf_atom['sequence'], f"Sequence mismatch at atom {i}"
#         assert top_atom['name'] == psf_atom['name'], f"Name mismatch at atom {i}"

#         # Topology information
#         assert top_atom['topo_type'] == psf_atom['topo_type'], f"Topology type mismatch at atom {i}"
#         assert abs(top_atom['topo_charge'] - psf_atom['topo_charge']) < 1e-6, f"Charge mismatch at atom {i}"
#         assert abs(top_atom['topo_mass'] - psf_atom['topo_mass']) < 1e-6, f"Mass mismatch at atom {i}"

#         # Coordinates (handle NaN values)
#         for coord in ['x', 'y', 'z']:
#             top_val = top_atom[coord]
#             psf_val = psf_atom[coord]
#             if math.isnan(top_val) and math.isnan(psf_val):
#                 continue
#             assert abs(top_val - psf_val) < 1e-6, f"{coord} coordinate mismatch at atom {i}"

#     # Compare box parameters
#     top_box = top_structure.get_box()
#     psf_box = psf_structure.get_box()
#     if top_box is None:
#         assert psf_box is None, "Box parameter mismatch"
#     else:
#         assert psf_box is not None, "Box parameter mismatch"
#         for i in range(len(top_box)):
#             assert abs(top_box[i] - psf_box[i]) < 1e-6, f"Box parameter {i} mismatch"

#     # Compare box vectors
#     top_box_vectors = top_structure.get_box_vectors()
#     psf_box_vectors = psf_structure.get_box_vectors()
#     assert len(top_box_vectors) == len(psf_box_vectors), "Box vectors length mismatch"
#     for i in range(len(top_box_vectors)):
#         for j in range(3):
#             assert abs(top_box_vectors[i][j] - psf_box_vectors[i][j]) < 1e-6, f"Box vector {i},{j} mismatch" 
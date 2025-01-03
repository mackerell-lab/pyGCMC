# # tests/core/test_system_structure.py

# import os
# import pytest
# from pygcmc import System, Project

# @pytest.fixture
# def test_data_dir():
#     """Get the path to the test data directory."""
#     current_dir = os.path.dirname(os.path.abspath(__file__))
#     return os.path.join(current_dir, "..", "data")

# @pytest.fixture
# def structure_files(test_data_dir):
#     """Get structure files."""
#     return {
#         'pdb': os.path.join(test_data_dir, "test.pdb"),
#         'psf': os.path.join(test_data_dir, "test_proa.psf"),
#         'top': os.path.join(test_data_dir, "test.top"),
#         'sol_pdb': os.path.join(test_data_dir, "mols", "sol.pdb"),
#         'sol_itp': os.path.join(test_data_dir, "mols", "sol.itp"),
#     }

# @pytest.fixture
# def structure(structure_files):
#     """Create a structure with PSF topology."""
#     project = Project("test_project")
#     structure = project.create_structure()
#     structure.read_pdb(structure_files['pdb'])
#     structure.read_psf(structure_files['psf'])
#     return structure

# @pytest.fixture
# def structure_top(structure_files):
#     """Create a structure with TOP topology."""
#     project = Project("test_project")
#     structure = project.create_structure()
#     structure.read_pdb(structure_files['pdb'])
#     structure.read_top(structure_files['top'])
#     return structure

# def test_system_creation_from_structure(structure):
#     """Test system creation from a Structure object with PSF topology."""
#     system = System()
#     system.load_structure(structure)
#     verify_system_content(system)

# def test_system_creation_from_structure_top(structure_top):
#     """Test system creation from a Structure object with TOP topology."""
#     system = System()
#     system.load_structure(structure_top)
#     verify_system_content(system)

# def verify_system_content(system):
#     """Helper function to verify system content."""
#     # Verify system has content
#     assert system.get_residue_count() > 0
    
#     # Check first residue has particles
#     first_residue = system.get_residue(0)
#     assert len(first_residue.particles) > 0
    
#     # Verify particle properties
#     particle = first_residue.particles[0]
#     assert particle.mass > 0
#     assert not particle.is_virtual
#     assert all(isinstance(x, float) for x in particle.position)
#     assert all(isinstance(x, float) for x in particle.velocity)
#     assert isinstance(particle.charge, float) 
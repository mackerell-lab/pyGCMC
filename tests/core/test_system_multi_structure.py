# tests/core/test_system_multi_structure.py

import os
import pytest
from pygcmc import System, Project

@pytest.fixture
def test_data_dir():
    """Get the path to the test data directory."""
    current_dir = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(current_dir, "..", "data")

@pytest.fixture
def structure_files(test_data_dir):
    """Get structure files."""
    return {
        'pdb': os.path.join(test_data_dir, "test.pdb"),
        'psf': os.path.join(test_data_dir, "test_proa.psf"),
        'benx_psf': os.path.join(test_data_dir, "mols", "benx.psf"),
        'prpx_psf': os.path.join(test_data_dir, "mols", "prpx.psf"),
        'sol_psf': os.path.join(test_data_dir, "mols", "sol.psf"),
        'top': os.path.join(test_data_dir, "test.top"),
    }

def verify_system_content(system):
    """Helper function to verify system content."""
    # Verify system has content
    assert system.get_residue_count() > 0
    
    # Check first residue has particles
    first_residue = system.get_residue(0)
    assert len(first_residue.particles) > 0
    
    # Verify particle properties
    particle = first_residue.particles[0]
    assert particle.mass > 0
    assert not particle.is_virtual
    assert all(isinstance(x, float) for x in particle.position)
    assert all(isinstance(x, float) for x in particle.velocity)
    assert isinstance(particle.charge, float)

def test_load_structure_psf_multi(structure_files):
    """Test loading structure using multi-residue PSF approach."""
    system = System()
    system.load_structure_psf_multi(structure_files['pdb'], structure_files['psf'])
    verify_system_content(system)

def test_load_structure_psf_single(structure_files):
    """Test loading structure using single-residue PSF approach."""
    system = System()
    system.load_structure_psf_single(structure_files['pdb'], structure_files['benx_psf'], "BENX")
    verify_system_content(system)

def test_load_structure_psf_auto_multi(structure_files):
    """Test auto PSF loading with multi-residue PSF file."""
    system = System()
    system.load_structure_psf_auto(structure_files['pdb'], structure_files['psf'])
    verify_system_content(system)

def test_load_structure_psf_auto_single(structure_files):
    """Test auto PSF loading with single-residue PSF file."""
    system = System()
    system.load_structure_psf_auto(structure_files['pdb'], structure_files['benx_psf'])
    verify_system_content(system)

def test_load_structure_psf_auto_fail(structure_files):
    """Test auto PSF loading with invalid PSF file."""
    system = System()
    with pytest.raises(RuntimeError, match="No atoms found in PSF file"):
        # Use TOP file as invalid PSF file to trigger error
        system.load_structure_psf_auto(structure_files['pdb'], structure_files['top'])

def test_load_structure_psf_multi_fail(structure_files):
    """Test multi-residue PSF loading with invalid PSF file."""
    system = System()
    with pytest.raises(RuntimeError, match="No atoms found in PSF file"):
        # Use TOP file as invalid PSF file to trigger error
        system.load_structure_psf_multi(structure_files['pdb'], structure_files['top'])

def test_load_structure_psf_single_fail(structure_files):
    """Test single-residue PSF loading with invalid PSF file."""
    system = System()
    with pytest.raises(RuntimeError, match="No atoms found in PSF file"):
        # Use TOP file as invalid PSF file to trigger error
        system.load_structure_psf_single(structure_files['pdb'], structure_files['top'], "BENX")

# def test_system_multiple_psf_auto(structure_files):
#     """Test loading multiple PSF files using auto-detection method."""
#     # Create system with multiple PSF files
#     psf_files = [
#         structure_files['psf'],
#         structure_files['benx_psf'],
#         structure_files['prpx_psf'],
#         structure_files['sol_psf']
#     ]
    
#     # Test constructor method
#     system = System(pdb=structure_files['pdb'], psf=psf_files)
#     verify_system_content(system)
    
#     # Test manual loading method
#     system2 = System()
#     for psf_file in psf_files:
#         system2.load_structure_psf_auto(structure_files['pdb'], psf_file)
#     verify_system_content(system2)
    
#     # Print debug information
#     print("\nSystem 1 (Constructor method):")
#     for i in range(system.get_residue_count()):
#         res = system.get_residue(i)
#         print(f"Residue {i}: {res.name}, {len(res.particles)} particles")
#         for j, p in enumerate(res.particles):
#             print(f"  Particle {j}: mass={p.mass}, charge={p.charge}")
            
#     print("\nSystem 2 (Manual loading):")
#     for i in range(system2.get_residue_count()):
#         res = system2.get_residue(i)
#         print(f"Residue {i}: {res.name}, {len(res.particles)} particles")
#         for j, p in enumerate(res.particles):
#             print(f"  Particle {j}: mass={p.mass}, charge={p.charge}")
    
#     # Verify both methods produce same results
#     assert system.get_residue_count() == system2.get_residue_count()
#     for i in range(system.get_residue_count()):
#         res1 = system.get_residue(i)
#         res2 = system2.get_residue(i)
#         print(f"\nComparing residue {i}:")
#         print(f"Residue 1: {res1.name}, {len(res1.particles)} particles")
#         print(f"Residue 2: {res2.name}, {len(res2.particles)} particles")
#         assert len(res1.particles) == len(res2.particles)
#         for j, (p1, p2) in enumerate(zip(res1.particles, res2.particles)):
#             print(f"  Particle {j}: mass1={p1.mass}, mass2={p2.mass}")
#             assert p1.mass == p2.mass
#             assert p1.charge == p2.charge
#             assert p1.is_virtual == p2.is_virtual
#             assert all(x1 == x2 for x1, x2 in zip(p1.position, p2.position))
#             assert all(x1 == x2 for x1, x2 in zip(p1.velocity, p2.velocity))




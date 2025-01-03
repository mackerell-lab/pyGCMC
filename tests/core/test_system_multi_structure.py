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


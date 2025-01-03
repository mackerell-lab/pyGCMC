# tests/core/test_system.py

import os
import pytest
from pygcmc import System

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
        'top': os.path.join(test_data_dir, "test.top"),
        'sol_pdb': os.path.join(test_data_dir, "mols", "sol.pdb"),
        'sol_itp': os.path.join(test_data_dir, "mols", "sol.itp"),
    }

def test_system_creation():
    """Test basic system creation."""
    system = System()
    assert system is not None
    assert system.get_residue_count() == 0

def test_system_creation_with_psf(structure_files):
    """Test system creation with PDB and PSF files."""
    # Using constructor
    system1 = System(pdb=structure_files['pdb'], top=structure_files['psf'])
    verify_system_content(system1)
    
    # Using load_structure
    system2 = System()
    system2.load_structure(pdb_file=structure_files['pdb'], top_file=structure_files['psf'])
    verify_system_content(system2)

def test_system_creation_with_top(structure_files):
    """Test system creation with PDB and TOP files."""
    # Using constructor
    system1 = System(pdb=structure_files['pdb'], top=structure_files['top'])
    verify_system_content(system1)
    
    # Using load_structure
    system2 = System()
    system2.load_structure(pdb_file=structure_files['pdb'], top_file=structure_files['top'])
    verify_system_content(system2)

def test_system_creation_with_nonexistent_files():
    """Test system creation with nonexistent files."""
    with pytest.raises(RuntimeError, match="无法打开文件"):
        System(pdb="nonexistent.pdb", top="nonexistent.psf")
    with pytest.raises(RuntimeError, match="无法打开文件"):
        System(pdb="nonexistent.pdb", top="nonexistent.top")

def test_system_creation_with_invalid_topology(structure_files):
    """Test system creation with invalid topology file."""
    with pytest.raises(RuntimeError, match="Unsupported topology file format"):
        System(pdb=structure_files['pdb'], top="test.xyz")

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


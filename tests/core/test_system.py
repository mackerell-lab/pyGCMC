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

def test_load_structure_psf(structure_files):
    """Test loading structure with PSF topology."""
    system = System()
    system.load_structure(pdb_file=structure_files['pdb'], 
                         top_file=structure_files['psf'])
    
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

def test_load_structure_itp(structure_files):
    """Test loading structure with ITP topology."""
    system = System()
    system.load_structure(pdb_file=structure_files['sol_pdb'], 
                         top_file=structure_files['sol_itp'])
    
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

def test_load_structure_invalid_topology(structure_files):
    """Test loading structure with invalid topology file."""
    system = System()
    with pytest.raises(RuntimeError, match="Unsupported topology file format"):
        system.load_structure(structure_files['pdb'], "test.xyz")

def test_load_structure_nonexistent_file():
    """Test loading structure with nonexistent files."""
    system = System()
    with pytest.raises(RuntimeError, match="无法打开文件"):
        system.load_structure("nonexistent.pdb", "nonexistent.psf")


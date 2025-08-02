"""
Pytest configuration and fixtures for Drude tests
"""

import pytest
import pygcmc
import logging

# Configure logging
logger = logging.getLogger(__name__)


@pytest.fixture(autouse=True)
def drude_cleanup():
    """Automatically cleanup Drude and PGP state before and after each test"""
    # Setup: Clear before test
    pygcmc.DrudeComplete.clear()
    try:
        pygcmc.resetPGPState()
    except AttributeError:
        # resetPGPState might not be available in all versions
        pass
    
    yield
    
    # Teardown: Clear after test
    pygcmc.DrudeComplete.clear()
    try:
        pygcmc.resetPGPState()
    except AttributeError:
        pass


@pytest.fixture
def simple_drude_system():
    """Create a simple parent-Drude system for basic tests"""
    state = pygcmc.MCState()
    state.info.box = [5.0, 5.0, 5.0]
    
    # Parent atom (neutral)
    parent = pygcmc.MCAtom()
    parent.x = parent.y = parent.z = 0.0
    parent.charge = 0.0
    parent.type = 0
    
    # Drude particle
    drude = pygcmc.MCAtom()
    drude.x = drude.y = drude.z = 0.0
    drude.charge = -1.0
    drude.type = 1
    
    state.atoms = [parent, drude]
    state.activeAtomCount = 2
    
    # Create residue
    res = pygcmc.MCResidue()
    res.atomStart = 0
    res.atomCount = 2
    res.active = True
    res.type = 0
    state.residues = [res]
    state.activeResidueCount = 1
    
    return state


@pytest.fixture
def drude_scf_params():
    """Standard SCF parameters for Drude tests"""
    params = pygcmc.DrudeSCFParams()
    params.tolerance = 0.01
    params.maxIterations = 100
    params.enableHardWall = False
    params.dampingFactor = 0.5
    return params
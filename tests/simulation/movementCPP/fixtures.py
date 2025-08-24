"""
Shared fixtures and helper functions for movement C++ tests
"""

try:
    import pygcmc
    MOVEMENT_AVAILABLE = hasattr(pygcmc, 'movement')
except ImportError:
    pygcmc = None
    MOVEMENT_AVAILABLE = False


def create_movement_module():
    """Create a MovementModule instance for testing"""
    if not MOVEMENT_AVAILABLE:
        return None
    params = pygcmc.movement.MovementParams(298.15)
    params.maxAtoms = 10000
    params.maxResidues = 1000
    params.chemicalPotential = -15.7  # kJ/mol for water
    params.useCavityBias = False  # Disable for simpler testing
    params.useConfigBias = False
    return pygcmc.movement.MovementModule(params)


def create_mock_state():
    """Create a mock MCState for testing"""
    if not MOVEMENT_AVAILABLE:
        return None
    state = pygcmc.MCState()
    # Set box dimensions (in nm) - must set as array
    state.info.box = [3.0, 3.0, 3.0]
    # Set molecule type info
    state.info.max_types = 1  # We have one molecule type
    state.forcefield.numTotalTypes = 1
    state.forcefield.numMovementTypes = 1
    # Set basic LJ parameters (1x1 matrix for 1 type)
    state.forcefield.ljEps = [0.0]  # No LJ interaction for simple test
    state.forcefield.ljSigma = [0.0]  # water
    return state


def create_active_pool():
    """Create an ActivePool instance"""
    if MOVEMENT_AVAILABLE:
        return pygcmc.movement.ActivePool(1000, 100)
    return None


def create_water_molecule():
    """Create a simple water molecule for testing"""
    if not MOVEMENT_AVAILABLE:
        return []
    
    atoms = []
    # Oxygen
    o = pygcmc.MCAtom()
    o.x, o.y, o.z = 0.0, 0.0, 0.0
    o.charge = -0.834
    o.type = 0
    atoms.append(o)
    
    # Hydrogen 1  
    h1 = pygcmc.MCAtom()
    h1.x, h1.y, h1.z = 0.0957, 0.0, 0.0  # in nm
    h1.charge = 0.417
    h1.type = 1
    atoms.append(h1)
    
    # Hydrogen 2
    h2 = pygcmc.MCAtom()
    h2.x, h2.y, h2.z = -0.0239, 0.0927, 0.0  # in nm
    h2.charge = 0.417  
    h2.type = 1
    atoms.append(h2)
    
    return atoms


def create_gcmc_system():
    """Create a complete GCMC system for testing"""
    if not MOVEMENT_AVAILABLE:
        return None, None, None
    
    params = pygcmc.movement.MovementParams(298.15)
    params.chemicalPotential = -15.7  # Water
    params.useCavityBias = False  # Disable for simpler testing
    params.useConfigBias = False
    params.maxAtoms = 10000
    params.maxResidues = 1000
    
    movement = pygcmc.movement.MovementModule(params)
    state = pygcmc.MCState()
    # Set box dimensions (in nm) - must set as array
    state.info.box = [3.0, 3.0, 3.0]
    # Set molecule type info
    state.info.max_types = 1  # We have one molecule type
    state.forcefield.numTotalTypes = 1
    state.forcefield.numMovementTypes = 1
    # Set basic LJ parameters (1x1 matrix for 1 type)
    state.forcefield.ljEps = [0.0]  # No LJ interaction for simple test
    state.forcefield.ljSigma = [0.0]
    
    return movement, state, params
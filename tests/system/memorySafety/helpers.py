"""
Helper functions for memory safety tests
"""
import pygcmc
from pygcmc import MCState, MCForceField, MCAtom, MCResidue, MCMovementResidueInfo


def create_test_system(num_atoms=2, distance=1.0):
    """Create a simple test system with specified number of atoms"""
    state = MCState()
    state.info.box = [10.0, 10.0, 10.0]
    state.info.cutoff = 2.0
    state.info.setTemperature(300.0)
    
    # Create forcefield
    ff = MCForceField()
    ff.numTotalTypes = 2
    ff.numMovementTypes = 1
    ff.ljSigma = [0.3, 0.3, 0.3, 0.3]  # 2x2 matrix flattened
    ff.ljEps = [0.5, 0.5, 0.5, 0.5]    # 2x2 matrix flattened
    state.forcefield = ff
    
    atoms = []
    residues = []
    
    # Create atoms
    for i in range(num_atoms):
        atom = MCAtom()
        atom.x = 5.0 + i * distance * 0.5
        atom.y = 5.0
        atom.z = 5.0
        atom.charge = 1.0 if i % 2 == 0 else -1.0
        atom.type = i % 2
        atoms.append(atom)
        
        # Create residue for each atom
        res = MCResidue()
        res.atomStart = i
        res.atomCount = 1
        res.active = True
        res.fixed = (i < num_atoms // 2)  # First half are fixed
        res.type = i % 2
        residues.append(res)
    
    # Set state properties
    state.atoms = atoms
    state.activeAtomCount = num_atoms
    state.residues = residues
    state.activeResidueCount = num_atoms
    
    # Create movement info for moveable residues
    movement_info = MCMovementResidueInfo()
    movement_info.startIndex = num_atoms // 2
    movement_info.activeCount = num_atoms // 2
    state.movementResidues = [movement_info]
    
    return state


def init_pme_pgp_parameters(state, grid_size=64):
    """Initialize PME and PGP parameters"""
    # PME parameters
    alpha = 3.0
    mesh_size = [grid_size, grid_size, grid_size]
    spline_order = 4
    tolerance = 1e-6
    
    # Set PME parameters
    pygcmc.setPMEParameters(alpha, mesh_size, spline_order, tolerance)
    pygcmc.initializePMEParameters(state.info.cutoff, state.info.box, 
                                  alpha, mesh_size, spline_order, tolerance)
    
    # Set PGP parameters
    potential_cutoff = state.info.cutoff
    potential_grid_size = mesh_size
    pygcmc.setPGPParameters(alpha, mesh_size, potential_cutoff, 
                          potential_grid_size, spline_order, tolerance)


def compute_all_energies(state):
    """Compute all energy types and return as dict"""
    # Precompute grid potential
    pygcmc.precomputeGridPotential(state, True)
    
    # Compute movement energies
    pgp_movement = pygcmc.computeMovementEnergyPGPCompleteCorrect(state)[0]
    pme_movement = pygcmc.computeMovementEnergyPME(state)[0]
    
    # Compute system energies
    pgp_system = pygcmc.computeSystemEnergyPGPComplete(state)[0]
    pme_system = pygcmc.computeSystemEnergyPMEComplete(state)[0]
    
    return {
        'pgp_movement': pgp_movement,
        'pme_movement': pme_movement,
        'pgp_system': pgp_system,
        'pme_system': pme_system
    }
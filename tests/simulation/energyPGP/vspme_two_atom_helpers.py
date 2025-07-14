# tests/simulation/energyPGP/vspme_two_atom_helpers.py
"""Helper functions for two-atom system tests."""

import math
import pygcmc
from . import pgp_wrapper
from pygcmc import MCState, MCAtom, MCResidue
from pygcmc import MCForceField

def create_simple_two_atom_state(distance, box_size, cutoff):
    """
    Create a simple two-atom system for testing energy calculations.
    
    Args:
        distance: Distance between the two atoms (nm)
        box_size: Size of the simulation box (nm) 
        cutoff: Cutoff distance for interactions (nm)
    
    Returns:
        MCState: Configured system with two atoms
    """
    print(f"Creating two-atom system with distance {distance:.3f} nm...")
    
    # Create system
    system = MCState()
    system.info.box = [box_size, box_size, box_size]
    system.info.cutoff = cutoff
    
    # Set up force field parameters
    ff = MCForceField()
    ff.numTotalTypes = 2
    
    # LJ parameters - must match those used in theoretical calculations
    sigma = 0.4  # nm - matches theoretical calculation
    eps = 0.02   # kJ/mol - matches theoretical calculation
    
    ff.ljSigma = [sigma, sigma, sigma, sigma]
    ff.ljEps = [eps, eps, eps, eps]
    system.forcefield = ff
    
    # Create two atoms
    atoms = []
    
    # Fixed atom (at box center)
    fixed_atom = MCAtom()
    fixed_atom.x = box_size / 2.0
    fixed_atom.y = box_size / 2.0
    fixed_atom.z = box_size / 2.0
    fixed_atom.charge = 1.0
    fixed_atom.type = 0
    atoms.append(fixed_atom)
    
    # Moving atom (at specified distance from fixed atom)
    moving_atom = MCAtom()
    moving_atom.x = fixed_atom.x + distance
    moving_atom.y = fixed_atom.y
    moving_atom.z = fixed_atom.z
    moving_atom.charge = -1.0
    moving_atom.type = 1
    atoms.append(moving_atom)
    
    print(f"Fixed atom: ({fixed_atom.x:.3f}, {fixed_atom.y:.3f}, {fixed_atom.z:.3f}), charge={fixed_atom.charge}")
    print(f"Moving atom: ({moving_atom.x:.3f}, {moving_atom.y:.3f}, {moving_atom.z:.3f}), charge={moving_atom.charge}")
    
    # Verify distance
    dx = moving_atom.x - fixed_atom.x
    dy = moving_atom.y - fixed_atom.y
    dz = moving_atom.z - fixed_atom.z
    actual_distance = math.sqrt(dx*dx + dy*dy + dz*dz)
    print(f"Actual distance: {actual_distance:.6f} nm (requested: {distance:.6f} nm)")
    
    # Create residues
    residues = []
    
    # Fixed residue
    fixed_res = MCResidue()
    fixed_res.active = True
    fixed_res.fixed = True
    fixed_res.atomStart = 0
    fixed_res.atomCount = 1
    fixed_res.energy_vdw = 0.0
    fixed_res.energy_elec = 0.0
    residues.append(fixed_res)
    
    # Moving residue
    moving_res = MCResidue()
    moving_res.active = True
    moving_res.fixed = False
    moving_res.atomStart = 1
    moving_res.atomCount = 1
    moving_res.energy_vdw = 0.0
    moving_res.energy_elec = 0.0
    residues.append(moving_res)
    
    # Set system arrays
    system.atoms = atoms
    system.residues = residues
    system.activeAtomCount = len(atoms)
    system.activeResidueCount = len(residues)
    
    # Set PME/PGP parameters - use more conservative parameters
    box = [box_size, box_size, box_size]
    mesh_size = [16, 16, 16]  # Reduce grid size
    try:
        pgp_wrapper.setPMEParameters(alpha=0.2, meshSize=mesh_size, splineOrder=4, tolerance=1e-4)
        pgp_wrapper.initializePMEParameters(cutoff, box, 0.2)
        pgp_wrapper.setPGPParameters(0.2, mesh_size, state.info.cutoff, mesh_size, 4, 1e-6)
        pgp_wrapper.precomputeGridPotential(system)  # Only compute fixed atom potential
    except Exception as e:
        print(f"Warning: Error in PME/PGP setup: {e}")
    
    print(f"Two-atom system created: {len(atoms)} atoms, {len(residues)} residues")
    print(f"Within cutoff: {actual_distance < cutoff}")
    
    return system

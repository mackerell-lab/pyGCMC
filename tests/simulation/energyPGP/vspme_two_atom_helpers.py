# tests/simulation/energyPGP/vspme_two_atom_helpers.py
"""Helper functions for two-atom system tests."""

import math
import pygcmc
from pygcmc import MCState, MCAtom, MCResidue, MCForceField


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
    
    # Simple LJ parameters (neutral atoms for cleaner testing)
    sigma = 0.35  # nm
    eps = 0.1     # kJ/mol
    
    ff.ljSigma = [sigma, sigma, sigma, sigma]
    ff.ljEps = [eps, eps, eps, eps]
    system.forcefield = ff
    
    # Create two atoms
    atoms = []
    
    # First atom (fixed) at origin + offset
    atom1 = MCAtom()
    atom1.x = box_size / 2.0 - distance / 2.0
    atom1.y = box_size / 2.0
    atom1.z = box_size / 2.0
    atom1.charge = 1.0  # +1e charge
    atom1.type = 0
    atoms.append(atom1)
    
    # Second atom (moving) at the specified distance
    atom2 = MCAtom()
    atom2.x = box_size / 2.0 + distance / 2.0
    atom2.y = box_size / 2.0
    atom2.z = box_size / 2.0
    atom2.charge = -1.0  # -1e charge (opposite to create attractive interaction)
    atom2.type = 1
    atoms.append(atom2)
    
    print(f"Atom 1: ({atom1.x:.3f}, {atom1.y:.3f}, {atom1.z:.3f}), charge={atom1.charge}")
    print(f"Atom 2: ({atom2.x:.3f}, {atom2.y:.3f}, {atom2.z:.3f}), charge={atom2.charge}")
    
    # Verify distance
    dx = atom2.x - atom1.x
    dy = atom2.y - atom1.y
    dz = atom2.z - atom1.z
    actual_distance = math.sqrt(dx*dx + dy*dy + dz*dz)
    print(f"Actual distance: {actual_distance:.6f} nm (requested: {distance:.6f} nm)")
    
    # Create residues
    residues = []
    
    # Fixed residue containing first atom
    fixed_res = MCResidue()
    fixed_res.atomStart = 0
    fixed_res.atomCount = 1
    fixed_res.active = True
    fixed_res.fixed = True
    residues.append(fixed_res)
    
    # Moving residue containing second atom
    moving_res = MCResidue()
    moving_res.atomStart = 1
    moving_res.atomCount = 1
    moving_res.active = True
    moving_res.fixed = False
    residues.append(moving_res)
    
    # Set system arrays
    system.atoms = atoms
    system.residues = residues
    system.activeAtomCount = len(atoms)
    system.activeResidueCount = len(residues)
    
    print(f"Two-atom system created: {len(atoms)} atoms, {len(residues)} residues")
    print(f"Within cutoff: {actual_distance < cutoff}")
    
    return system

# tests/simulation/energyPGP/helpers.py
"""PGP test helper functions and utilities."""

import math
import random
import pygcmc
from pygcmc import MCState, MCInfo, MCAtom, MCResidue, MCForceField, MCMovementResidueInfo


def configure_energypgp_debug_logging(enabled: bool = False) -> None:
    """Configure logging for local debugging without polluting global test state."""
    pygcmc.System.set_log_level(pygcmc.LogLevel.INFO)
    pygcmc.System.set_verbose(enabled)
    pygcmc.set_platform_verbose(enabled)
    pygcmc.set_platform_log_level(pygcmc.PlatformLogLevel.INFO)
    pygcmc.set_platform_debug_mode(enabled)


def create_nacl_crystal(box_size, n_cells):
    """
    Create a NaCl crystal model
    
    Args:
        box_size: box size (nm)
        n_cells: number of unit cells in each dimension
    """
    state = MCState()
    
    # Set box size and temperature
    state.info.box = [box_size, box_size, box_size]
    state.info.setTemperature(300.0)  # 300K
    state.info.cutoff = 1.2  # 1.2 nm cutoff
    
    # Set force field parameters
    ff = MCForceField()
    ff.numTotalTypes = 2  # Na+ and Cl-
    
    # LJ parameters (from OPLS-AA force field)
    sigma_na = 0.333  # nm
    sigma_cl = 0.442  # nm
    eps_na = 0.0115  # kJ/mol
    eps_cl = 0.4184  # kJ/mol
    
    # Set LJ parameter matrix
    ff.ljSigma = [
        sigma_na, (sigma_na + sigma_cl)/2.0,
        (sigma_na + sigma_cl)/2.0, sigma_cl
    ]
    ff.ljEps = [
        eps_na, math.sqrt(eps_na * eps_cl),
        math.sqrt(eps_na * eps_cl), eps_cl
    ]
    
    state.forcefield = ff
    
    # NaCl lattice constant (0.564 nm)
    a = 0.564  
    atoms = []
    residues = []
    
    # Create NaCl lattice
    print(f"\nCreating {n_cells}x{n_cells}x{n_cells} NaCl crystal...")
    for i in range(n_cells):
        for j in range(n_cells):
            for k in range(n_cells):
                # Na+ ion
                na = MCAtom()
                na.x = i * a
                na.y = j * a
                na.z = k * a
                na.charge = 1.0
                na.type = 0
                atoms.append(na)
                
                # Cl- ion
                cl = MCAtom()
                cl.x = i * a + a/2
                cl.y = j * a + a/2
                cl.z = k * a + a/2
                cl.charge = -1.0
                cl.type = 1
                atoms.append(cl)
                
                # Create a residue for each ion pair
                res = MCResidue()
                res.atomStart = len(atoms) - 2
                res.atomCount = 2
                res.active = True
                res.fixed = False
                residues.append(res)
                
    print(f"Creation complete, added a total of {len(atoms)} atoms and {len(residues)} residues.")
    state.atoms = atoms
    state.residues = residues
    state.activeAtomCount = len(atoms)
    state.activeResidueCount = len(residues)
    
    return state


def create_long_distance_system(box_size):
    """
    Create a specialized system for testing reciprocal space calculations
    
    Atoms are placed far apart, beyond the real space cutoff distance, 
    so that reciprocal space calculations will dominate
    
    Args:
        box_size: Box size (nm)
    """
    state = MCState()
    
    # Set box size and temperature
    state.info.box = [box_size, box_size, box_size]
    state.info.setTemperature(300.0)  # 300K
    state.info.cutoff = 1.2  # 1.2 nm cutoff
    
    # Set force field parameters
    ff = MCForceField()
    ff.numTotalTypes = 2  # Two ion types
    
    # LJ parameters
    sigma_na = 0.333  # nm
    sigma_cl = 0.442  # nm
    eps_na = 0.0115  # kJ/mol
    eps_cl = 0.4184  # kJ/mol
    
    # Set LJ parameter matrix
    ff.ljSigma = [
        sigma_na, (sigma_na + sigma_cl)/2.0,
        (sigma_na + sigma_cl)/2.0, sigma_cl
    ]
    ff.ljEps = [
        eps_na, math.sqrt(eps_na * eps_cl),
        math.sqrt(eps_na * eps_cl), eps_cl
    ]
    
    state.forcefield = ff
    
    atoms = []
    residues = []
    
    # Create fixed part: two charged ions placed at opposite corners of the box
    print(f"\nCreating system with long-distance interactions...")
    
    # Ion 1: placed at one corner of the box
    ion1 = MCAtom()
    ion1.x = 0.1
    ion1.y = 0.1
    ion1.z = 0.1
    ion1.charge = 1.0
    ion1.type = 0
    atoms.append(ion1)
    
    # Ion 2: placed at the opposite corner
    ion2 = MCAtom()
    ion2.x = box_size - 0.1
    ion2.y = box_size - 0.1
    ion2.z = box_size - 0.1
    ion2.charge = -1.0
    ion2.type = 1
    atoms.append(ion2)
    
    # Create fixed residue
    fixed_res = MCResidue()
    fixed_res.atomStart = 0
    fixed_res.atomCount = 2
    fixed_res.active = True
    fixed_res.fixed = True
    residues.append(fixed_res)
    
    # Create an ion for movement, placed in the center of the box
    ion3 = MCAtom()
    ion3.x = box_size / 2.0
    ion3.y = box_size / 2.0
    ion3.z = box_size / 2.0
    ion3.charge = 1.0
    ion3.type = 0
    atoms.append(ion3)
    
    # Create moving residue
    move_res = MCResidue()
    move_res.atomStart = 2
    move_res.atomCount = 1
    move_res.active = True
    move_res.fixed = False
    residues.append(move_res)
    
    print(f"System creation complete, total of {len(atoms)} atoms and {len(residues)} residues.")
    state.atoms = atoms
    state.residues = residues
    state.activeAtomCount = len(atoms)
    state.activeResidueCount = len(residues)
    
    return state


def calculate_pbc_distance(pos1, pos2, box_size):
    """Calculate distance between two points in periodic boundary conditions"""
    dx = abs(pos1[0] - pos2[0])
    dy = abs(pos1[1] - pos2[1])
    dz = abs(pos1[2] - pos2[2])

    # Apply periodic boundary conditions
    if dx > box_size/2:
        dx = box_size - dx
    if dy > box_size/2:
        dy = box_size - dy
    if dz > box_size/2:
        dz = box_size - dz

    return math.sqrt(dx*dx + dy*dy + dz*dz)


def is_safe_position(mobile_positions, fixed_positions, cutoff, box_size):
    """Check if all distances between mobile particles and fixed particles are > cutoff"""
    for mobile_pos in mobile_positions:
        for fixed_pos, _ in fixed_positions:
            distance = calculate_pbc_distance(mobile_pos, fixed_pos, box_size)
            if distance <= cutoff:
                return False, distance
    return True, None


def generate_safe_move(current_positions, fixed_positions, cutoff, box_size, max_step=0.3):
    """Generate a safe random movement vector, ensuring all particles are at distance > cutoff from fixed particles after movement"""
    for attempt in range(100):  # Try up to 100 times
        # Generate random displacement
        dx = (random.random() - 0.5) * 2 * max_step
        dy = (random.random() - 0.5) * 2 * max_step
        dz = (random.random() - 0.5) * 2 * max_step
        
        # Calculate new positions
        new_positions = []
        for pos in current_positions:
            new_pos = (
                (pos[0] + dx) % box_size,
                (pos[1] + dy) % box_size,
                (pos[2] + dz) % box_size
            )
            new_positions.append(new_pos)
        
        # Check if new positions are safe
        is_safe, min_distance = is_safe_position(new_positions, fixed_positions, cutoff, box_size)
        if is_safe:
            return (dx, dy, dz), new_positions
    
    # If 100 attempts fail, return a smaller movement
    print("Warning: 100 attempts failed to find a safe position, using smaller movement")
    dx = 0.05
    dy = 0.05
    dz = 0.05
    new_positions = []
    for pos in current_positions:
        new_pos = (
            (pos[0] + dx) % box_size,
            (pos[1] + dy) % box_size,
            (pos[2] + dz) % box_size
        )
        new_positions.append(new_pos)
    return (dx, dy, dz), new_positions

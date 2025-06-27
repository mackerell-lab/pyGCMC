# tests/simulation/pgp/helpers.py
"""PGP test helper functions and utilities."""

import math
import random
import pygcmc
from pygcmc import MCState, MCInfo, MCAtom, MCResidue, MCForceField, MCMovementResidueInfo
import sys

# Set log level to INFO or lower to ensure detailed log output
# System log settings
pygcmc.System.set_log_level(pygcmc.LogLevel.INFO)
pygcmc.System.set_verbose(True)

# Platform log settings (for log output in energyPGP.cpp)
pygcmc.set_platform_verbose(True)  # Enable platform log output
pygcmc.set_platform_log_level(pygcmc.PlatformLogLevel.INFO)
pygcmc.set_platform_debug_mode(True)  # Enable debug mode for testing

# Ensure output buffer is flushed immediately
sys.stdout.flush()
print("Log level settings completed")
sys.stdout.flush()


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


# PGPvsPME helper functions and constants
# 物理常数
BOLTZMANN = 0.00831446261815324  # kJ/mol/K


def create_test_system(box_size):
    """
    Creates a simple system with fixed and moving ions for testing.
    Similar to create_long_distance_system but simplified.

    Args:
        box_size: Box size (nm)
    """
    from pygcmc import MCState, MCAtom, MCResidue, MCForceField
    import math
    
    state = MCState()

    # Set box size and temperature
    state.info.box = [box_size, box_size, box_size]
    state.info.setTemperature(300.0)  # 300K
    state.info.cutoff = 1.2  # Default cutoff, ensure it's set later

    # Set force field parameters (example values)
    ff = MCForceField()
    ff.numTotalTypes = 2  # Two ion types
    sigma_na = 0.333
    sigma_cl = 0.442
    eps_na = 0.0115
    eps_cl = 0.4184
    ff.ljSigma = [sigma_na, (sigma_na + sigma_cl)/2.0, (sigma_na + sigma_cl)/2.0, sigma_cl]
    ff.ljEps = [eps_na, math.sqrt(eps_na * eps_cl), math.sqrt(eps_na * eps_cl), eps_cl]
    state.forcefield = ff

    atoms = []
    residues = []

    # Fixed part: two charged ions placed far apart
    print(f"Creating system with fixed and moving parts (box={box_size}nm)...")
    ion1 = MCAtom()
    ion1.x = 0.1
    ion1.y = 0.1
    ion1.z = 0.1
    ion1.charge = 1.0
    ion1.type = 0

    ion2 = MCAtom()
    ion2.x = box_size - 0.1
    ion2.y = box_size - 0.1
    ion2.z = box_size - 0.1
    ion2.charge = -1.0
    ion2.type = 1

    atoms.extend([ion1, ion2])

    # Create fixed residue
    fixed_res = MCResidue()
    fixed_res.atomStart = 0
    fixed_res.atomCount = 2
    fixed_res.active = True
    fixed_res.fixed = True
    residues.append(fixed_res)

    # Moving part: one ion in the center
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

    print(f"System creation complete, total {len(atoms)} atoms, {len(residues)} residues.")
    state.atoms = atoms
    state.residues = residues
    state.activeAtomCount = len(atoms)
    state.activeResidueCount = len(residues)

    return state


def create_close_interaction_system(box_size, cutoff=1.2):
    """
    Creates a test system with atoms close enough to have real-space interactions.
    This system has atoms within the cutoff to test direct space electrostatics and LJ.

    Args:
        box_size: Box size (nm)
        cutoff: Cutoff distance (nm)
    """
    from pygcmc import MCState, MCAtom, MCResidue, MCForceField
    import math
    
    state = MCState()

    # Set box size and temperature
    state.info.box = [box_size, box_size, box_size]
    state.info.setTemperature(300.0)  # 300K
    state.info.cutoff = cutoff  # Set cutoff

    # Set force field parameters (example values, enhanced for LJ testing)
    ff = MCForceField()
    ff.numTotalTypes = 2  # Two atom types
    sigma_na = 0.333  # nm
    sigma_cl = 0.442  # nm
    eps_na = 0.0115   # kJ/mol - INCREASED from original values
    eps_cl = 0.4184   # kJ/mol
    # Create LJ parameter matrices
    ff.ljSigma = [
        sigma_na, (sigma_na + sigma_cl)/2.0,
        (sigma_na + sigma_cl)/2.0, sigma_cl
    ]
    ff.ljEps = [
        eps_na * 5.0,  # Increase LJ well depth for stronger interactions
        math.sqrt(eps_na * eps_cl) * 5.0,
        math.sqrt(eps_na * eps_cl) * 5.0,
        eps_cl * 5.0
    ]
    state.forcefield = ff

    atoms = []
    residues = []

    # Fixed part: two ions
    print(f"Creating close interaction system with box={box_size}nm, cutoff={cutoff}nm...")
    
    # Fixed ion 1
    ion1 = MCAtom()
    ion1.x = 0.5
    ion1.y = 0.5
    ion1.z = 0.5
    ion1.charge = 2.0  # INCREASED charge
    ion1.type = 0
    atoms.append(ion1)
    
    # Fixed ion 2
    ion2 = MCAtom()
    ion2.x = box_size - 0.5
    ion2.y = box_size - 0.5
    ion2.z = box_size - 0.5
    ion2.charge = -2.0  # INCREASED charge
    ion2.type = 1
    atoms.append(ion2)

    # Create fixed residue
    fixed_res = MCResidue()
    fixed_res.atomStart = 0
    fixed_res.atomCount = 2
    fixed_res.active = True
    fixed_res.fixed = True
    residues.append(fixed_res)

    # Moving part: ion placed VERY close to the first fixed ion (well within cutoff)
    # Place it even closer (0.1 nm) to first ion to ensure strong interactions
    ion3 = MCAtom()
    ion3.x = 0.6  # Now just 0.1 nm from ion1 at 0.5
    ion3.y = 0.5  # Same y-coordinate
    ion3.z = 0.5  # Same z-coordinate
    ion3.charge = -2.0  # INCREASED charge
    ion3.type = 1
    atoms.append(ion3)

    # Calculate distance to check within cutoff
    dx = ion3.x - ion1.x
    dy = ion3.y - ion1.y
    dz = ion3.z - ion1.z
    dist = math.sqrt(dx*dx + dy*dy + dz*dz)
    
    print(f"Distance between moving ion and fixed ion 1: {dist:.3f} nm (within cutoff: {dist < cutoff})")

    # Create moving residue
    move_res = MCResidue()
    move_res.atomStart = 2
    move_res.atomCount = 1
    move_res.active = True
    move_res.fixed = False
    residues.append(move_res)

    print(f"Close interaction system: {len(atoms)} atoms, {len(residues)} residues.")
    state.atoms = atoms
    state.residues = residues
    state.activeAtomCount = len(atoms)
    state.activeResidueCount = len(residues)

    return state


def create_very_close_system(box_size, cutoff=1.2):
    """
    Create a simple test system with atoms positioned at more reasonable distances.
    Places ions close enough to interact but not so close as to cause extreme energies.
    
    Previous version had atoms at 0.02nm which resulted in unrealistic energy values.
    """
    from pygcmc import MCState, MCAtom, MCResidue, MCForceField
    import math
    
    print(f"Creating system with reasonable interaction distances (box={box_size}nm, cutoff={cutoff}nm)...")
    
    # Create state
    system = MCState()
    system.info.box = [box_size, box_size, box_size]
    system.info.cutoff = cutoff  # nm
    
    # Create forcefield with strong LJ interactions
    ff = MCForceField()
    ff.numTotalTypes = 2
    
    # Realistic LJ parameters
    sigma_1 = 0.333  # nm (sodium-like)
    sigma_2 = 0.442  # nm (chloride-like)
    combined_sigma = (sigma_1 + sigma_2) / 2  # For mixed interactions
    
    # Use more reasonable epsilon values to prevent extreme energies
    eps = 0.5  # kJ/mol (lower value to avoid extreme LJ energies)
    
    # Set LJ parameters
    ff.ljSigma = [sigma_1, combined_sigma, combined_sigma, sigma_2]
    ff.ljEps = [eps, eps, eps, eps]
    system.forcefield = ff
    
    # Create atoms
    atoms = []
    
    # Fixed central ion
    fixed_ion = MCAtom()
    fixed_ion.x = 1.0
    fixed_ion.y = 1.0
    fixed_ion.z = 1.0
    fixed_ion.charge = 1.0  # +1 charge
    fixed_ion.type = 0
    atoms.append(fixed_ion)
    
    # Fixed counterion at a reasonable distance (well outside extreme LJ region)
    fixed_counterion = MCAtom()
    fixed_counterion.x = 4.0  # Far enough away to not affect the test
    fixed_counterion.y = 4.0
    fixed_counterion.z = 4.0
    fixed_counterion.charge = -1.0  # -1 charge
    fixed_counterion.type = 1
    atoms.append(fixed_counterion)
    
    # Moving ion at a distance where interactions are meaningful but not extreme
    # About 0.7nm is a good balance - within cutoff but not in extreme repulsion
    moving_ion = MCAtom()
    moving_ion.x = 1.7  # Distance of ~0.7nm from fixed_ion
    moving_ion.y = 1.0
    moving_ion.z = 1.0
    moving_ion.charge = -1.0  # -1 charge
    moving_ion.type = 1
    atoms.append(moving_ion)
    
    # Calculate actual distance for verification
    dx = moving_ion.x - fixed_ion.x
    dy = moving_ion.y - fixed_ion.y
    dz = moving_ion.z - fixed_ion.z
    dist = math.sqrt(dx*dx + dy*dy + dz*dz)
    print(f"Distance between moving ion and fixed ion 1: {dist:.3f} nm (within cutoff: {dist < cutoff})")
    print(f"This distance allows meaningful interactions without extreme energy values")
    
    # Create residues
    residues = []
    
    # Fixed ions residue
    fixed_res = MCResidue()
    fixed_res.atomStart = 0
    fixed_res.atomCount = 2  # Both fixed ions in one residue
    fixed_res.active = True
    fixed_res.fixed = True
    residues.append(fixed_res)
    
    # Moving ion residue
    moving_res = MCResidue()
    moving_res.atomStart = 2
    moving_res.atomCount = 1
    moving_res.active = True
    moving_res.fixed = False
    residues.append(moving_res)
    
    # Set system state
    system.atoms = atoms
    system.residues = residues
    system.activeAtomCount = len(atoms)
    system.activeResidueCount = len(residues)
    
    print(f"Close interaction system created: {len(atoms)} atoms, {len(residues)} residues.")
    return system


def create_simple_two_atom_state(distance, box_size, cutoff):
    """Create a simple system with exactly two atoms separated by the given distance."""
    from pygcmc import MCState, MCAtom, MCResidue, MCForceField
    import math
    
    print(f"Creating two-atom system with distance {distance:.3f} nm in {box_size:.1f} nm box")
    
    # Create the system
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
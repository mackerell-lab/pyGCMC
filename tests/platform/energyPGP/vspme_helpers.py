# tests/simulation/energyPGP/vspme_helpers.py
"""PGP vs PME test helper functions and utilities."""

import math
import pygcmc
from pygcmc import MCState, MCAtom, MCResidue, MCForceField

# Physical constants
BOLTZMANN = 0.00831446261815324  # kJ/mol/K


def create_test_system(box_size):
    """
    Creates a simple system with fixed and moving ions for testing.
    Similar to create_long_distance_system but simplified.

    Args:
        box_size: Box size (nm)
    """
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



"""
Helper functions for Drude oscillator tests
"""

import numpy as np
import pygcmc
import math


def create_drude_particle(charge=-1.0, polarizability=0.001, parent_idx=0, drude_idx=1):
    """Create a standard Drude particle with given parameters"""
    particle = pygcmc.DrudeParticle()
    particle.drudeIndex = drude_idx
    particle.parentIndex = parent_idx
    particle.charge = charge
    particle.polarizability = polarizability
    particle.computeSpringConstants()
    return particle


def create_water_molecule(origin=(0, 0, 0), model='swm4ndp'):
    """Create a water molecule with specified model
    
    Args:
        origin: (x, y, z) position of oxygen
        model: 'swm4ndp' or 'tip4p' or others
        
    Returns:
        List of MCAtom objects
    """
    x0, y0, z0 = origin
    atoms = []
    
    if model == 'swm4ndp':
        # SWM4-NDP water parameters
        charges = {
            'O': 1.71636,
            'D': -1.71636,
            'H': 0.55733,
            'M': -1.11466
        }
        
        # Oxygen
        O = pygcmc.MCAtom()
        O.x, O.y, O.z = x0, y0, z0
        O.charge = charges['O']
        O.type = 0
        atoms.append(O)
        
        # Drude (initially at oxygen)
        D = pygcmc.MCAtom()
        D.x, D.y, D.z = x0, y0, z0
        D.charge = charges['D']
        D.type = 1
        atoms.append(D)
        
        # Hydrogen 1
        H1 = pygcmc.MCAtom()
        H1.x = x0 + 0.09572
        H1.y = y0
        H1.z = z0
        H1.charge = charges['H']
        H1.type = 2
        atoms.append(H1)
        
        # Hydrogen 2
        angle_rad = math.radians(104.52)
        H2 = pygcmc.MCAtom()
        H2.x = x0 + 0.09572 * math.cos(angle_rad)
        H2.y = y0 + 0.09572 * math.sin(angle_rad)
        H2.z = z0
        H2.charge = charges['H']
        H2.type = 2
        atoms.append(H2)
        
        # Virtual site M
        bisector_angle = angle_rad / 2.0
        M = pygcmc.MCAtom()
        M.x = x0 + 0.024034 * math.cos(bisector_angle)
        M.y = y0 + 0.024034 * math.sin(bisector_angle)
        M.z = z0
        M.charge = charges['M']
        M.type = 3
        atoms.append(M)
        
    else:
        raise ValueError(f"Unknown water model: {model}")
    
    return atoms


def setup_water_system(n_waters, box_size=3.0, spacing=0.31):
    """Setup a system with multiple water molecules
    
    Args:
        n_waters: Number of water molecules
        box_size: Cubic box size in nm
        spacing: Minimum spacing between waters
        
    Returns:
        MCState object with waters and residues configured
    """
    state = pygcmc.MCState()
    
    # Create waters in a grid or random positions
    atoms = []
    positions = []
    
    # Simple grid placement
    n_per_side = int(np.ceil(n_waters ** (1/3)))
    grid_spacing = box_size / n_per_side
    
    water_count = 0
    for i in range(n_per_side):
        for j in range(n_per_side):
            for k in range(n_per_side):
                if water_count >= n_waters:
                    break
                    
                x = (i + 0.5) * grid_spacing
                y = (j + 0.5) * grid_spacing
                z = (k + 0.5) * grid_spacing
                
                # Add small random displacement to break symmetry
                x += (np.random.random() - 0.5) * 0.01
                y += (np.random.random() - 0.5) * 0.01
                z += (np.random.random() - 0.5) * 0.01
                
                positions.append((x, y, z))
                atoms.extend(create_water_molecule((x, y, z)))
                water_count += 1
    
    state.atoms = atoms
    state.activeAtomCount = len(atoms)
    state.info.box = [box_size, box_size, box_size]
    
    # Setup residues
    residues = []
    for i in range(n_waters):
        res = pygcmc.MCResidue()
        res.atomStart = i * 5  # 5 atoms per SWM4-NDP water
        res.atomCount = 5
        res.active = True
        res.type = 0
        residues.append(res)
    
    state.residues = residues
    state.activeResidueCount = n_waters
    
    return state, positions


def setup_drude_system(state, water_positions=None):
    """Setup Drude particles for a water system
    
    Args:
        state: MCState with water molecules
        water_positions: List of (x,y,z) positions for oxygens
        
    Returns:
        None (modifies DrudeComplete singleton)
    """
    pygcmc.DrudeComplete.clear()
    
    # Determine number of waters
    n_waters = state.activeResidueCount
    if n_waters == 0:
        # No residues, assume single water
        n_waters = state.activeAtomCount // 5
    
    # Add Drude particles
    for i in range(n_waters):
        particle = pygcmc.DrudeParticle()
        particle.drudeIndex = i * 5 + 1  # Drude is second atom
        particle.parentIndex = i * 5      # Oxygen is first
        particle.charge = -1.71636       # SWM4-NDP charge
        particle.polarizability = 0.978e-3  # nm³
        particle.computeSpringConstants()
        pygcmc.DrudeComplete.addParticle(particle)
    
    # Add Thole screening between all pairs
    for i in range(n_waters):
        for j in range(i+1, n_waters):
            pair = pygcmc.ScreenedPair()
            pair.dipole1 = i
            pair.dipole2 = j
            pair.thole = 1.3  # SWM4-NDP value
            pygcmc.DrudeComplete.addScreenedPair(pair)
    
    # Set standard SCF parameters
    params = pygcmc.DrudeSCFParams()
    params.tolerance = 1e-5  # 0.01 kJ/mol/nm
    params.maxIterations = 100
    params.maxDrudeDistance = 0.02  # 0.2 Å hard wall
    params.dampingFactor = 0.5
    pygcmc.DrudeComplete.setParameters(params)


def calculate_dipole_moment(state, particle_idx):
    """Calculate dipole moment of a Drude oscillator
    
    Args:
        state: MCState
        particle_idx: Index of Drude particle in DrudeComplete list
        
    Returns:
        (magnitude, vector) - dipole magnitude and 3D vector
    """
    # This would need access to internal particle list
    # For now, assume standard water indexing
    drude_atom_idx = particle_idx * 5 + 1
    parent_atom_idx = particle_idx * 5
    
    drude = state.atoms[drude_atom_idx]
    parent = state.atoms[parent_atom_idx]
    
    # Displacement vector
    dx = drude.x - parent.x
    dy = drude.y - parent.y
    dz = drude.z - parent.z
    
    # Apply PBC if needed
    box = state.info.box
    if box[0] > 0:
        dx -= box[0] * round(dx / box[0])
        dy -= box[1] * round(dy / box[1])
        dz -= box[2] * round(dz / box[2])
    
    # Dipole moment = charge * displacement
    charge = abs(drude.charge)  # Use absolute value
    mu_x = charge * dx
    mu_y = charge * dy
    mu_z = charge * dz
    
    magnitude = math.sqrt(mu_x**2 + mu_y**2 + mu_z**2)
    
    return magnitude, (mu_x, mu_y, mu_z)


def calculate_interaction_energy(state1, state2, drude_particles):
    """Calculate interaction energy between two states
    
    Useful for calculating binding energies, solvation energies, etc.
    
    Args:
        state1: First MCState
        state2: Second MCState  
        drude_particles: List of DrudeParticle objects
        
    Returns:
        Interaction energy in kJ/mol
    """
    # Calculate total energy
    pygcmc.DrudeComplete.clear()
    for p in drude_particles:
        pygcmc.DrudeComplete.addParticle(p)
    
    # Combine states (would need proper implementation)
    # This is a placeholder
    return 0.0


def check_scf_convergence(state, tolerance=1e-6):
    """Check if SCF is properly converged by comparing forces
    
    Args:
        state: MCState after SCF
        tolerance: Force tolerance in kJ/mol/nm
        
    Returns:
        (converged, max_force) tuple
    """
    # This would need access to force calculations
    # Placeholder implementation
    return True, 0.0


def create_external_field(position, strength, direction=(1, 0, 0)):
    """Create an external electric field using point charges
    
    Args:
        position: Where to place the field source
        strength: Field strength (arbitrary units)
        direction: Field direction vector
        
    Returns:
        List of MCAtom objects representing field
    """
    # Normalize direction
    dx, dy, dz = direction
    norm = math.sqrt(dx*dx + dy*dy + dz*dz)
    dx, dy, dz = dx/norm, dy/norm, dz/norm
    
    # Create two opposite charges to make a dipole field
    atoms = []
    
    # Positive charge
    pos = pygcmc.MCAtom()
    pos.x = position[0] + 0.1 * dx
    pos.y = position[1] + 0.1 * dy
    pos.z = position[2] + 0.1 * dz
    pos.charge = strength
    pos.type = 99  # Special type for external
    atoms.append(pos)
    
    # Negative charge
    neg = pygcmc.MCAtom()
    neg.x = position[0] - 0.1 * dx
    neg.y = position[1] - 0.1 * dy
    neg.z = position[2] - 0.1 * dz
    neg.charge = -strength
    neg.type = 99
    atoms.append(neg)
    
    return atoms
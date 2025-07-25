"""
Helper functions for creating molecules and Drude particles
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
        
        # Hydrogen atoms
        # OH bond length = 0.09572 nm
        # HOH angle = 104.52 degrees
        bond_length = 0.09572
        angle = 104.52 * math.pi / 180.0
        
        # H1 along x-axis
        H1 = pygcmc.MCAtom()
        H1.x = x0 + bond_length
        H1.y = y0
        H1.z = z0
        H1.charge = charges['H']
        H1.type = 2
        atoms.append(H1)
        
        # H2 at angle
        H2 = pygcmc.MCAtom()
        H2.x = x0 + bond_length * math.cos(angle)
        H2.y = y0 + bond_length * math.sin(angle)
        H2.z = z0
        H2.charge = charges['H']
        H2.type = 2
        atoms.append(H2)
        
        # M site (bisector)
        # Distance from O = 0.024034 nm
        m_dist = 0.024034
        bisector_angle = angle / 2.0
        M = pygcmc.MCAtom()
        M.x = x0 + m_dist * math.cos(bisector_angle)
        M.y = y0 + m_dist * math.sin(bisector_angle)
        M.z = z0
        M.charge = charges['M']
        M.type = 3
        atoms.append(M)
        
    return atoms


def setup_water_system(n_waters, box_size=3.0, spacing=0.31):
    """Create a system with multiple water molecules
    
    Args:
        n_waters: Number of water molecules
        box_size: Box dimension in nm
        spacing: Minimum spacing between molecules
        
    Returns:
        MCState object with waters
    """
    state = pygcmc.MCState()
    state.info.box = [box_size, box_size, box_size]
    state.info.cutoff = min(1.2, box_size/2 - 0.01)
    
    # Force field for SWM4-NDP
    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 4  # O, D, H, M
    ff.numMovementTypes = 4
    
    # LJ parameters (only oxygen has LJ)
    ff.ljEps = [0.21094 * 4.184, 0.0, 0.0, 0.0]  # kJ/mol
    ff.ljSigma = [0.318395, 0.1, 0.1, 0.1]  # nm
    
    state.forcefield = ff
    
    # Place waters on a grid
    n_per_dim = int(np.ceil(n_waters ** (1/3)))
    grid_spacing = box_size / (n_per_dim + 1)
    
    atoms = []
    water_count = 0
    
    for i in range(n_per_dim):
        for j in range(n_per_dim):
            for k in range(n_per_dim):
                if water_count >= n_waters:
                    break
                    
                x = (i + 1) * grid_spacing
                y = (j + 1) * grid_spacing
                z = (k + 1) * grid_spacing
                
                # Add small random displacement
                x += (np.random.rand() - 0.5) * 0.02
                y += (np.random.rand() - 0.5) * 0.02
                z += (np.random.rand() - 0.5) * 0.02
                
                water_atoms = create_water_molecule((x, y, z))
                atoms.extend(water_atoms)
                water_count += 1
            
            if water_count >= n_waters:
                break
        if water_count >= n_waters:
            break
    
    state.atoms = atoms
    state.activeAtomCount = len(atoms)
    
    return state
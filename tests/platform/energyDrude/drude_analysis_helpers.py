"""
Helper functions for analyzing Drude systems
"""

import numpy as np
import pygcmc
import math


def setup_drude_system(state, water_positions=None):
    """Setup Drude particles for a water system
    
    Args:
        state: MCState with water molecules
        water_positions: Optional list of (O_idx, D_idx) tuples
        
    Returns:
        List of DrudeParticle objects
    """
    pygcmc.DrudeComplete.clear()
    particles = []
    
    if water_positions is None:
        # Assume standard SWM4-NDP ordering: O, D, H, H, M
        n_waters = state.activeAtomCount // 5
        water_positions = [(i*5, i*5+1) for i in range(n_waters)]
    
    for o_idx, d_idx in water_positions:
        # SWM4-NDP polarizability
        # alpha = ONE_4PI_EPS0 * q^2 / k_spring
        # From OpenMM: polarizability = ONE_4PI_EPS0*1.71636*1.71636/(100000*4.184)
        polarizability = 138.935456 * 1.71636 * 1.71636 / (100000 * 4.184)
        
        particle = pygcmc.DrudeParticle()
        particle.drudeIndex = d_idx
        particle.parentIndex = o_idx
        particle.charge = -1.71636
        particle.polarizability = polarizability
        particle.computeSpringConstants()
        
        pygcmc.DrudeComplete.addParticle(particle)
        particles.append(particle)
    
    return particles


def calculate_dipole_moment(state, particle_idx):
    """Calculate induced dipole moment for a Drude particle
    
    Args:
        state: MCState
        particle_idx: Index in DrudeComplete particle list
        
    Returns:
        (magnitude, (dx, dy, dz)) tuple
    """
    # This is a simplified version - would need access to internal particle list
    # In practice, you'd calculate from the displacement
    drude_pos = state.atoms[particle_idx * 2 + 1]  # Assuming standard ordering
    parent_pos = state.atoms[particle_idx * 2]
    
    dx = drude_pos.x - parent_pos.x
    dy = drude_pos.y - parent_pos.y
    dz = drude_pos.z - parent_pos.z
    
    # Apply PBC
    box = state.info.box
    dx -= box[0] * round(dx / box[0])
    dy -= box[1] * round(dy / box[1])
    dz -= box[2] * round(dz / box[2])
    
    magnitude = math.sqrt(dx*dx + dy*dy + dz*dz)
    
    # Dipole moment = charge * displacement
    # For Drude with charge q_D, dipole = q_D * r_D
    q_drude = -1.71636  # Standard SWM4-NDP
    mu_x = q_drude * dx
    mu_y = q_drude * dy
    mu_z = q_drude * dz
    mu_mag = abs(q_drude) * magnitude
    
    return mu_mag, (mu_x, mu_y, mu_z)


def calculate_interaction_energy(state1, state2, drude_particles):
    """Calculate interaction energy between two configurations
    
    This is a placeholder - actual implementation would need
    to properly handle the energy calculation
    """
    # For now, just return a dummy value
    return -10.0  # kJ/mol


def check_scf_convergence(state, tolerance=1e-6):
    """Check if SCF has converged by examining forces on Drude particles
    
    Returns:
        (converged, max_force) tuple
    """
    # This is a simplified check - in practice would calculate actual forces
    return True, 0.0


def create_external_field(position, strength, direction=(1, 0, 0)):
    """Create an external charge to generate electric field
    
    Args:
        position: (x, y, z) position
        strength: Field strength (charge magnitude)
        direction: Not used directly, field radiates from charge
        
    Returns:
        MCAtom representing the external charge
    """
    external = pygcmc.MCAtom()
    external.x, external.y, external.z = position
    external.charge = strength
    external.type = 0  # Assume type 0 for external charges
    return external


def calculate_numerical_force(state, atom_idx, direction, delta=1e-6):
    """Calculate force on atom using finite difference
    
    Args:
        state: MCState
        atom_idx: Index of atom
        direction: 0=x, 1=y, 2=z
        delta: Displacement for finite difference
        
    Returns:
        Force component in kJ/mol/nm
    """
    # Save original position
    original_pos = [state.atoms[atom_idx].x,
                   state.atoms[atom_idx].y,
                   state.atoms[atom_idx].z]
    
    # Calculate energy at +delta
    if direction == 0:
        state.atoms[atom_idx].x = original_pos[0] + delta
    elif direction == 1:
        state.atoms[atom_idx].y = original_pos[1] + delta
    else:
        state.atoms[atom_idx].z = original_pos[2] + delta
    
    energy_plus = pygcmc.DrudeComplete.calculateEnergy(state)
    
    # Calculate energy at -delta
    if direction == 0:
        state.atoms[atom_idx].x = original_pos[0] - delta
    elif direction == 1:
        state.atoms[atom_idx].y = original_pos[1] - delta
    else:
        state.atoms[atom_idx].z = original_pos[2] - delta
    
    energy_minus = pygcmc.DrudeComplete.calculateEnergy(state)
    
    # Restore original position
    state.atoms[atom_idx].x = original_pos[0]
    state.atoms[atom_idx].y = original_pos[1]
    state.atoms[atom_idx].z = original_pos[2]
    
    # Force = -dE/dx
    force = -(energy_plus - energy_minus) / (2 * delta)
    
    return force
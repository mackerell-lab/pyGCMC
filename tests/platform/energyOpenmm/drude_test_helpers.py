"""
Helper functions for Drude tests
"""

import pygcmc


def calculate_numerical_force(state, particle_index, direction, delta=1e-6):
    """Calculate force using finite difference of energy"""
    # Save original position
    original_pos = [state.atoms[particle_index].x, 
                   state.atoms[particle_index].y,
                   state.atoms[particle_index].z]
    
    # Calculate energy at +delta
    if direction == 0:  # x
        state.atoms[particle_index].x = original_pos[0] + delta
    elif direction == 1:  # y
        state.atoms[particle_index].y = original_pos[1] + delta
    else:  # z
        state.atoms[particle_index].z = original_pos[2] + delta
    
    energy_plus = pygcmc.DrudeComplete.calculateEnergy(state)
    
    # Calculate energy at -delta
    if direction == 0:
        state.atoms[particle_index].x = original_pos[0] - delta
    elif direction == 1:
        state.atoms[particle_index].y = original_pos[1] - delta
    else:
        state.atoms[particle_index].z = original_pos[2] - delta
    
    energy_minus = pygcmc.DrudeComplete.calculateEnergy(state)
    
    # Restore original position
    state.atoms[particle_index].x = original_pos[0]
    state.atoms[particle_index].y = original_pos[1]
    state.atoms[particle_index].z = original_pos[2]
    
    # Force = -dE/dx
    force = -(energy_plus - energy_minus) / (2 * delta)
    
    return force
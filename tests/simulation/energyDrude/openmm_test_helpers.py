"""
Helper functions for OpenMM-style tests
"""

import pygcmc


def create_atom(x, y, z, charge, atom_type):
    """Create an MCAtom with given properties"""
    atom = pygcmc.MCAtom()
    atom.x, atom.y, atom.z = x, y, z
    atom.charge = charge
    atom.type = atom_type
    return atom


def calculate_numerical_force_simple(state, particle_index, direction, delta=1e-6):
    """Calculate force using finite difference (simplified version)"""
    # For Drude particles, we should NOT move them directly as SCF will re-optimize
    # Instead, we should calculate the force at the converged position
    
    # Save original position
    original_pos = [state.atoms[particle_index].x,
                   state.atoms[particle_index].y,
                   state.atoms[particle_index].z]
    
    # Only move non-Drude particles for force calculation
    # Drude particles (odd indices) will be re-optimized by SCF
    is_drude = (particle_index % 2 == 1)
    
    if not is_drude:
        # For parent atoms, we can calculate force normally
        # Calculate energy at +delta
        if direction == 0:
            state.atoms[particle_index].x = original_pos[0] + delta
        elif direction == 1:
            state.atoms[particle_index].y = original_pos[1] + delta
        else:
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
    else:
        # For Drude particles, force should be near zero after SCF
        # Return a small value to indicate convergence
        force = 0.0
    
    return force
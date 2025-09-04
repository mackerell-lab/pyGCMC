# tests/simulation/movementInsert/gcmc_insertion_helpers_advanced.py
"""
Advanced GCMC insertion helper functions (continued from basic)
"""

import math
import random
import pygcmc
import numpy as np

# Import basic helpers
from .gcmc_insertion_helpers_basic import (
    create_benzene_forcefield,
    create_water_forcefield
)

# Continue the water molecule creation that was split
def complete_water_molecule(state, x, y, z):
    """Complete water molecule creation (continuation from basic file)"""
    # Oxygen
    o_atom = pygcmc.MCAtom()
    o_atom.x = x
    o_atom.y = y
    o_atom.z = z
    o_atom.charge = -0.834
    o_atom.type = 0
    state.atoms.append(o_atom)
    
    # Hydrogen 1
    h1_atom = pygcmc.MCAtom()
    h1_atom.x = x + 0.0957
    h1_atom.y = y
    h1_atom.z = z
    h1_atom.charge = 0.417
    h1_atom.type = 1
    state.atoms.append(h1_atom)
    
    # Hydrogen 2
    h2_atom = pygcmc.MCAtom()
    h2_atom.x = x - 0.0957 * 0.5
    h2_atom.y = y + 0.0957 * 0.866
    h2_atom.z = z
    h2_atom.charge = 0.417
    h2_atom.type = 1
    state.atoms.append(h2_atom)
    
    # Create residue
    res = pygcmc.MCResidue()
    res.atom_start = len(state.atoms) - 3
    res.atom_count = 3
    res.active = True
    res.fixed = False
    res.type = 0
    state.residues.append(res)
    
    return state


def setup_water_box(nx=3, ny=3, nz=3, spacing=0.5):
    """Setup a box of water molecules for testing"""
    state = pygcmc.MCState()
    
    # Set box dimensions
    box_size = [nx * spacing, ny * spacing, nz * spacing]
    state.info.box = box_size
    state.info.cutoff = min(box_size) / 2.0 - 0.1
    
    # Create water forcefield
    state.forcefield = create_water_forcefield()
    
    # Initialize empty lists
    state.atoms = []
    state.residues = []
    
    # Add water molecules on a grid
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                x = (i + 0.5) * spacing
                y = (j + 0.5) * spacing
                z = (k + 0.5) * spacing
                
                state = complete_water_molecule(state, x, y, z)
    
    state.activeAtomCount = len(state.atoms)
    state.activeResidueCount = len(state.residues)
    
    return state


def calculate_insertion_acceptance_probability(
    state,
    energy_before,
    energy_after,
    n_molecules_before,
    temperature=300.0,
    chemical_potential=-5.0,
    use_cavity_bias=False,
    cavity_volume_fraction=1.0
):
    """
    Calculate GCMC insertion acceptance probability
    
    Parameters:
    -----------
    state : MCState
        System state
    energy_before : float
        System energy before insertion (kJ/mol)
    energy_after : float
        System energy after insertion (kJ/mol)
    n_molecules_before : int
        Number of molecules before insertion
    temperature : float
        Temperature in K
    chemical_potential : float
        Chemical potential in kJ/mol
    use_cavity_bias : bool
        Whether cavity bias was used
    cavity_volume_fraction : float
        Fraction of volume available for insertion
        
    Returns:
    --------
    float : Acceptance probability [0, 1]
    """
    kB = 0.008314463  # kJ/(mol·K)
    beta = 1.0 / (kB * temperature)
    
    # Energy change
    delta_E = energy_after - energy_before
    
    # Volume
    V = state.info.box[0] * state.info.box[1] * state.info.box[2]
    
    # Number of molecules
    N = n_molecules_before
    
    # Standard GCMC acceptance for insertion
    # A = min(1, (V/(N+1)) * exp(β(μ - ΔE)))
    
    # Cavity bias factor
    if use_cavity_bias:
        V_eff = V * cavity_volume_fraction
    else:
        V_eff = V
    
    # Calculate acceptance
    log_acc = math.log(V_eff / (N + 1)) + beta * (chemical_potential - delta_E)
    
    if log_acc > 0:
        return 1.0
    else:
        return math.exp(log_acc)


def calculate_deletion_acceptance_probability(
    state,
    energy_before,
    energy_after,
    n_molecules_before,
    temperature=300.0,
    chemical_potential=-5.0,
    use_cavity_bias=False,
    cavity_volume_fraction=1.0
):
    """
    Calculate GCMC deletion acceptance probability
    
    Similar to insertion but with reversed formula
    """
    kB = 0.008314463  # kJ/(mol·K)
    beta = 1.0 / (kB * temperature)
    
    # Energy change (note: this is E_after - E_before where after has fewer molecules)
    delta_E = energy_after - energy_before
    
    # Volume
    V = state.info.box[0] * state.info.box[1] * state.info.box[2]
    
    # Number of molecules before deletion
    N = n_molecules_before
    
    # Cavity bias factor
    if use_cavity_bias:
        V_eff = V * cavity_volume_fraction
    else:
        V_eff = V
    
    # Calculate acceptance for deletion
    # A = min(1, (N/V) * exp(-β(μ + ΔE)))
    log_acc = math.log(N / V_eff) - beta * (chemical_potential + delta_E)
    
    if log_acc > 0:
        return 1.0
    else:
        return math.exp(log_acc)
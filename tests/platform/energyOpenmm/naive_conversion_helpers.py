# tests/simulation/energyOpenmm/naive_conversion_helpers.py

import pytest
import math
import pygcmc
import os
import warnings

# Filter SWIG-related DeprecationWarning in advance
warnings.filterwarnings("ignore", category=DeprecationWarning,
                        message="builtin type SwigPyPacked has no __module__ attribute")
warnings.filterwarnings("ignore", category=DeprecationWarning,
                        message="builtin type SwigPyObject has no __module__ attribute")
warnings.filterwarnings("ignore", category=DeprecationWarning,
                        message="builtin type swigvarlink has no __module__ attribute")
from openmm import *
from openmm.app import *
from openmm.unit import *
from .naive_system_helpers import create_test_system

def convert_openmm_state_to_mcstate():
    """Convert OpenMM test system to MCState for naive implementation."""
    # Create OpenMM test system
    system, topology, positions = create_test_system()
    
    # Create MCState
    state = pygcmc.MCState()
    
    # Set cutoff distance to 1.0nm, consistent with OpenMM
    state.info.cutoff = 1.0
    
    # 1. Set up force field
    state.forcefield.numTotalTypes = 3  # Three types: C, O, and H
    state.forcefield.numMovementTypes = 1  # C is the movement type
    
    # Get force field parameters from OpenMM system
    nb_force = None
    for force in system.getForces():
        if isinstance(force, NonbondedForce):
            nb_force = force
            break
    
    # Get parameters
    c_params = nb_force.getParticleParameters(0)  # Carbon parameters
    o_params = nb_force.getParticleParameters(6)  # Oxygen parameters
    h_params = nb_force.getParticleParameters(7)  # Hydrogen parameters
    
    # Note: OpenMM's epsilon already includes the factor of 4, so no need to divide by 4
    c_c_eps = c_params[2].value_in_unit(kilojoules_per_mole)
    c_c_sigma = c_params[1].value_in_unit(nanometers)
    o_o_eps = o_params[2].value_in_unit(kilojoules_per_mole)
    o_o_sigma = o_params[1].value_in_unit(nanometers)
    h_h_eps = h_params[2].value_in_unit(kilojoules_per_mole)  # Should be 0
    h_h_sigma = h_params[1].value_in_unit(nanometers)         # Should be 0
    
    # Use OpenMM mixing rules
    def mix_params(eps1, sigma1, eps2, sigma2):
        if eps1 == 0 or eps2 == 0 or sigma1 == 0 or sigma2 == 0:
            return 0.0, 0.0
        # OpenMM mixing rules:
        # - sigma: arithmetic mean 0.5*(sigma1 + sigma2)
        # - epsilon: geometric mean sqrt(eps1 * eps2)
        mixed_sigma = 0.5 * (sigma1 + sigma2)  # Modified here, using 0.5*(sigma1 + sigma2)
        mixed_eps = math.sqrt(eps1 * eps2)
        return mixed_eps, mixed_sigma
    
    # Calculate mixed parameters
    c_o_eps, c_o_sigma = mix_params(c_c_eps, c_c_sigma, o_o_eps, o_o_sigma)
    c_h_eps, c_h_sigma = mix_params(c_c_eps, c_c_sigma, h_h_eps, h_h_sigma)
    o_h_eps, o_h_sigma = mix_params(o_o_eps, o_o_sigma, h_h_eps, h_h_sigma)
    
    # Set force field parameter matrices (numTotalTypes * numTotalTypes = 3 * 3)
    # Complete interaction matrix:
    # [C-C, C-O, C-H]
    # [O-C, O-O, O-H]
    # [H-C, H-O, H-H]
    state.forcefield.ljEps = [
        c_c_eps, c_o_eps, c_h_eps,    # C with (C,O,H) interactions
        c_o_eps, o_o_eps, o_h_eps,    # O with (C,O,H) interactions
        c_h_eps, o_h_eps, h_h_eps     # H with (C,O,H) interactions
    ]
    state.forcefield.ljSigma = [
        c_c_sigma, c_o_sigma, c_h_sigma,    # C with (C,O,H) interactions
        c_o_sigma, o_o_sigma, o_h_sigma,    # O with (C,O,H) interactions
        c_h_sigma, o_h_sigma, h_h_sigma     # H with (C,O,H) interactions
    ]
    
    # Set movement types
    state.movementAtomTypes = [0]  # Type 0 (C) is movement type
    state.numMovementAtomTypes = 1
    
    # 2. Set atoms
    atoms = []
    # Add benzene atoms
    for i in range(6):
        pos = positions[i].value_in_unit(nanometers)
        atom = pygcmc.MCAtom()
        atom.x = pos[0]
        atom.y = pos[1]
        atom.z = pos[2]
        atom.charge = c_params[0].value_in_unit(elementary_charge)
        atom.type = 0  # Carbon type
        atoms.append(atom)
    
    # Add water molecule atoms
    for i in range(6, 9):
        pos = positions[i].value_in_unit(nanometers)
        atom = pygcmc.MCAtom()
        atom.x = pos[0]
        atom.y = pos[1]
        atom.z = pos[2]
        if i == 6:  # Oxygen
            atom.charge = o_params[0].value_in_unit(elementary_charge)
            atom.type = 1  # Oxygen type
        else:  # Hydrogen
            atom.charge = h_params[0].value_in_unit(elementary_charge)
            atom.type = 2  # Hydrogen type (separate from oxygen)
        atoms.append(atom)
    
    state.atoms = atoms
    state.activeAtomCount = len(atoms)
    
    # 3. Set residues
    # Movement residue (benzene)
    movement_res = pygcmc.MCResidue()
    movement_res.active = True
    movement_res.type = 0  # Movement type
    movement_res.atomStart = 0
    movement_res.atomCount = 6
    
    # Fixed residue (water)
    fixed_res = pygcmc.MCResidue()
    fixed_res.active = True
    fixed_res.type = 1  # Fixed type
    fixed_res.atomStart = 6
    fixed_res.atomCount = 3
    
    state.residues = [movement_res, fixed_res]
    state.activeResidueCount = 2
    
    # 4. Set movement residue info
    movement_info = pygcmc.MCMovementResidueInfo()
    movement_info.startIndex = 0
    movement_info.activeCount = 1
    movement_info.totalCount = 1
    
    state.movementResidues = [movement_info]
    
    # 5. Set box and cutoff
    box_vectors = system.getDefaultPeriodicBoxVectors()
    state.info.box = [
        box_vectors[0][0].value_in_unit(nanometers),
        box_vectors[1][1].value_in_unit(nanometers),
        box_vectors[2][2].value_in_unit(nanometers)
    ]
    state.info.cutoff = nb_force.getCutoffDistance().value_in_unit(nanometers)
    
    return state, system, positions

def print_force_field_params():
    """Print force field parameters for both implementations."""
    state, system, positions = convert_openmm_state_to_mcstate()
    
    # Print OpenMM parameters
    nb_force = None
    for force in system.getForces():
        if isinstance(force, NonbondedForce):
            nb_force = force
            break
            
    print("\nOpenMM parameters:")
    for i in range(nb_force.getNumParticles()):
        charge, sigma, epsilon = nb_force.getParticleParameters(i)
        print(f"Atom {i}: q={charge.value_in_unit(elementary_charge):.3f}e, "
              f"sigma={sigma.value_in_unit(nanometers):.3f}nm, "
              f"epsilon={epsilon.value_in_unit(kilojoules_per_mole):.3f}kJ/mol")
        
    print("\nNaive implementation parameters:")
    print("LJ Epsilon matrix [kJ/mol]:")
    n = int(math.sqrt(len(state.forcefield.ljEps)))
    for i in range(n):
        row = state.forcefield.ljEps[i*n:(i+1)*n]
        print(f"Type {i}: {[f'{x:.3f}' for x in row]}")
    
    print("\nLJ Sigma matrix [nm]:")
    for i in range(n):
        row = state.forcefield.ljSigma[i*n:(i+1)*n]
        print(f"Type {i}: {[f'{x:.3f}' for x in row]}")


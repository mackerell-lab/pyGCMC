# tests/simulation/energyOpenmm/nbfix_openmm_helpers.py
"""
OpenMM helper functions for NBFIX tests

This module provides OpenMM system setup functions for NBFIX testing.
"""

import pytest
import math
import pygcmc
import os
import warnings

# Filter SWIG-related warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)

try:
    from openmm import *
    from openmm.app import *
    from openmm.unit import *
    HAS_OPENMM = True
except ImportError:
    HAS_OPENMM = False

# Constants
ANGSTROM_TO_NM = 0.1
KCAL_TO_KJ = 4.184
kC = 138.935456  # Coulomb constant in kJ·nm/mol/e^2

def setup_openmm_system_with_nbfix(mc_state, forcefield, separate_forces=False):
    """Create OpenMM system with NBFIX parameters matching MCState
    
    Args:
        mc_state: PyGCMC MCState object
        forcefield: PyGCMC ForceField object with NBFIX parameters
        separate_forces: If True, create separate VDW and electrostatic forces for detailed comparison
    
    Returns:
        system: OpenMM System object
        positions: List of atomic positions
    """
    system = System()
    
    # Add particles with masses (use standard masses for now)
    masses = {"SOD": 22.98977, "OT": 15.9994, "HT": 1.008, "CLA": 35.45}
    for i in range(mc_state.activeAtomCount):
        atom_type = mc_state.atomTypes.atomTypes[mc_state.atoms[i].type]
        mass = masses.get(atom_type, 1.0)
        system.addParticle(mass)
    
    # Set box
    box = mc_state.info.box
    system.setDefaultPeriodicBoxVectors(
        Vec3(box[0], 0, 0) * nanometers,
        Vec3(0, box[1], 0) * nanometers,
        Vec3(0, 0, box[2]) * nanometers
    )
    
    # Build parameter tables from MCState forcefield
    n_types = mc_state.forcefield.numTotalTypes
    sigma_table = []
    epsilon_table = []
    
    # Copy parameters from MCState forcefield matrix
    for idx in range(n_types * n_types):
        sigma_table.append(mc_state.forcefield.ljSigma[idx])
        epsilon_table.append(mc_state.forcefield.ljEps[idx])
    
    if separate_forces:
        # Create separate VDW and electrostatic forces for detailed energy comparison
        
        # VDW Force (Force Group 1)
        vdw_expr = """4*epsilon*((sigma/r)^12-(sigma/r)^6);
                      sigma=sigma_table(type1, type2);
                      epsilon=epsilon_table(type1, type2)"""
        vdw_force = CustomNonbondedForce(vdw_expr)
        vdw_force.addPerParticleParameter("type")
        
        # Add tables for VDW
        vdw_force.addTabulatedFunction("sigma_table", 
                                      Discrete2DFunction(n_types, n_types, sigma_table))
        vdw_force.addTabulatedFunction("epsilon_table", 
                                      Discrete2DFunction(n_types, n_types, epsilon_table))
        
        # Add particles to VDW force
        for i in range(mc_state.activeAtomCount):
            vdw_force.addParticle([float(mc_state.atoms[i].type)])
        
        # Set VDW cutoff and method
        vdw_force.setNonbondedMethod(CustomNonbondedForce.CutoffPeriodic)
        vdw_force.setCutoffDistance(mc_state.info.cutoff * nanometers)
        vdw_force.setForceGroup(1)
        
        # Electrostatic Force (Force Group 2)
        elec_expr = "kC*q1*q2/r"
        elec_force = CustomNonbondedForce(elec_expr)
        elec_force.addGlobalParameter("kC", kC)  # kC = 138.935456 kJ·nm/mol/e²
        elec_force.addPerParticleParameter("q")
        
        # Add particles to electrostatic force
        for i in range(mc_state.activeAtomCount):
            elec_force.addParticle([mc_state.atoms[i].charge])
        
        # Set electrostatic cutoff and method
        elec_force.setNonbondedMethod(CustomNonbondedForce.CutoffPeriodic)
        elec_force.setCutoffDistance(mc_state.info.cutoff * nanometers)
        elec_force.setForceGroup(2)
        
        # Add exclusions to BOTH forces for bonded atoms in water
        # Water has atoms at indices 1, 2, 3 (OT, HT, HT)
        vdw_force.addExclusion(1, 2)   # O-H1
        vdw_force.addExclusion(1, 3)   # O-H2
        vdw_force.addExclusion(2, 3)   # H1-H2
        elec_force.addExclusion(1, 2)  # O-H1
        elec_force.addExclusion(1, 3)  # O-H2
        elec_force.addExclusion(2, 3)  # H1-H2
        
        system.addForce(vdw_force)
        system.addForce(elec_force)
    else:
        # Create combined force (original behavior)
        energy_expr = """4*epsilon*((sigma/r)^12-(sigma/r)^6) + kC*q1*q2/r;
                        sigma=sigma_table(type1, type2);
                        epsilon=epsilon_table(type1, type2)"""
        
        nbforce = CustomNonbondedForce(energy_expr)
        nbforce.addGlobalParameter("kC", kC)  # kC = 138.935456 kJ·nm/mol/e²
        nbforce.addPerParticleParameter("q")
        nbforce.addPerParticleParameter("type")
        
        # Add tables to force
        nbforce.addTabulatedFunction("sigma_table", 
                                    Discrete2DFunction(n_types, n_types, sigma_table))
        nbforce.addTabulatedFunction("epsilon_table", 
                                    Discrete2DFunction(n_types, n_types, epsilon_table))
        
        # Add particles with parameters - ensure type is float to avoid interpolation issues
        for i in range(mc_state.activeAtomCount):
            atom = mc_state.atoms[i]
            nbforce.addParticle([atom.charge, float(atom.type)])
        
        # Set method and cutoff
        nbforce.setNonbondedMethod(CustomNonbondedForce.CutoffPeriodic)
        nbforce.setCutoffDistance(mc_state.info.cutoff * nanometers)
        
        # Add exclusions for bonded atoms in water (TIP3P has 3 atoms)
        nbforce.addExclusion(1, 2)  # O-H1
        nbforce.addExclusion(1, 3)  # O-H2
        nbforce.addExclusion(2, 3)  # H1-H2
        
        system.addForce(nbforce)
    
    # Create positions
    positions = []
    for i in range(mc_state.activeAtomCount):
        atom = mc_state.atoms[i]
        positions.append(Vec3(atom.x, atom.y, atom.z) * nanometers)
    
    return system, positions


# tests/simulation/movementInsert/advanced_insertion_systems.py
"""
System creation functions for advanced insertion tests

This module provides various system creation utilities for testing
advanced molecule insertion strategies.
"""

import pytest
import random
import math
import pygcmc
# Helper functions

def create_empty_system():
    """Create an empty system for insertions"""
    state = pygcmc.MCState()
    state.info.box = [5.0, 5.0, 5.0]
    state.info.cutoff = 2.5
    
    # Define atom types
    state.atomTypes.get_or_add_type("FRAMEWORK")
    state.atomTypes.get_or_add_type("GUEST")
    state.atomTypes.get_or_add_type("OT")
    state.atomTypes.get_or_add_type("HT")
    
    # Set up force field
    state.forcefield = create_mixed_forcefield()
    
    state.atoms = []
    state.residues = []
    state.activeAtomCount = 0
    state.activeResidueCount = 0
    
    return state


def create_cavity_system():
    """Create a system with a central cavity"""
    system = create_empty_system()
    
    # Pre-add GUEST atom type to avoid forcefield mismatch
    system.atomTypes.get_or_add_type("GUEST")
    
    # Update forcefield to handle 5 atom types
    system.forcefield.numTotalTypes = 5
    system.forcefield.numMovementTypes = 5
    
    # Extend LJ parameters to 5x5 matrix = 25 values
    # Keep existing 4x4 and add guest interactions
    old_eps = list(system.forcefield.ljEps)
    old_sigma = list(system.forcefield.ljSigma)
    
    # Create new 5x5 arrays
    new_eps = []
    new_sigma = []
    
    # Check if we have the expected 16 values
    if len(old_eps) != 16 or len(old_sigma) != 16:
        # If not enough parameters, just create simple 5x5 matrix
        for i in range(25):
            new_eps.append(0.5)
            new_sigma.append(0.3)
    else:
        # Copy old 4x4 values and add guest column for each row
        for i in range(4):
            for j in range(4):
                new_eps.append(old_eps[i*4 + j])
                new_sigma.append(old_sigma[i*4 + j])
            # Add guest interaction for this row
            new_eps.append(0.5)  # Guest interaction
            new_sigma.append(0.3)
        
        # Add guest row
        for i in range(5):
            new_eps.append(0.5)  # Guest interactions
            new_sigma.append(0.3)
    
    system.forcefield.ljEps = new_eps
    system.forcefield.ljSigma = new_sigma
    
    # Create framework atoms around cavity
    framework_positions = [
        # Corners of a cube with cavity in center
        (1.0, 1.0, 1.0), (1.0, 1.0, 4.0),
        (1.0, 4.0, 1.0), (1.0, 4.0, 4.0),
        (4.0, 1.0, 1.0), (4.0, 1.0, 4.0),
        (4.0, 4.0, 1.0), (4.0, 4.0, 4.0),
    ]
    
    for x, y, z in framework_positions:
        atom = pygcmc.MCAtom()
        atom.x = x
        atom.y = y
        atom.z = z
        atom.charge = -0.2
        atom.type = 0  # FRAMEWORK (OT type)
        
        res = pygcmc.MCResidue()
        res.active = True
        res.atomStart = len(system.atoms)
        res.atomCount = 1
        res.type = 1  # Framework type
        
        system.atoms.append(atom)
        system.residues.append(res)
    
    system.activeAtomCount = len(system.atoms)
    system.activeResidueCount = len(system.residues)
    
    return system


def create_multi_cavity_system():
    """Create system with multiple cavities of different sizes"""
    system = create_empty_system()
    
    # Small cavity (tight framework)
    small_cavity_atoms = [
        (1.0, 1.0, 2.5), (1.0, 2.0, 2.5),
        (2.0, 1.0, 2.5), (2.0, 2.0, 2.5),
    ]
    
    # Large cavity (spacious)
    large_cavity_atoms = [
        (3.0, 3.0, 2.0), (3.0, 4.0, 2.0),
        (4.0, 3.0, 2.0), (4.0, 4.0, 2.0),
        (3.0, 3.0, 3.0), (3.0, 4.0, 3.0),
        (4.0, 3.0, 3.0), (4.0, 4.0, 3.0),
    ]
    
    # Edge cavity
    edge_cavity_atoms = [
        (2.0, 2.0, 0.5), (2.0, 3.0, 0.5),
        (3.0, 2.0, 0.5), (3.0, 3.0, 0.5),
    ]
    
    all_positions = small_cavity_atoms + large_cavity_atoms + edge_cavity_atoms
    
    for x, y, z in all_positions:
        atom = pygcmc.MCAtom()
        atom.x = x
        atom.y = y
        atom.z = z
        atom.charge = -0.1
        atom.type = 0  # FRAMEWORK
        
        res = pygcmc.MCResidue()
        res.active = True
        res.atomStart = len(system.atoms)
        res.atomCount = 1
        res.type = 1  # Framework type
        
        system.atoms.append(atom)
        system.residues.append(res)
    
    system.activeAtomCount = len(system.atoms)
    system.activeResidueCount = len(system.residues)
    
    return system


def create_mixed_system():
    """Create system with framework and existing guests"""
    system = create_cavity_system()
    
    # Add some existing guest molecules
    existing_guests = [
        (1.5, 1.5, 2.5, -0.5),
        (3.5, 3.5, 2.5, 0.5),
    ]
    
    for x, y, z, charge in existing_guests:
        atom = pygcmc.MCAtom()
        atom.x = x
        atom.y = y
        atom.z = z
        atom.charge = charge
        atom.type = 1  # GUEST
        
        res = pygcmc.MCResidue()
        res.active = True
        res.atomStart = len(system.atoms)
        res.atomCount = 1
        res.type = 0  # Guest type
        
        system.atoms.append(atom)
        system.residues.append(res)
    
    system.activeAtomCount = len(system.atoms)
    system.activeResidueCount = len(system.residues)
    
    return system


def create_ranked_cavity_system():
    """Create system with cavities of different favorability"""
    system = create_empty_system()
    
    # Create framework with asymmetric charge distribution
    framework_configs = [
        # Best cavity - surrounded by attractive charges
        [(1.0, 2.0, 2.5, -0.3), (1.0, 3.0, 2.5, -0.3),
         (2.0, 2.0, 2.5, -0.3), (2.0, 3.0, 2.5, -0.3)],
        
        # Second best - mixed charges
        [(3.0, 2.0, 2.5, -0.2), (3.0, 3.0, 2.5, -0.1),
         (4.0, 2.0, 2.5, -0.1), (4.0, 3.0, 2.5, -0.2)],
        
        # Third - mostly neutral
        [(2.0, 1.0, 2.5, -0.05), (3.0, 1.0, 2.5, -0.05),
         (2.0, 1.0, 3.0, -0.05), (3.0, 1.0, 3.0, -0.05)],
        
        # Fourth - some repulsion
        [(2.0, 4.0, 2.5, 0.1), (3.0, 4.0, 2.5, -0.1),
         (2.0, 4.0, 3.0, -0.1), (3.0, 4.0, 3.0, 0.1)],
    ]
    
    for cavity_atoms in framework_configs:
        for x, y, z, charge in cavity_atoms:
            atom = pygcmc.MCAtom()
            atom.x = x
            atom.y = y
            atom.z = z
            atom.charge = charge
            atom.type = 0  # FRAMEWORK
            
            res = pygcmc.MCResidue()
            res.active = True
            res.atomStart = len(system.atoms)
            res.atomCount = 1
            res.type = 1  # Framework type
            
            system.atoms.append(atom)
            system.residues.append(res)
    
    system.activeAtomCount = len(system.atoms)
    system.activeResidueCount = len(system.residues)
    
    return system


def create_mixed_forcefield():
    """Create force field for framework/guest/water system"""
    ff = pygcmc.MCForceField()
    
    # 4 types: FRAMEWORK, GUEST, OT, HT
    ff.numTotalTypes = 4
    ff.numMovementTypes = 3  # GUEST, OT, HT can move
    
    # LJ parameters - need full 4x4 matrix = 16 values
    # Simple parameters for testing
    ff.ljEps = [
        0.5, 0.5, 0.5, 0.0,    # FRAMEWORK with all
        0.5, 0.5, 0.5, 0.0,    # GUEST with all
        0.5, 0.5, 0.6364, 0.0, # OT (water oxygen) with all
        0.0, 0.0, 0.0, 0.0     # HT (water hydrogen) with all
    ]
    ff.ljSigma = [
        0.35, 0.35, 0.35, 0.0,    # FRAMEWORK
        0.35, 0.35, 0.35, 0.0,    # GUEST
        0.35, 0.35, 0.3166, 0.0,  # OT
        0.0, 0.0, 0.0, 0.0        # HT
    ]
    
    return ff


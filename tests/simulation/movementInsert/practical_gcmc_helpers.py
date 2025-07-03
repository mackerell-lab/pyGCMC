# tests/simulation/movementInsert/practical_gcmc_helpers.py
"""
Helper functions for practical GCMC insertion tests

This module provides utility functions for cavity finding, molecule creation,
and energy calculations.
"""

import pytest
import math
import random
import numpy as np
import pygcmc
import os

# Constants
kB = 0.008314463  # Boltzmann constant in kJ/mol/K
kC = 138.935456   # Coulomb constant in kJ·nm/mol/e²


def find_best_cavity_position(system):
    """Find best cavity position by grid search"""
    best_pos = None
    best_score = -float('inf')
    
    # Grid search in cavity region
    for _ in range(50):
        x = random.uniform(2.0, 3.0)
        y = random.uniform(2.0, 3.0)
        z = random.uniform(2.0, 3.0)
        
        # Score based on distance to protein atoms
        score = 0.0
        for atom in system.atoms:
            if atom.type > 0:  # Protein atom
                dist = math.sqrt((x-atom.x)**2 + (y-atom.y)**2 + (z-atom.z)**2)
                if dist < 0.25:  # Too close
                    score = -float('inf')
                    break
                elif dist < 0.5:  # Good distance for VDW
                    score += 1.0
                else:
                    score += 0.1 / dist
        
        if score > best_score:
            best_score = score
            best_pos = (x, y, z)
    
    return best_pos if best_score > 0 else None


def create_oriented_benzene(x, y, z, theta, phi):
    """Create benzene molecule with specified position and orientation"""
    # Benzene ring coordinates (centered at origin)
    ring_coords = []
    for i in range(6):
        angle = i * math.pi / 3
        rx = 0.14 * math.cos(angle)  # 1.4 Å radius
        ry = 0.14 * math.sin(angle)
        ring_coords.append((rx, ry, 0.0))
    
    # Apply rotation
    atoms = []
    for rx, ry, rz in ring_coords:
        # Rotate around z-axis (theta)
        x1 = rx * math.cos(theta) - ry * math.sin(theta)
        y1 = rx * math.sin(theta) + ry * math.cos(theta)
        z1 = rz
        
        # Rotate around y-axis (phi)
        x2 = x1 * math.cos(phi) + z1 * math.sin(phi)
        y2 = y1
        z2 = -x1 * math.sin(phi) + z1 * math.cos(phi)
        
        # Translate to position
        atom = pygcmc.MCAtom()
        atom.x = x + x2
        atom.y = y + y2
        atom.z = z + z2
        atom.charge = -0.115  # Benzene carbon charge
        atom.type = 0  # Guest type
        atoms.append(atom)
    
    return atoms


def create_system_with_guest(original_system, guest_atoms):
    """Create new system with guest molecule added"""
    new_state = pygcmc.MCState()
    new_state.info.box = original_system.info.box
    new_state.info.cutoff = original_system.info.cutoff
    new_state.atomTypes = original_system.atomTypes
    new_state.forcefield = original_system.forcefield
    
    # Copy atoms
    new_atoms = []
    for atom in original_system.atoms:
        new_atom = pygcmc.MCAtom()
        new_atom.x = atom.x
        new_atom.y = atom.y
        new_atom.z = atom.z
        new_atom.charge = atom.charge
        new_atom.type = atom.type
        new_atoms.append(new_atom)
    
    # Add guest atoms
    guest_start = len(new_atoms)
    for atom in guest_atoms:
        new_atoms.append(atom)
    
    # Copy residues
    new_residues = []
    for res in original_system.residues:
        new_res = pygcmc.MCResidue()
        new_res.active = res.active
        new_res.atomStart = res.atomStart
        new_res.atomCount = res.atomCount
        new_res.type = res.type
        new_residues.append(new_res)
    
    # Add guest residue
    guest_res = pygcmc.MCResidue()
    guest_res.active = True
    guest_res.atomStart = guest_start
    guest_res.atomCount = len(guest_atoms)
    guest_res.type = 0  # Guest type
    new_residues.append(guest_res)
    
    new_state.atoms = new_atoms
    new_state.residues = new_residues
    new_state.activeAtomCount = len(new_atoms)
    new_state.activeResidueCount = len(new_residues)
    
    return new_state


def calculate_total_energy(state):
    """Calculate total system energy"""
    total = 0.0
    for res in state.residues:
        if res.active:
            total += res.energy_vdw + res.energy_elec
    return total / 2.0  # Correct for double counting


def count_guest_molecules(system):
    """Count guest molecules (type 0 residues)"""
    return sum(1 for res in system.residues if res.type == 0)


def estimate_cavity_fraction(system):
    """Estimate fraction of unoccupied volume"""
    # Simplified estimation based on number of atoms
    n_atoms = len(system.atoms)
    box_volume = system.info.box[0] * system.info.box[1] * system.info.box[2]
    atom_volume = n_atoms * 0.05  # Approximate volume per atom
    return max(0.1, 1.0 - atom_volume / box_volume)


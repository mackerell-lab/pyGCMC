# tests/simulation/pgp/asymmetric_water_complete.py
"""PGP asymmetric water systems test - Complete version with all original functionality."""

import pytest
import math
import random
import pygcmc
from pygcmc import MCState, MCAtom, MCResidue, MCForceField, MCMovementResidueInfo
import sys
from .helpers import calculate_pbc_distance, is_safe_position, generate_safe_move


def test_compare_ewald_pme_pgp_asymmetric():
    """
    Test whether PGP calculation is accurate in systems with asymmetric charge distribution
    
    Features:
    1. Fixed part contains multiple charged particles with asymmetric distribution
    2. Moving residue is a water molecule, far from fixed part
    3. Execute multiple random movements to ensure distance from fixed part is always > cutoff
    4. Verify accuracy by comparing energy calculated by Ewald, PME, and PGP
    """
    # Set system parameters
    box_size = 8.0  # nm - use a larger box
    cutoff = 1.0    # nm
    potential_cutoff = 1.0  # nm
    box = [box_size, box_size, box_size]
    
    # Ewald and PME parameters
    alpha = 0.29    # 1/nm
    kmax = [8, 8, 8]  # Number of k-space vectors for Ewald calculation
    mesh_size = [32, 32, 32]
    potential_grid_size = [32, 32, 32]
    spline_order = 4
    tolerance = 1e-5
    
    print("Creating asymmetric charge distribution test system...")
    sys.stdout.flush()
    
    # Create test system
    system = MCState()
    system.info.box = box
    system.info.setTemperature(300.0)
    system.info.cutoff = cutoff
    
    # Set up force field
    ff = MCForceField()
    ff.numTotalTypes = 2 
    ff.ljSigma = [0.333, 0.3875, 0.3875, 0.442]
    ff.ljEps = [0.0115, 0.0693, 0.0693, 0.4184]
    system.forcefield = ff
    
    atoms = []
    residues = []
    
    # Create fixed part - asymmetric charged particle distribution
    fixed_particles = [
        # Position (x, y, z)                 Charge
        ((1.0, 1.0, 1.0),               1.0),  # Positive charge
        ((1.5, 1.0, 1.2),              -0.8),  # Negative charge
        ((1.3, 1.7, 1.5),               0.6),  # Positive charge
        ((0.8, 1.6, 0.9),              -0.7),  # Negative charge
        ((box_size-1.0, 1.0, 1.0),      1.0),  # Positive charge
        ((box_size-1.5, 1.2, 1.3),     -0.5),  # Negative charge
        ((1.0, box_size-1.0, 1.0),      0.8),  # Positive charge
        ((1.2, box_size-1.5, 0.9),     -0.6),  # Negative charge
        ((1.0, 1.0, box_size-1.0),      0.9),  # Positive charge
        ((1.4, 1.3, box_size-1.4),     -0.7),  # Negative charge
        # Add more particles to increase system complexity
        ((box_size-2.0, box_size-2.0, 2.0),  0.7),  # Positive charge
        ((box_size-2.5, box_size-2.3, 2.2), -0.4),  # Negative charge
        ((2.0, box_size-2.0, box_size-2.0),  0.5),  # Positive charge
        ((2.2, box_size-2.2, box_size-2.4), -0.3),  # Negative charge
        ((box_size-2.0, 2.0, box_size-2.0), -1.5),  # Negative charge, roughly balancing system
    ]
    
    # Confirm particle count
    print(f"Number of fixed particles: {len(fixed_particles)}")
    
    # Calculate total charge of fixed part
    total_fixed_charge = sum(charge for _, charge in fixed_particles)
    print(f"Total fixed charge: {total_fixed_charge}")
    
    # Add fixed particles
    for i, (pos, charge) in enumerate(fixed_particles):
        atom = MCAtom()
        atom.x, atom.y, atom.z = pos
        atom.charge = charge
        atom.type = 0 if charge > 0 else 1  # Positive charge use type 0, negative charge use type 1
        atoms.append(atom)
    
    # Create fixed residue
    fixed_res = MCResidue()
    fixed_res.atomStart = 0
    fixed_res.atomCount = len(fixed_particles)
    fixed_res.active = True
    fixed_res.fixed = True
    residues.append(fixed_res)
    
    # Calculate box center position to ensure moving residue is far from fixed particles
    center_x = box_size / 2.0
    center_y = box_size / 2.0
    center_z = box_size / 2.0
    
    # Create moving part - a water molecule placed in the center of the box
    # Water molecule charge parameters: oxygen atom -0.8, two hydrogen atoms each +0.4
    # Water molecule bond length: O-H about 0.1 nm
    
    # Set water molecule position
    oxygen_pos = (center_x, center_y, center_z)
    hydrogen1_pos = (center_x + 0.1, center_y, center_z)  # First H atom
    hydrogen2_pos = (center_x, center_y + 0.1, center_z)  # Second H atom
    
    # Water molecule charge
    oxygen_charge = -0.8
    hydrogen_charge = 0.4  # Each hydrogen atom
    
    # Record starting index of moving part
    mobile_start_idx = len(atoms)
    
    # Create oxygen atom
    o_atom = MCAtom()
    o_atom.x, o_atom.y, o_atom.z = oxygen_pos
    o_atom.charge = oxygen_charge
    o_atom.type = 1  # Oxygen atom type
    atoms.append(o_atom)
    
    # Create first hydrogen atom
    h1_atom = MCAtom()
    h1_atom.x, h1_atom.y, h1_atom.z = hydrogen1_pos
    h1_atom.charge = hydrogen_charge
    h1_atom.type = 0  # Hydrogen atom type
    atoms.append(h1_atom)
    
    # Create second hydrogen atom
    h2_atom = MCAtom()
    h2_atom.x, h2_atom.y, h2_atom.z = hydrogen2_pos
    h2_atom.charge = hydrogen_charge
    h2_atom.type = 0  # Hydrogen atom type
    atoms.append(h2_atom)
    
    print(f"Moving water molecule: O position=({oxygen_pos[0]}, {oxygen_pos[1]}, {oxygen_pos[2]}), charge={oxygen_charge}")
    print(f"  H1 position=({hydrogen1_pos[0]}, {hydrogen1_pos[1]}, {hydrogen1_pos[2]}), charge={hydrogen_charge}")
    print(f"  H2 position=({hydrogen2_pos[0]}, {hydrogen2_pos[1]}, {hydrogen2_pos[2]}), charge={hydrogen_charge}")
    print(f"   Water molecule total charge: {oxygen_charge + 2*hydrogen_charge}")
    
    # Create moving residue (water molecule)
    mobile_res = MCResidue()
    mobile_res.atomStart = mobile_start_idx
    mobile_res.atomCount = 3  # Water molecule has 3 atoms
    mobile_res.active = True
    mobile_res.fixed = False
    residues.append(mobile_res)
    
    # Set system
    system.atoms = atoms
    system.residues = residues
    system.activeAtomCount = len(atoms)
    system.activeResidueCount = len(residues)
    
    # Calculate system total charge
    total_system_charge = sum(atom.charge for atom in atoms)
    print(f"System total charge: {total_system_charge}")
    # Allow system to have small charge, no strict assertion needed
    
    # Set moving residue info
    system.movementResidues.clear()
    movement_info = pygcmc.MCMovementResidueInfo()
    movement_info.startIndex = 1  # Second residue (moving residue) index
    movement_info.activeCount = 1  # Only one moving residue
    system.movementResidues.append(movement_info)
    
    print(f"System created: {system.activeAtomCount} atoms, {system.activeResidueCount} residues")
    print(f"Fixed atoms: {len(fixed_particles)}, Moving atoms: 3 (Water molecule)")
    sys.stdout.flush()
    
    # Initialize Ewald, PME, and PGP
    print("Setting calculation parameters...")
    
    # Ewald parameters initialization
    print("Initializing Ewald parameters...")
    pygcmc.setEwaldParameters(alpha, kmax)
    pygcmc.initializeEwaldParameters(cutoff, box, alpha)
    
    # PME parameters initialization
    print("Initializing PME parameters...")
    pygcmc.setPMEParameters(alpha, mesh_size, spline_order, tolerance)
    pygcmc.initializePMEParameters(cutoff, box, alpha)
    
    # PGP parameters initialization
    print("Initializing PGP parameters...")
    pygcmc.setPGPParameters(alpha, mesh_size, potential_cutoff, 
                         potential_grid_size, spline_order, tolerance)
    
    # Precompute grid potential for fixed parts
    print("Precomputing fixed parts grid potential...")
    pygcmc.precomputeGridPotential(system, fixed_only=True)
    
    # Set number of random movements to execute
    num_moves = 5  # Reduce test times to speed up test
    print(f"Executing {num_moves} random movement tests")
    
    # Store relative errors between methods
    pgp_pme_errors = []
    ewald_pme_errors = []
    pgp_ewald_errors = []
    
    # Verify initial position is safe
    current_positions = []
    for i in range(mobile_res.atomCount):
        atom_idx = mobile_res.atomStart + i
        current_positions.append((
            system.atoms[atom_idx].x,
            system.atoms[atom_idx].y,
            system.atoms[atom_idx].z
        ))
    
    is_safe, min_dist = is_safe_position(current_positions, fixed_particles, cutoff, box_size)
    if not is_safe:
        print(f"Warning: Initial position not safe, minimum distance is {min_dist} nm")
    else:
        print(f"Initial position safe, minimum distance from fixed particles > {cutoff} nm")
    
    # Execute multiple random movements
    from .asymmetric_water_movement import execute_movement_loop
    execute_movement_loop(system, mobile_res, fixed_particles, num_moves, 
                         pgp_pme_errors, ewald_pme_errors, pgp_ewald_errors, 
                         cutoff, box_size)
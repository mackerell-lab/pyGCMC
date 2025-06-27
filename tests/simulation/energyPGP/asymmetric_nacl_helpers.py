# tests/simulation/energyPGP/asymmetric_nacl_helpers.py
"""Helper functions for asymmetric NaCl test with complex movement algorithm."""

import math
import pygcmc
from .helpers import calculate_pbc_distance, is_safe_position


def create_nacl_crystal_system(n_cells, a, box_size, cutoff):
    """Create NaCl crystal system with moving residue."""
    from pygcmc import MCState, MCAtom, MCResidue, MCForceField, MCMovementResidueInfo
    
    # Create test system
    system = MCState()
    system.info.box = [box_size, box_size, box_size]
    system.info.setTemperature(300.0)
    system.info.cutoff = cutoff
    
    # Set force field parameters
    ff = MCForceField()
    ff.numTotalTypes = 2  # Na+ and Cl-
    
    # LJ parameters (from OPLS-AA force field)
    sigma_na = 0.333  # nm
    sigma_cl = 0.442  # nm
    eps_na = 0.0115  # kJ/mol
    eps_cl = 0.4184  # kJ/mol
    
    # Set LJ parameter matrix
    ff.ljSigma = [
        sigma_na, (sigma_na + sigma_cl)/2.0,
        (sigma_na + sigma_cl)/2.0, sigma_cl
    ]
    ff.ljEps = [
        eps_na, math.sqrt(eps_na * eps_cl),
        math.sqrt(eps_na * eps_cl), eps_cl
    ]
    
    system.forcefield = ff
    
    # Initialize NaCl lattice
    atoms = []
    residues = []
    fixed_positions = []
    
    # Calculate offset to center the crystal in the box
    offset = (box_size - n_cells * a) / 2.0
    print(f"Crystal centered in box with offset: {offset} nm")
    
    for i in range(n_cells):
        for j in range(n_cells):
            for k in range(n_cells):
                # Na+ ion position
                na_pos = (i * a + offset, j * a + offset, k * a + offset)
                
                # Cl- ion position
                cl_pos = (i * a + a/2 + offset, j * a + a/2 + offset, k * a + a/2 + offset)
                
                # Create residue for each ion pair (Na+, Cl-)
                if i == n_cells-1 and j == n_cells-1 and k == n_cells-1:
                    # Last ion pair will be the moving residue
                    # Add as separate residue later
                    continue
                else:
                    # Add Na+ ion
                    na = MCAtom()
                    na.x, na.y, na.z = na_pos
                    na.charge = 1.0
                    na.type = 0
                    atoms.append(na)
                    fixed_positions.append((na_pos, 1.0))  # Include charge
                    
                    # Add Cl- ion
                    cl = MCAtom()
                    cl.x, cl.y, cl.z = cl_pos
                    cl.charge = -1.0
                    cl.type = 1
                    atoms.append(cl)
                    fixed_positions.append((cl_pos, -1.0))  # Include charge
                    
                    # Create a residue for this ion pair
                    res = MCResidue()
                    res.atomStart = len(atoms) - 2
                    res.atomCount = 2
                    res.active = True
                    res.fixed = True
                    residues.append(res)
    
    # Create moving NaCl residue (the last ion pair)
    mobile_start_idx = len(atoms)
    
    # Last Na+ ion at the corner of the supercell
    na_pos = ((n_cells-1) * a + offset, (n_cells-1) * a + offset, (n_cells-1) * a + offset)
    
    # Last Cl- ion at the corresponding position
    cl_pos = ((n_cells-1) * a + a/2 + offset, (n_cells-1) * a + a/2 + offset, (n_cells-1) * a + a/2 + offset)
    
    # Add Na+ ion for moving residue
    na = MCAtom()
    na.x, na.y, na.z = na_pos
    na.charge = 1.0
    na.type = 0
    atoms.append(na)
    
    # Add Cl- ion for moving residue
    cl = MCAtom()
    cl.x, cl.y, cl.z = cl_pos
    cl.charge = -1.0
    cl.type = 1
    atoms.append(cl)
    
    # Create moving residue
    mobile_res = MCResidue()
    mobile_res.atomStart = mobile_start_idx
    mobile_res.atomCount = 2  # NaCl has 2 atoms
    mobile_res.active = True
    mobile_res.fixed = False
    residues.append(mobile_res)
    
    print(f"Moving NaCl residue: Na+ position=({na_pos[0]:.3f}, {na_pos[1]:.3f}, {na_pos[2]:.3f}), charge=1.0")
    print(f"  Cl- position=({cl_pos[0]:.3f}, {cl_pos[1]:.3f}, {cl_pos[2]:.3f}), charge=-1.0")
    print(f"  NaCl total charge: 0.0")
    
    # Set system
    system.atoms = atoms
    system.residues = residues
    system.activeAtomCount = len(atoms)
    system.activeResidueCount = len(residues)
    
    # Set moving residue info
    system.movementResidues.clear()
    movement_info = pygcmc.MCMovementResidueInfo()
    movement_info.startIndex = len(residues) - 1  # Last residue is the moving residue
    movement_info.activeCount = 1  # Only one moving residue
    system.movementResidues.append(movement_info)
    
    return system, mobile_res, fixed_positions


def handle_initial_position_safety(system, mobile_res, fixed_positions, cutoff, box_size):
    """Handle initial position safety checks and corrections."""
    # Get initial mobile atom positions
    current_positions = []
    for i in range(mobile_res.atomCount):
        atom_idx = mobile_res.atomStart + i
        current_positions.append((
            system.atoms[atom_idx].x,
            system.atoms[atom_idx].y,
            system.atoms[atom_idx].z
        ))
    
    # Check if initial position is safe
    is_safe, min_dist = is_safe_position(current_positions, fixed_positions, cutoff, box_size)
    if not is_safe:
        print(f"Warning: Initial position not safe, minimum distance is {min_dist} nm")
        print("Moving mobile ions to center to avoid immediate conflicts...")
        
        # Move mobile residue to box center to avoid conflicts
        center = box_size / 2.0
        na_atom_idx = mobile_res.atomStart
        cl_atom_idx = mobile_res.atomStart + 1
        
        # Save relative position relationship
        rel_x = system.atoms[cl_atom_idx].x - system.atoms[na_atom_idx].x
        rel_y = system.atoms[cl_atom_idx].y - system.atoms[na_atom_idx].y
        rel_z = system.atoms[cl_atom_idx].z - system.atoms[na_atom_idx].z
        
        # Set Na+ ion to center
        system.atoms[na_atom_idx].x = center
        system.atoms[na_atom_idx].y = center
        system.atoms[na_atom_idx].z = center
        
        # Maintain Cl- relative position
        system.atoms[cl_atom_idx].x = center + rel_x
        system.atoms[cl_atom_idx].y = center + rel_y
        system.atoms[cl_atom_idx].z = center + rel_z
        
        print(f"Moved Na+ to: ({system.atoms[na_atom_idx].x:.4f}, {system.atoms[na_atom_idx].y:.4f}, {system.atoms[na_atom_idx].z:.4f})")
        print(f"Moved Cl- to: ({system.atoms[cl_atom_idx].x:.4f}, {system.atoms[cl_atom_idx].y:.4f}, {system.atoms[cl_atom_idx].z:.4f})")
        
        # Check if position is safe after movement
        current_positions = [
            (system.atoms[na_atom_idx].x, system.atoms[na_atom_idx].y, system.atoms[na_atom_idx].z),
            (system.atoms[cl_atom_idx].x, system.atoms[cl_atom_idx].y, system.atoms[cl_atom_idx].z)
        ]
        is_safe, min_dist = is_safe_position(current_positions, fixed_positions, cutoff, box_size)
        if not is_safe:
            print(f"Warning: Position still not safe after centering, minimum distance is {min_dist} nm")
            print("Proceeding with the test anyway, disabling distance checks")
    else:
        print(f"Initial position safe, minimum distance from fixed particles > {cutoff} nm")


def calculate_complex_movement_vector(system, mobile_res, residues, box_size):
    """Calculate complex directional movement vector based on nearest residue."""
    # Step 2: Calculate moving residue center
    residue_center = [0, 0, 0]
    for i in range(mobile_res.atomCount):
        atom_idx = mobile_res.atomStart + i
        residue_center[0] += system.atoms[atom_idx].x / mobile_res.atomCount
        residue_center[1] += system.atoms[atom_idx].y / mobile_res.atomCount
        residue_center[2] += system.atoms[atom_idx].z / mobile_res.atomCount
    
    # Find nearest fixed residue center
    nearest_distance = float('inf')
    direction_vector = [0, 0, 0]
    
    for i in range(len(residues) - 1):  # Skip last residue (moving residue)
        fixed_res = residues[i]
        fixed_center = [0, 0, 0]
        
        for j in range(fixed_res.atomCount):
            atom_idx = fixed_res.atomStart + j
            fixed_center[0] += system.atoms[atom_idx].x / fixed_res.atomCount
            fixed_center[1] += system.atoms[atom_idx].y / fixed_res.atomCount
            fixed_center[2] += system.atoms[atom_idx].z / fixed_res.atomCount
        
        # Calculate PBC distance
        dist = calculate_pbc_distance(residue_center, fixed_center, box_size)
        if dist < nearest_distance:
            nearest_distance = dist
            
            # Calculate direction vector with PBC
            for k in range(3):
                diff = fixed_center[k] - residue_center[k]
                # Apply PBC to get shortest distance
                if abs(diff) > box_size/2:
                    if diff > 0:
                        diff -= box_size
                    else:
                        diff += box_size
                direction_vector[k] = diff
    
    # To avoid moving toward fixed residues (which would decrease distance), we move in the opposite direction
    direction_vector = [-d for d in direction_vector]
    
    # Normalize direction vector
    magnitude = math.sqrt(sum(d*d for d in direction_vector))
    if magnitude > 0:
        direction_vector = [d/magnitude for d in direction_vector]
    else:
        # Default to diagonal direction if magnitude is zero (unlikely)
        direction_vector = [0.577, 0.577, 0.577]  # 1/sqrt(3) in each dimension
    
    # Set movement distance
    movement_distance = 0.3  # nm - Reduce movement distance to avoid getting too close to box boundaries
    
    # Calculate movement vector
    movement_vector = [d * movement_distance for d in direction_vector]
    print(f"Movement vector: {movement_vector}")
    
    return movement_vector


def apply_movement_to_residue(system, mobile_res, movement_vector, box_size):
    """Apply movement vector to all atoms in the moving residue."""
    # Apply movement to all atoms in the moving residue
    for i in range(mobile_res.atomCount):
        atom_idx = mobile_res.atomStart + i
        system.atoms[atom_idx].x += movement_vector[0]
        system.atoms[atom_idx].y += movement_vector[1]
        system.atoms[atom_idx].z += movement_vector[2]
        
        # Apply PBC
        system.atoms[atom_idx].x %= box_size
        system.atoms[atom_idx].y %= box_size
        system.atoms[atom_idx].z %= box_size
        
        print(f"Moved atom {i} to: ({system.atoms[atom_idx].x:.4f}, {system.atoms[atom_idx].y:.4f}, {system.atoms[atom_idx].z:.4f})")


def verify_moved_position_safety(system, mobile_res, fixed_positions, cutoff, box_size):
    """Verify if moved position is safe and report distances."""
    # Get current positions after movement
    current_positions = []
    for i in range(mobile_res.atomCount):
        atom_idx = mobile_res.atomStart + i
        current_positions.append((
            system.atoms[atom_idx].x,
            system.atoms[atom_idx].y,
            system.atoms[atom_idx].z
        ))
    
    # Verify moved position is safe
    is_safe, min_dist = is_safe_position(current_positions, fixed_positions, cutoff, box_size)
    if not is_safe:
        print(f"Warning: Moved position not safe, minimum distance is {min_dist} nm")
        # Even if distance is smaller than cutoff, we only compare reciprocal space energy changes, so continue testing
        print("Proceeding with the test - we only compare reciprocal space energy changes")
    else:
        min_dist = float('inf')
        for mobile_pos in current_positions:
            for fixed_pos in fixed_positions:
                dist = calculate_pbc_distance(mobile_pos, fixed_pos, box_size)
                min_dist = min(min_dist, dist)
        print(f"Moved position safe, minimum distance from fixed particles: {min_dist:.4f} nm (cutoff={cutoff} nm)")
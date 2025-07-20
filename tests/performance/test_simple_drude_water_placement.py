#!/usr/bin/env python
"""Simple approach: optimize positions with basic model, then replace with Drude water"""

import sys
import numpy as np
import math
sys.path.append('/home/zhaomt/gcmc/test107/pygcmc_dev/build/modules/bindings')

import pygcmc

def create_simple_water_box(n_waters=10, box_size=2.0):
    """Create water box with simple point charges for initial optimization"""
    state = pygcmc.MCState()
    
    # Set box size (nm)
    state.info.box = [box_size, box_size, box_size]
    state.info.cutoff = min(box_size/2 - 0.1, 1.2)  # Reasonable cutoff
    
    # Simple TIP3P-like parameters for initial placement
    # Just oxygen positions, we'll add Drude later
    n_per_dim = int(np.ceil(n_waters**(1/3)))
    spacing = box_size / n_per_dim
    
    atoms = []
    water_id = 0
    for i in range(n_per_dim):
        for j in range(n_per_dim):
            for k in range(n_per_dim):
                if water_id >= n_waters:
                    break
                
                # Place oxygen
                x = (i + 0.5) * spacing
                y = (j + 0.5) * spacing  
                z = (k + 0.5) * spacing
                
                o = pygcmc.MCAtom()
                o.x = x
                o.y = y
                o.z = z
                o.charge = -0.834  # Simple TIP3P charge
                o.type = 0
                atoms.append(o)
                
                water_id += 1
    
    state.atoms = atoms
    state.activeAtomCount = len(atoms)
    return state

def optimize_simple_positions(state, n_steps=100):
    """Simple MC optimization to avoid overlaps"""
    print(f"Optimizing {len(state.atoms)} oxygen positions...")
    
    # Just do some random moves to avoid overlaps
    for step in range(n_steps):
        atom_idx = np.random.randint(0, len(state.atoms))
        
        # Small random displacement
        dx = (np.random.random() - 0.5) * 0.1
        dy = (np.random.random() - 0.5) * 0.1
        dz = (np.random.random() - 0.5) * 0.1
        
        # Move atom
        atom = state.atoms[atom_idx]
        old_x, old_y, old_z = atom.x, atom.y, atom.z
        
        atom.x += dx
        atom.y += dy
        atom.z += dz
        
        # Apply PBC
        box = state.info.box
        atom.x = atom.x % box[0]
        atom.y = atom.y % box[1]
        atom.z = atom.z % box[2]
        
        # Simple check for overlaps (could calculate energy instead)
        overlap = False
        for j in range(len(state.atoms)):
            if j == atom_idx:
                continue
            other = state.atoms[j]
            dx = atom.x - other.x
            dy = atom.y - other.y
            dz = atom.z - other.z
            
            # Minimum image
            if dx > box[0]/2: dx -= box[0]
            if dx < -box[0]/2: dx += box[0]
            if dy > box[1]/2: dy -= box[1]
            if dy < -box[1]/2: dy += box[1]
            if dz > box[2]/2: dz -= box[2]
            if dz < -box[2]/2: dz += box[2]
            
            r2 = dx*dx + dy*dy + dz*dz
            if r2 < 0.09:  # 0.3 nm minimum distance
                overlap = True
                break
        
        if overlap:
            # Reject move
            atom.x = old_x
            atom.y = old_y
            atom.z = old_z
    
    print("Optimization complete")
    return state

def create_swm4_water_at_position(x, y, z):
    """Create SWM4-NDP water atoms at given position"""
    atoms = []
    
    # SWM4-NDP parameters
    qO = 1.71636
    qD = -1.71636
    qH = 0.55733
    qM = -1.11466
    rOH = 0.09572  # nm
    aHOH = 104.52 * math.pi / 180  # radians
    
    # Oxygen
    o = pygcmc.MCAtom()
    o.x, o.y, o.z = x, y, z
    o.charge = qO
    o.type = 0
    atoms.append(o)
    
    # Drude (initially at parent position)
    d = pygcmc.MCAtom()
    d.x, d.y, d.z = x, y, z
    d.charge = qD
    d.type = 1
    atoms.append(d)
    
    # Hydrogen 1
    h1 = pygcmc.MCAtom()
    h1.x = x + rOH
    h1.y = y
    h1.z = z
    h1.charge = qH
    h1.type = 2
    atoms.append(h1)
    
    # Hydrogen 2
    h2 = pygcmc.MCAtom()
    h2.x = x + rOH * math.cos(aHOH)
    h2.y = y + rOH * math.sin(aHOH)
    h2.z = z
    h2.charge = qH
    h2.type = 2
    atoms.append(h2)
    
    # M-site (virtual site)
    m = pygcmc.MCAtom()
    w_O = 0.786646558
    w_H = 0.106676721
    m.x = w_O * o.x + w_H * h1.x + w_H * h2.x
    m.y = w_O * o.y + w_H * h1.y + w_H * h2.y
    m.z = w_O * o.z + w_H * h1.z + w_H * h2.z
    m.charge = qM
    m.type = 3
    atoms.append(m)
    
    return atoms

def convert_to_drude_water(simple_state):
    """Convert simple oxygen positions to full Drude water molecules"""
    drude_state = pygcmc.MCState()
    
    # Copy box info
    drude_state.info.box = simple_state.info.box
    drude_state.info.cutoff = simple_state.info.cutoff
    
    # Initialize Drude force
    pygcmc.initializeDrudeForce()
    
    n_waters = len(simple_state.atoms)
    all_atoms = []
    # For each oxygen, create full SWM4-NDP water
    for i in range(len(simple_state.atoms)):
        o_simple = simple_state.atoms[i]
        
        # Create SWM4-NDP water at this position
        atoms = create_swm4_water_at_position(o_simple.x, o_simple.y, o_simple.z)
        
        # Add all atoms to list
        all_atoms.extend(atoms)
        
        # Add Drude particle (oxygen is parent at index 5*i, drude at 5*i+1)
        pygcmc.addDrudeParticle(
            drudeIndex=5*i + 1,
            parentIndex=5*i,
            charge=-1.71636,
            polarizability=0.000978253  # nm^3
        )
    
    drude_state.atoms = all_atoms
    drude_state.activeAtomCount = len(all_atoms)
    
    # Create residues for each water molecule
    residues = []
    for i in range(n_waters):
        res = pygcmc.MCResidue()
        res.atomStart = 5 * i
        res.atomCount = 5
        res.active = True
        res.type = 0
        residues.append(res)
    
    drude_state.residues = residues
    drude_state.activeResidueCount = len(residues)
    
    # Set up Thole screening between all Drude pairs
    n_waters = len(simple_state.atoms)
    for i in range(n_waters):
        for j in range(i+1, n_waters):
            pygcmc.addDrudeScreenedPair(i, j, 1.3)  # Thole parameter
    
    return drude_state

def main():
    print("=== Simple Drude Water Placement Test ===\n")
    
    # Step 1: Create and optimize simple water positions
    n_waters = 27  # 3x3x3 grid
    box_size = 2.0  # nm
    
    print(f"Creating {n_waters} waters in {box_size} nm box")
    simple_state = create_simple_water_box(n_waters, box_size)
    
    # Optional: optimize positions
    simple_state = optimize_simple_positions(simple_state, n_steps=100)
    
    # Step 2: Convert to Drude water
    print("\nConverting to Drude water molecules...")
    drude_state = convert_to_drude_water(simple_state)
    
    print(f"Created {len(drude_state.atoms)} atoms total")
    print(f"Number of Drude particles: {pygcmc.getNumDrudeParticles()}")
    
    # Step 3: Calculate energy with Drude model
    print("\nCalculating Drude system energy...")
    result = pygcmc.computeSystemEnergyDrude(drude_state)
    
    # Check if result is tuple or single value
    if isinstance(result, tuple):
        energy = result[0]  # Assume first element is total energy
        print(f"Result tuple: {result}")
    else:
        energy = result
    
    print(f"Total energy: {energy:.2f} kJ/mol")
    print(f"Energy per water: {energy/n_waters:.2f} kJ/mol")
    
    # Check if energy is reasonable
    if abs(energy/n_waters) < 1000:
        print("\n✓ Energy looks reasonable!")
    else:
        print("\n✗ Energy seems too high/low, may need parameter adjustment")
    
    # Optional: Save coordinates for visualization
    print("\nFirst few water molecules:")
    for i in range(min(3, n_waters)):
        print(f"\nWater {i+1}:")
        for j in range(5):
            atom = drude_state.atoms[5*i + j]
            atom_type = ["O", "D", "H1", "H2", "M"][j]
            print(f"  {atom_type}: ({atom.x:.3f}, {atom.y:.3f}, {atom.z:.3f}) q={atom.charge:.3f}")

if __name__ == "__main__":
    main()
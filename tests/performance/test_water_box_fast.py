#!/usr/bin/env python
"""Fast performance test for water box - optimized version"""

import sys
import time
import numpy as np
sys.path.append('/home/zhaomt/gcmc/test107/pygcmc_dev/build/modules/bindings')

import pygcmc

def create_water_box_fast(box_size, water_spacing=0.31):
    """Create a box of water molecules - optimized version"""
    state = pygcmc.MCState()
    
    # Set box
    state.info.box = [box_size, box_size, box_size]
    state.info.cutoff = min(box_size/2 - 0.1, 1.2)
    state.info.setTemperature(300.0)
    
    # Calculate number of waters per dimension
    n_per_dim = int(box_size / water_spacing)
    total_waters = n_per_dim ** 3
    
    print(f"Creating water box:")
    print(f"  Box size: {box_size} nm")
    print(f"  Waters per dimension: {n_per_dim}")
    print(f"  Total water molecules: {total_waters}")
    print(f"  Total atoms: {total_waters * 3}")
    print(f"  Cutoff: {state.info.cutoff} nm")
    
    # Pre-allocate all atoms
    atoms = []
    residues = []
    
    # Water geometry
    angle = 104.52 * np.pi / 180
    h1_offset = [0.09572, 0, 0]
    h2_offset = [0.09572 * np.cos(angle), 0.09572 * np.sin(angle), 0]
    
    # Create all waters at once
    mol_id = 0
    for ix in range(n_per_dim):
        for iy in range(n_per_dim):
            for iz in range(n_per_dim):
                x = (ix + 0.5) * water_spacing
                y = (iy + 0.5) * water_spacing
                z = (iz + 0.5) * water_spacing
                
                # Oxygen
                o = pygcmc.MCAtom()
                o.x, o.y, o.z = x, y, z
                o.charge = -0.834
                o.type = 0
                atoms.append(o)
                
                # H1
                h1 = pygcmc.MCAtom()
                h1.x, h1.y, h1.z = x + h1_offset[0], y + h1_offset[1], z + h1_offset[2]
                h1.charge = 0.417
                h1.type = 1
                atoms.append(h1)
                
                # H2
                h2 = pygcmc.MCAtom()
                h2.x, h2.y, h2.z = x + h2_offset[0], y + h2_offset[1], z + h2_offset[2]
                h2.charge = 0.417
                h2.type = 1
                atoms.append(h2)
                
                # Residue
                res = pygcmc.MCResidue()
                res.atomStart = mol_id * 3
                res.atomCount = 3
                res.active = True
                res.fixed = False
                res.type = 0
                residues.append(res)
                
                mol_id += 1
    
    state.atoms = atoms
    state.activeAtomCount = len(atoms)
    state.residues = residues
    state.activeResidueCount = len(residues)
    
    # Set force field
    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 2
    ff.numMovementTypes = 2
    ff.ljSigma = [0.315, 0.0, 0.0, 0.0]
    ff.ljEps = [0.636, 0.0, 0.0, 0.0]
    state.forcefield = ff
    
    return state, total_waters

def main():
    print("=" * 60)
    print("Water Box Performance Test - Optimized Non-Drude")
    print("=" * 60)
    
    # Test different box sizes
    box_sizes = [3.0, 5.0, 10.0]  # nm
    
    for box_size in box_sizes:
        print(f"\n\n=== Testing {box_size} nm box ===")
        
        # Create water box
        state, n_waters = create_water_box_fast(box_size)
        
        # Select center molecule
        target_mol = n_waters // 2
        start_atom = state.residues[target_mol].atomStart
        
        print(f"\nTesting molecule {target_mol}")
        
        # Warmup
        print("\nWarmup (100 energy calculations)...")
        for _ in range(100):
            pygcmc.computeSystemEnergyCutoff(state)
        
        # Time 1000 energy calculations
        print("\nTiming 1000 energy calculations...")
        start_time = time.time()
        
        for _ in range(1000):
            pygcmc.computeSystemEnergyCutoff(state)
            
        end_time = time.time()
        total_time = end_time - start_time
        
        # Get energy
        energy = sum(r.energy_elec + r.energy_vdw for r in state.residues)
        
        print(f"\nResults for {box_size} nm box:")
        print(f"  Total waters: {n_waters}")
        print(f"  Total atoms: {n_waters * 3}")
        print(f"  Energy: {energy:.3f} kJ/mol")
        print(f"  Time for 1000 calculations: {total_time:.3f} s")
        print(f"  Time per calculation: {total_time/1000*1000:.3f} ms")
        print(f"  Calculations per second: {1000/total_time:.1f}")
        
        # Quick movement test
        print(f"\nQuick movement test (1000 steps)...")
        mol_start_time = time.time()
        
        for step in range(1000):
            # Small random displacement
            dx = (np.random.random() - 0.5) * 0.002
            dy = (np.random.random() - 0.5) * 0.002
            dz = (np.random.random() - 0.5) * 0.002
            
            # Move molecule
            for i in range(3):
                state.atoms[start_atom + i].x += dx
                state.atoms[start_atom + i].y += dy
                state.atoms[start_atom + i].z += dz
                
                # Simple PBC
                state.atoms[start_atom + i].x = state.atoms[start_atom + i].x % box_size
                state.atoms[start_atom + i].y = state.atoms[start_atom + i].y % box_size
                state.atoms[start_atom + i].z = state.atoms[start_atom + i].z % box_size
            
            # Calculate energy
            pygcmc.computeSystemEnergyCutoff(state)
        
        mol_end_time = time.time()
        mol_time = mol_end_time - mol_start_time
        
        print(f"  Time for 1000 move+energy steps: {mol_time:.3f} s")
        print(f"  Time per step: {mol_time/1000*1000:.3f} ms")
        print(f"  Steps per second: {1000/mol_time:.1f}")

if __name__ == "__main__":
    main()
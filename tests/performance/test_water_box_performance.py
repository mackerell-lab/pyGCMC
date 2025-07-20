#!/usr/bin/env python
"""Performance test for water box - non-Drude version"""

import sys
import time
import numpy as np
import math
sys.path.append('/home/zhaomt/gcmc/test107/pygcmc_dev/build/modules/bindings')

import pygcmc

def create_water_molecule(x, y, z, mol_id):
    """Create a single water molecule at given position"""
    atoms = []
    
    # Oxygen
    o = pygcmc.MCAtom()
    o.x = x
    o.y = y
    o.z = z
    o.charge = -0.834  # TIP3P charge
    o.type = 0
    atoms.append(o)
    
    # Hydrogen 1
    h1 = pygcmc.MCAtom()
    h1.x = x + 0.09572
    h1.y = y
    h1.z = z
    h1.charge = 0.417
    h1.type = 1
    atoms.append(h1)
    
    # Hydrogen 2
    h2 = pygcmc.MCAtom()
    angle = 104.52 * math.pi / 180
    h2.x = x + 0.09572 * math.cos(angle)
    h2.y = y + 0.09572 * math.sin(angle)
    h2.z = z
    h2.charge = 0.417
    h2.type = 1
    atoms.append(h2)
    
    return atoms

def create_water_box(box_size, water_spacing=0.31):
    """Create a box of water molecules
    
    Args:
        box_size: Box size in nm
        water_spacing: Spacing between water molecules in nm
    """
    state = pygcmc.MCState()
    
    # Set box
    state.info.box = [box_size, box_size, box_size]
    state.info.cutoff = min(box_size/2 - 0.1, 1.2)  # Cutoff at most 1.2 nm
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
    
    # Create water molecules
    all_atoms = []
    residues = []
    mol_id = 0
    
    for ix in range(n_per_dim):
        for iy in range(n_per_dim):
            for iz in range(n_per_dim):
                x = (ix + 0.5) * water_spacing
                y = (iy + 0.5) * water_spacing
                z = (iz + 0.5) * water_spacing
                
                # Add small random displacement to avoid perfect grid
                x += (np.random.random() - 0.5) * 0.02
                y += (np.random.random() - 0.5) * 0.02
                z += (np.random.random() - 0.5) * 0.02
                
                # Create water
                water_atoms = create_water_molecule(x, y, z, mol_id)
                
                # Create residue
                res = pygcmc.MCResidue()
                res.atomStart = len(all_atoms)
                res.atomCount = 3
                res.active = True
                res.fixed = False
                res.type = 0
                residues.append(res)
                
                all_atoms.extend(water_atoms)
                mol_id += 1
    
    state.atoms = all_atoms
    state.activeAtomCount = len(all_atoms)
    state.residues = residues
    state.activeResidueCount = len(residues)
    
    # Set force field - TIP3P water
    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 2  # O and H
    ff.numMovementTypes = 2
    
    # LJ parameters for TIP3P
    # O-O: sigma = 0.315 nm, epsilon = 0.636 kJ/mol
    # H has no LJ
    ljSigma = [0.315, 0.0, 0.0, 0.0]  # O-O, O-H, H-O, H-H
    ljEps = [0.636, 0.0, 0.0, 0.0]
    
    ff.ljSigma = ljSigma
    ff.ljEps = ljEps
    state.forcefield = ff
    
    return state, total_waters

def test_single_molecule_movement(state, mol_index, n_steps=20000):
    """Test moving a single water molecule and calculating energy"""
    
    # Get the residue for this molecule
    res = state.residues[mol_index]
    start_atom = res.atomStart
    
    # Save initial positions
    initial_pos = []
    for i in range(3):
        atom = state.atoms[start_atom + i]
        initial_pos.append([atom.x, atom.y, atom.z])
    
    print(f"\nTesting molecule {mol_index} movement:")
    print(f"  Initial O position: ({initial_pos[0][0]:.3f}, {initial_pos[0][1]:.3f}, {initial_pos[0][2]:.3f})")
    
    # Arrays to store results
    energies = []
    times = []
    
    # Movement parameters
    max_displacement = 0.001  # 0.001 nm = 0.01 Å per step
    max_rotation = 0.01  # radians per step
    
    print(f"\nRunning {n_steps} steps...")
    
    # Initial energy
    start_time = time.time()
    pygcmc.computeSystemEnergyCutoff(state)
    initial_energy = sum(r.energy_elec + r.energy_vdw for r in state.residues)
    end_time = time.time()
    
    print(f"Initial energy: {initial_energy:.3f} kJ/mol")
    print(f"Initial calculation time: {(end_time - start_time)*1000:.3f} ms")
    
    # Progress markers
    progress_steps = [int(n_steps * p) for p in [0.1, 0.25, 0.5, 0.75, 0.9, 1.0]]
    
    # Main loop
    total_start = time.time()
    
    for step in range(n_steps):
        # Random translation
        dx = (np.random.random() - 0.5) * 2 * max_displacement
        dy = (np.random.random() - 0.5) * 2 * max_displacement
        dz = (np.random.random() - 0.5) * 2 * max_displacement
        
        # Random rotation around center of mass
        # For simplicity, just do small rotations around z-axis
        angle = (np.random.random() - 0.5) * 2 * max_rotation
        
        # Get current COM of water
        com_x = sum(state.atoms[start_atom + i].x for i in range(3)) / 3
        com_y = sum(state.atoms[start_atom + i].y for i in range(3)) / 3
        com_z = sum(state.atoms[start_atom + i].z for i in range(3)) / 3
        
        # Move and rotate water
        cos_a = math.cos(angle)
        sin_a = math.sin(angle)
        
        for i in range(3):
            atom = state.atoms[start_atom + i]
            
            # Translate to origin
            rel_x = atom.x - com_x
            rel_y = atom.y - com_y
            
            # Rotate
            new_rel_x = rel_x * cos_a - rel_y * sin_a
            new_rel_y = rel_x * sin_a + rel_y * cos_a
            
            # Translate back and add displacement
            atom.x = new_rel_x + com_x + dx
            atom.y = new_rel_y + com_y + dy
            atom.z += dz
            
            # Apply periodic boundary conditions
            atom.x = atom.x - state.info.box[0] * math.floor(atom.x / state.info.box[0])
            atom.y = atom.y - state.info.box[1] * math.floor(atom.y / state.info.box[1])
            atom.z = atom.z - state.info.box[2] * math.floor(atom.z / state.info.box[2])
        
        # Calculate energy
        step_start = time.time()
        pygcmc.computeSystemEnergyCutoff(state)
        energy = sum(r.energy_elec + r.energy_vdw for r in state.residues)
        step_end = time.time()
        
        energies.append(energy)
        times.append(step_end - step_start)
        
        # Progress update
        if (step + 1) in progress_steps:
            progress = (step + 1) / n_steps * 100
            avg_time = np.mean(times) * 1000  # ms
            print(f"  Progress: {progress:.0f}% ({step+1}/{n_steps} steps)")
            print(f"    Current energy: {energy:.3f} kJ/mol")
            print(f"    Avg time per step: {avg_time:.3f} ms")
    
    total_end = time.time()
    total_time = total_end - total_start
    
    # Analysis
    energies = np.array(energies)
    times = np.array(times)
    
    print(f"\n=== Performance Summary ===")
    print(f"Total simulation time: {total_time:.2f} seconds")
    print(f"Average time per step: {np.mean(times)*1000:.3f} ± {np.std(times)*1000:.3f} ms")
    print(f"Min/Max time per step: {np.min(times)*1000:.3f} / {np.max(times)*1000:.3f} ms")
    print(f"Steps per second: {n_steps/total_time:.1f}")
    
    print(f"\n=== Energy Statistics ===")
    print(f"Average energy: {np.mean(energies):.3f} ± {np.std(energies):.3f} kJ/mol")
    print(f"Min/Max energy: {np.min(energies):.3f} / {np.max(energies):.3f} kJ/mol")
    print(f"Energy drift: {energies[-1] - energies[0]:.3f} kJ/mol")
    
    # Final position
    final_o = state.atoms[start_atom]
    print(f"\nFinal O position: ({final_o.x:.3f}, {final_o.y:.3f}, {final_o.z:.3f})")
    
    return energies, times

def main():
    print("=" * 60)
    print("Water Box Performance Test - Non-Drude Model")
    print("=" * 60)
    
    # Set random seed for reproducibility
    np.random.seed(42)
    
    # Create water box
    box_size = 3.0  # 3 nm box for faster testing, can increase to 10-20 nm
    state, n_waters = create_water_box(box_size)
    
    # Select a water molecule near the center
    center_mol = n_waters // 2
    print(f"\nSelected molecule {center_mol} for movement test")
    
    # Run performance test
    energies, times = test_single_molecule_movement(state, center_mol, n_steps=20000)
    
    # Save some statistics for comparison
    with open('water_box_performance_nondrude.txt', 'w') as f:
        f.write(f"Box size: {box_size} nm\n")
        f.write(f"Total waters: {n_waters}\n")
        f.write(f"Total atoms: {n_waters * 3}\n")
        f.write(f"Steps: 20000\n")
        f.write(f"Avg time per step: {np.mean(times)*1000:.3f} ms\n")
        f.write(f"Total time: {np.sum(times):.2f} s\n")
        f.write(f"Avg energy: {np.mean(energies):.3f} kJ/mol\n")
    
    print(f"\nResults saved to water_box_performance_nondrude.txt")

if __name__ == "__main__":
    main()
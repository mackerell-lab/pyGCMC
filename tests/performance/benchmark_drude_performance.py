#!/usr/bin/env python
"""Benchmark Drude implementation performance"""

import sys
import time
import numpy as np
import math
sys.path.append('/home/zhaomt/gcmc/test107/pygcmc_dev/build/modules/bindings')

import pygcmc

def create_water_box(n_waters, density=1000):  # kg/m^3
    """Create a box of SWM4-NDP water molecules"""
    # Calculate box size from density
    # Mass of one water = 18.01528 g/mol
    # For n waters: mass = n * 18.01528 / 6.022e23 g = n * 2.992e-23 g = n * 2.992e-26 kg
    # Volume = mass / density
    # For water density ~1000 kg/m^3
    mass_per_water = 18.01528 / 6.022e23 * 1e-3  # kg
    total_mass = n_waters * mass_per_water  # kg
    volume = total_mass / density  # m^3
    box_size = (volume * 1e27) ** (1/3)  # nm
    
    print(f"Creating {n_waters} waters in {box_size:.3f} nm box (density={density} kg/m³)")
    
    state = pygcmc.MCState()
    state.info.box = [box_size, box_size, box_size]
    state.info.cutoff = min(box_size/2 - 0.1, 1.2)
    
    # Initialize Drude force
    pygcmc.initializeDrudeForce()
    
    # Use OpenMM default SCF parameters
    scf_params = pygcmc.DrudeSCFParams()
    # OpenMM defaults: tolerance=1.0, maxIterations=50, dampingFactor=0.5
    # We keep the defaults, just set maxDrudeDistance for safety
    scf_params.maxDrudeDistance = 0.02  # OpenMM default is 0.0, but 0.02 is safer
    pygcmc.setDrudeSCFParameters(scf_params)
    
    atoms = []
    residues = []
    
    # SWM4-NDP parameters
    qO = 1.71636
    qD = -1.71636
    qH = 0.55733
    qM = -1.11466
    rOH = 0.09572  # nm
    aHOH = 104.52 * math.pi / 180
    
    # Place waters on a grid
    n_per_dim = int(np.ceil(n_waters**(1/3)))
    spacing = box_size / n_per_dim
    
    water_id = 0
    for i in range(n_per_dim):
        for j in range(n_per_dim):
            for k in range(n_per_dim):
                if water_id >= n_waters:
                    break
                
                # Position with small random displacement to avoid perfect grid
                x = (i + 0.5) * spacing + (np.random.random() - 0.5) * 0.05
                y = (j + 0.5) * spacing + (np.random.random() - 0.5) * 0.05
                z = (k + 0.5) * spacing + (np.random.random() - 0.5) * 0.05
                
                # Random orientation
                theta = np.random.random() * 2 * math.pi
                phi = np.random.random() * math.pi
                psi = np.random.random() * 2 * math.pi
                
                # Rotation matrix
                ct, st = math.cos(theta), math.sin(theta)
                cp, sp = math.cos(phi), math.sin(phi)
                cps, sps = math.cos(psi), math.sin(psi)
                
                # Water geometry in local frame
                # O at origin, H1 along x, H2 in xy plane
                h1_local = np.array([rOH, 0, 0])
                h2_local = np.array([rOH * math.cos(aHOH), rOH * math.sin(aHOH), 0])
                
                # Rotate to global frame
                rot = np.array([
                    [ct*cp, ct*sp*sps - st*cps, ct*sp*cps + st*sps],
                    [st*cp, st*sp*sps + ct*cps, st*sp*cps - ct*sps],
                    [-sp, cp*sps, cp*cps]
                ])
                
                h1_global = rot @ h1_local
                h2_global = rot @ h2_local
                
                # Create atoms
                # Oxygen
                o = pygcmc.MCAtom()
                o.x, o.y, o.z = x, y, z
                o.charge = qO
                o.type = 0
                atoms.append(o)
                
                # Drude
                d = pygcmc.MCAtom()
                d.x, d.y, d.z = x, y, z
                d.charge = qD
                d.type = 1
                atoms.append(d)
                
                # Hydrogens
                h1 = pygcmc.MCAtom()
                h1.x = x + h1_global[0]
                h1.y = y + h1_global[1]
                h1.z = z + h1_global[2]
                h1.charge = qH
                h1.type = 2
                atoms.append(h1)
                
                h2 = pygcmc.MCAtom()
                h2.x = x + h2_global[0]
                h2.y = y + h2_global[1]
                h2.z = z + h2_global[2]
                h2.charge = qH
                h2.type = 2
                atoms.append(h2)
                
                # M-site
                m = pygcmc.MCAtom()
                w_O = 0.786646558
                w_H = 0.106676721
                m.x = w_O * o.x + w_H * h1.x + w_H * h2.x
                m.y = w_O * o.y + w_H * h1.y + w_H * h2.y
                m.z = w_O * o.z + w_H * h1.z + w_H * h2.z
                m.charge = qM
                m.type = 3
                atoms.append(m)
                
                # Create residue
                res = pygcmc.MCResidue()
                res.atomStart = 5 * water_id
                res.atomCount = 5
                res.active = True
                res.type = 0
                residues.append(res)
                
                # Add Drude particle
                pygcmc.addDrudeParticle(
                    drudeIndex=5*water_id + 1,
                    parentIndex=5*water_id,
                    charge=qD,
                    polarizability=0.000978253
                )
                
                water_id += 1
    
    state.atoms = atoms
    state.activeAtomCount = len(atoms)
    state.residues = residues
    state.activeResidueCount = len(residues)
    
    # Add Thole screening
    print(f"Adding Thole screening for {n_waters} waters...")
    for i in range(n_waters):
        for j in range(i+1, n_waters):
            pygcmc.addDrudeScreenedPair(i, j, 1.3)
    
    return state

def benchmark_energy_calculation(state, n_waters, n_iterations=10):
    """Benchmark energy calculation speed"""
    print(f"\nBenchmarking {n_iterations} energy calculations...")
    
    times = []
    energies = []
    
    # Warm up
    pygcmc.computeSystemEnergyDrude(state)
    
    for i in range(n_iterations):
        start = time.time()
        result = pygcmc.computeSystemEnergyDrude(state)
        elapsed = time.time() - start
        times.append(elapsed)
        
        if isinstance(result, tuple):
            energies.append(result[0])
        else:
            energies.append(result)
    
    times = np.array(times)
    energies = np.array(energies)
    
    print(f"  Average time: {times.mean()*1000:.2f} ms")
    print(f"  Std dev: {times.std()*1000:.2f} ms")
    print(f"  Min time: {times.min()*1000:.2f} ms")
    print(f"  Max time: {times.max()*1000:.2f} ms")
    print(f"  Average energy: {energies.mean():.2f} kJ/mol ({energies.mean()/n_waters:.2f} per water)")
    print(f"  Energy std dev: {energies.std():.2f} kJ/mol")
    
    return times.mean()

def benchmark_scaling():
    """Test scaling with system size"""
    print("=== Drude Performance Benchmark ===\n")
    
    system_sizes = [8, 27, 64]  # Perfect cubes for easier setup
    results = []
    
    for n_waters in system_sizes:
        print(f"\n--- System: {n_waters} waters ---")
        
        # Create system
        start = time.time()
        state = create_water_box(n_waters)
        setup_time = time.time() - start
        print(f"Setup time: {setup_time:.2f} s")
        
        # Benchmark
        avg_time = benchmark_energy_calculation(state, n_waters, n_iterations=5)
        
        results.append({
            'n_waters': n_waters,
            'n_atoms': n_waters * 5,
            'setup_time': setup_time,
            'energy_time': avg_time,
            'time_per_water': avg_time / n_waters * 1000  # ms
        })
    
    # Summary
    print("\n\n=== SCALING SUMMARY ===")
    print("\n| Waters | Atoms | Setup(s) | Energy(ms) | ms/water |")
    print("|--------|-------|----------|------------|----------|")
    for r in results:
        print(f"| {r['n_waters']:6d} | {r['n_atoms']:5d} | {r['setup_time']:8.2f} | "
              f"{r['energy_time']*1000:10.2f} | {r['time_per_water']:8.3f} |")
    
    # Calculate scaling
    if len(results) > 1:
        # Log-log fit to determine scaling: time = a * N^b
        x = np.log([r['n_waters'] for r in results])
        y = np.log([r['energy_time'] for r in results])
        
        # Linear regression in log space
        A = np.vstack([x, np.ones(len(x))]).T
        b, log_a = np.linalg.lstsq(A, y, rcond=None)[0]
        
        print(f"\nScaling: time ∝ N^{b:.2f}")
        if b < 1.5:
            print("Good scaling! (better than O(N²))")
        elif b < 2.5:
            print("Acceptable scaling (approximately O(N²))")
        else:
            print("Poor scaling (worse than O(N²))")

def compare_with_reference():
    """Compare with reference values"""
    print("\n\n=== Reference Comparison ===\n")
    
    # Create standard test system
    n_waters = 64
    state = create_water_box(n_waters, density=997)  # Standard water density at 25°C
    
    # Calculate energy
    result = pygcmc.computeSystemEnergyDrude(state)
    if isinstance(result, tuple):
        energy = result[0]
        components = result[1]
    else:
        energy = result
        components = {}
    
    energy_per_water = energy / n_waters
    
    print(f"PyGCMC Drude Results:")
    print(f"  Total energy: {energy:.2f} kJ/mol")
    print(f"  Energy per water: {energy_per_water:.2f} kJ/mol")
    print(f"  Components: {components}")
    
    print(f"\nExpected ranges (from literature):")
    print(f"  Liquid water at 298K: -41 to -46 kJ/mol per molecule")
    print(f"  Gas phase water dimer: -20 to -25 kJ/mol")
    
    if -50 < energy_per_water < -35:
        print(f"\n✓ Energy is in reasonable range for liquid water")
    else:
        print(f"\n✗ Energy seems off - may need parameter adjustment")

if __name__ == "__main__":
    benchmark_scaling()
    compare_with_reference()
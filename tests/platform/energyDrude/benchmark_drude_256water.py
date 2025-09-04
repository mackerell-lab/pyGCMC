#!/usr/bin/env python3
"""
Benchmark comparison of three Drude algorithms on 256 water molecules
Tests SCF, OPT3, and FBP algorithms
"""

import time
import numpy as np
import pygcmc

def create_simple_water_system(n_waters=256):
    """Create a simple water system for benchmarking"""
    state = pygcmc.MCState()
    
    # Box size calculation
    # Density ~ 1 g/cm³, each water ~18 g/mol
    # Volume = n_waters * 18 / (6.022e23 * 1.0) cm³
    volume = n_waters * 18.0 / (6.022e23 * 1.0) * 1e21  # in nm³
    box_size = volume ** (1.0/3.0)
    box_size *= 1.1  # Add 10% for comfort
    
    state.info.box = [box_size, box_size, box_size]
    # ang might not be needed or have different name
    
    # Setup minimal forcefield
    state.forcefield.numTotalTypes = 4  # O, D, H, M
    state.forcefield.ljSigma = [0.0, 0.0, 0.0, 0.0]
    state.forcefield.ljEps = [0.0, 0.0, 0.0, 0.0]
    
    # Simple cubic grid layout
    n_per_side = int(np.ceil(n_waters ** (1.0/3.0)))
    spacing = box_size / n_per_side
    
    water_id = 0
    for ix in range(n_per_side):
        for iy in range(n_per_side):
            for iz in range(n_per_side):
                if water_id >= n_waters:
                    break
                    
                # Base position
                x = (ix + 0.5) * spacing
                y = (iy + 0.5) * spacing
                z = (iz + 0.5) * spacing
                
                # Create simplified water (just O and D for Drude test)
                # Oxygen
                atom_O = pygcmc.MCAtom()
                atom_O.x, atom_O.y, atom_O.z = x, y, z
                atom_O.type = 0
                atom_O.charge = 1.71636  # SWM4-NDP
                
                # Drude
                atom_D = pygcmc.MCAtom()
                atom_D.x = x + 0.001  # Slight displacement
                atom_D.y = y
                atom_D.z = z
                atom_D.type = 1
                atom_D.charge = -1.71636
                # Mark as Drude through type system instead
                
                # Simplified H atoms (for charge balance)
                atom_H1 = pygcmc.MCAtom()
                atom_H1.x = x + 0.1
                atom_H1.y = y
                atom_H1.z = z
                atom_H1.type = 2
                atom_H1.charge = 0.55733
                
                atom_H2 = pygcmc.MCAtom()
                atom_H2.x = x
                atom_H2.y = y + 0.1
                atom_H2.z = z
                atom_H2.type = 2
                atom_H2.charge = 0.55733
                
                # M-site
                atom_M = pygcmc.MCAtom()
                atom_M.x = x
                atom_M.y = y
                atom_M.z = z + 0.05
                atom_M.type = 3
                atom_M.charge = -1.11466
                
                # Add to state
                state.atoms.append(atom_O)
                state.atoms.append(atom_D)
                state.atoms.append(atom_H1)
                state.atoms.append(atom_H2)
                state.atoms.append(atom_M)
                
                # Create residue
                residue = pygcmc.MCResidue()
                residue.name = "WAT"
                residue.atomStart = water_id * 5
                residue.atomCount = 5
                state.residues.append(residue)
                
                water_id += 1
    
    state.activeAtomCount = n_waters * 5
    state.activeResidueCount = n_waters
    
    return state

def setup_drude_complete():
    """Setup Drude system using DrudeComplete"""
    # Clear any existing particles
    pygcmc.DrudeComplete.clear()
    
    # Add Drude particle for each water oxygen
    # In our simplified system, every 5th atom starting from index 1 is a Drude
    particle_id = 0
    for water_id in range(256):
        particle = pygcmc.DrudeParticle()
        particle.drudeIndex = water_id * 5 + 1  # Drude is second atom in each water
        particle.parentIndex = water_id * 5     # Oxygen is first atom
        particle.aniso1Index = -1
        particle.aniso2Index = -1
        particle.aniso3Index = -1
        particle.aniso4Index = -1
        particle.charge = -1.71636
        particle.polarizability = 0.00097822  # nm³
        particle.aniso12 = 1.0
        particle.aniso34 = 1.0
        particle.computeSpringConstants()
        
        pygcmc.DrudeComplete.addParticle(particle)
        particle_id += 1
    
    return particle_id

def benchmark_algorithm(state, algorithm_type, n_iterations=50):
    """Benchmark a specific algorithm"""
    
    # Setup Drude particles
    n_particles = setup_drude_complete()
    
    # Set algorithm
    if algorithm_type == "SCF":
        pygcmc.DrudeComplete.setAlgorithm(pygcmc.DrudeAlgorithm.SCF)
        params = pygcmc.DrudeSCFParams()
        params.tolerance = 1.0
        params.maxIterations = 50
        params.dampingFactor = 0.5
        params.maxDrudeDistance = 0.02
        pygcmc.DrudeComplete.setParameters(params)
    elif algorithm_type == "OPT3":
        pygcmc.DrudeComplete.setAlgorithm(pygcmc.DrudeAlgorithm.OPT3)
    elif algorithm_type == "FBP":
        pygcmc.DrudeComplete.setAlgorithm(pygcmc.DrudeAlgorithm.FBP)
    
    # Warmup
    print(f"\n{algorithm_type}: Warming up...")
    for _ in range(3):
        energy = pygcmc.DrudeComplete.calculateEnergy(state)
    
    # Benchmark
    print(f"{algorithm_type}: Running {n_iterations} iterations...")
    energies = []
    
    start_time = time.time()
    
    for i in range(n_iterations):
        # Small perturbation every 10 iterations
        if i > 0 and i % 10 == 0:
            idx = np.random.randint(0, state.activeAtomCount)
            # Skip Drude atoms (type 1) and only perturb real atoms
            if state.atoms[idx].type != 1:
                state.atoms[idx].x += 0.0001 * (np.random.rand() - 0.5)
        
        energy = pygcmc.DrudeComplete.calculateEnergy(state)
        energies.append(energy)
    
    end_time = time.time()
    
    # Calculate stats
    total_time = end_time - start_time
    time_per_iter = total_time / n_iterations
    
    return {
        'algorithm': algorithm_type,
        'total_time': total_time,
        'time_per_iter_ms': time_per_iter * 1000,
        'iterations_per_sec': 1.0 / time_per_iter if time_per_iter > 0 else 0,
        'avg_energy': np.mean(energies),
        'std_energy': np.std(energies)
    }

def main():
    """Run the benchmark comparison"""
    print("="*60)
    print("Drude Algorithm Performance Comparison")
    print("System: 256 water molecules")
    print("="*60)
    
    # Create system
    print("\nCreating water system...")
    state = create_simple_water_system(256)
    print(f"Created: {state.activeAtomCount} atoms ({state.activeResidueCount} waters)")
    print(f"Box: {state.info.box[0]:.2f} nm³")
    
    # Test algorithms
    algorithms = ["SCF", "OPT3", "FBP"]
    results = []
    
    for algo in algorithms:
        try:
            result = benchmark_algorithm(state.copy(), algo)
            results.append(result)
        except Exception as e:
            print(f"Error with {algo}: {e}")
            # Try with reduced parameters
            if algo == "OPT3" or algo == "FBP":
                print(f"Note: {algo} may not be fully implemented yet")
    
    # Display results
    print("\n" + "="*60)
    print("RESULTS")
    print("="*60)
    print(f"{'Algorithm':<10} {'Time/iter (ms)':<15} {'Iter/sec':<12} {'Energy (kJ/mol)':<15}")
    print("-"*60)
    
    if results:
        # Find fastest
        min_time = min(r['time_per_iter_ms'] for r in results)
        
        for r in results:
            speedup = min_time / r['time_per_iter_ms'] if r['time_per_iter_ms'] > 0 else 0
            print(f"{r['algorithm']:<10} "
                  f"{r['time_per_iter_ms']:>13.2f} "
                  f"{r['iterations_per_sec']:>11.1f} "
                  f"{r['avg_energy']:>14.2f} "
                  f"(x{speedup:.2f})")
    
    # Performance comparison
    if len(results) > 1:
        print("\n" + "="*60)
        print("RELATIVE PERFORMANCE")
        print("="*60)
        
        scf_result = next((r for r in results if r['algorithm'] == 'SCF'), None)
        if scf_result and scf_result['time_per_iter_ms'] > 0:
            for r in results:
                if r['algorithm'] != 'SCF':
                    speedup = scf_result['time_per_iter_ms'] / r['time_per_iter_ms']
                    print(f"{r['algorithm']} is {speedup:.2f}x faster than SCF")

if __name__ == "__main__":
    main()
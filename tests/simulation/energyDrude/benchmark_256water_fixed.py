#!/usr/bin/env python3
"""
Benchmark comparison of three Drude algorithms on 256 water molecules
"""

import time
import numpy as np
import pygcmc
from helpers import create_water_molecule, setup_water_system

def benchmark_algorithm(algorithm_type, n_iterations=50):
    """Benchmark a specific algorithm"""
    
    # Create 256 water system
    print(f"\n{algorithm_type}: Creating 256 water system...")
    state, positions = setup_water_system(256, box_size=2.5)
    
    # Setup Drude system
    pygcmc.DrudeComplete.clear()
    
    # Add Drude particles for each water
    n_waters = 256
    for i in range(n_waters):
        particle = pygcmc.DrudeParticle()
        particle.drudeIndex = i * 5 + 1  # Drude is second atom
        particle.parentIndex = i * 5      # Oxygen is first
        particle.charge = -1.71636       # SWM4-NDP charge
        particle.polarizability = 0.00097822  # nm³
        particle.aniso1Index = -1
        particle.aniso2Index = -1
        particle.aniso3Index = -1
        particle.aniso4Index = -1
        particle.aniso12 = 1.0
        particle.aniso34 = 1.0
        particle.computeSpringConstants()
        pygcmc.DrudeComplete.addParticle(particle)
    
    # Add some Thole screening pairs (subset for performance testing)
    # In real simulation, would add all pairs
    for i in range(min(50, n_waters)):
        for j in range(i+1, min(i+3, n_waters)):
            pair = pygcmc.ScreenedPair()
            pair.dipole1 = i
            pair.dipole2 = j
            pair.thole = 1.3
            pygcmc.DrudeComplete.addScreenedPair(pair)
    
    # Set SCF parameters (used by all algorithms)
    params = pygcmc.DrudeSCFParams()
    params.tolerance = 1.0
    params.maxIterations = 50
    params.dampingFactor = 0.5
    params.maxDrudeDistance = 0.02
    pygcmc.DrudeComplete.setParameters(params)
    
    # Get algorithm enum
    if algorithm_type == "SCF":
        algo = pygcmc.DrudeAlgorithm.SCF
    elif algorithm_type == "OPT3":
        algo = pygcmc.DrudeAlgorithm.OPT3
    elif algorithm_type == "FBP":
        algo = pygcmc.DrudeAlgorithm.FBP
    else:
        raise ValueError(f"Unknown algorithm: {algorithm_type}")
    
    # Warmup
    print(f"{algorithm_type}: Warming up...")
    for _ in range(3):
        if algorithm_type == "SCF":
            # Use default SCF
            energy = pygcmc.DrudeComplete.calculateEnergy(state)
        else:
            # Use specific algorithm
            energy = pygcmc.DrudeComplete.calculateEnergy(state, algo)
    
    # Benchmark
    print(f"{algorithm_type}: Running {n_iterations} iterations...")
    energies = []
    
    start_time = time.time()
    
    for i in range(n_iterations):
        # Small perturbation every 10 iterations
        if i > 0 and i % 10 == 0:
            idx = np.random.randint(0, state.activeAtomCount)
            # Skip Drude atoms (type 1)
            if state.atoms[idx].type != 1:
                state.atoms[idx].x += 0.0001 * (np.random.rand() - 0.5)
        
        if algorithm_type == "SCF":
            energy = pygcmc.DrudeComplete.calculateEnergy(state)
        else:
            energy = pygcmc.DrudeComplete.calculateEnergy(state, algo)
        energies.append(energy)
    
    end_time = time.time()
    
    # Calculate stats
    total_time = end_time - start_time
    time_per_iter = total_time / n_iterations
    
    # Clean up
    pygcmc.DrudeComplete.clear()
    
    return {
        'algorithm': algorithm_type,
        'total_time': total_time,
        'time_per_iter_ms': time_per_iter * 1000,
        'iterations_per_sec': 1.0 / time_per_iter if time_per_iter > 0 else 0,
        'avg_energy': np.mean(energies),
        'std_energy': np.std(energies),
        'min_energy': np.min(energies),
        'max_energy': np.max(energies)
    }

def main():
    """Run the benchmark comparison"""
    print("="*60)
    print("Drude Algorithm Performance Comparison")
    print("System: 256 water molecules (SWM4-NDP)")
    print("="*60)
    
    # Test algorithms
    algorithms = ["SCF", "OPT3", "FBP"]
    results = []
    
    for algo in algorithms:
        try:
            result = benchmark_algorithm(algo, n_iterations=50)
            results.append(result)
            print(f"{algo} completed successfully")
        except Exception as e:
            print(f"Error with {algo}: {e}")
            import traceback
            traceback.print_exc()
    
    # Display results
    if results:
        print("\n" + "="*80)
        print("RESULTS")
        print("="*80)
        print(f"{'Algorithm':<10} {'Time/iter (ms)':<15} {'Iter/sec':<12} {'Avg E (kJ/mol)':<15} {'Std E':<10}")
        print("-"*80)
        
        # Find fastest
        min_time = min(r['time_per_iter_ms'] for r in results)
        
        for r in results:
            speedup = min_time / r['time_per_iter_ms'] if r['time_per_iter_ms'] > 0 else 0
            print(f"{r['algorithm']:<10} "
                  f"{r['time_per_iter_ms']:>13.2f} "
                  f"{r['iterations_per_sec']:>11.1f} "
                  f"{r['avg_energy']:>14.2f} "
                  f"{r['std_energy']:>9.4f} "
                  f"(x{speedup:.2f})")
        
        # Performance comparison
        if len(results) > 1:
            print("\n" + "="*80)
            print("RELATIVE PERFORMANCE")
            print("="*80)
            
            scf_result = next((r for r in results if r['algorithm'] == 'SCF'), None)
            if scf_result and scf_result['time_per_iter_ms'] > 0:
                for r in results:
                    if r['algorithm'] != 'SCF':
                        speedup = scf_result['time_per_iter_ms'] / r['time_per_iter_ms']
                        if speedup > 1:
                            print(f"{r['algorithm']} is {speedup:.2f}x faster than SCF")
                        else:
                            print(f"{r['algorithm']} is {1/speedup:.2f}x slower than SCF")
            
            # Check energy consistency
            print("\n" + "="*80)
            print("ENERGY CONSISTENCY CHECK")
            print("="*80)
            energies = [r['avg_energy'] for r in results]
            energy_std = np.std(energies)
            print(f"Energy standard deviation across algorithms: {energy_std:.6f} kJ/mol")
            if energy_std > 1.0:
                print("WARNING: Large energy differences between algorithms!")
                for r in results:
                    print(f"  {r['algorithm']}: {r['avg_energy']:.2f} kJ/mol")
            else:
                print("Good: All algorithms converge to similar energies")

if __name__ == "__main__":
    main()
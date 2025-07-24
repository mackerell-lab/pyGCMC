#!/usr/bin/env python3
"""
Benchmark with full Thole screening to check energy consistency
"""

import time
import numpy as np
import pygcmc
from helpers import setup_water_system

def benchmark_with_full_thole(algorithm_type, n_waters=27):
    """Benchmark with complete Thole screening"""
    
    # Create smaller system for full Thole test
    print(f"\n{algorithm_type}: Creating {n_waters} water system...")
    state, positions = setup_water_system(n_waters, box_size=1.5)
    
    # Setup Drude system
    pygcmc.DrudeComplete.clear()
    
    # Add Drude particles
    for i in range(n_waters):
        particle = pygcmc.DrudeParticle()
        particle.drudeIndex = i * 5 + 1
        particle.parentIndex = i * 5
        particle.charge = -1.71636
        particle.polarizability = 0.00097822
        particle.aniso1Index = -1
        particle.aniso2Index = -1
        particle.aniso3Index = -1
        particle.aniso4Index = -1
        particle.aniso12 = 1.0
        particle.aniso34 = 1.0
        particle.computeSpringConstants()
        pygcmc.DrudeComplete.addParticle(particle)
    
    # Add ALL Thole screening pairs
    print(f"Adding {n_waters * (n_waters - 1) // 2} Thole pairs...")
    for i in range(n_waters):
        for j in range(i+1, n_waters):
            pair = pygcmc.ScreenedPair()
            pair.dipole1 = i
            pair.dipole2 = j
            pair.thole = 1.3
            pygcmc.DrudeComplete.addScreenedPair(pair)
    
    # Set parameters
    params = pygcmc.DrudeSCFParams()
    params.tolerance = 0.1  # Tighter tolerance
    params.maxIterations = 100
    params.dampingFactor = 0.5
    params.maxDrudeDistance = 0.02
    pygcmc.DrudeComplete.setParameters(params)
    
    # Get algorithm
    if algorithm_type == "SCF":
        algo = pygcmc.DrudeAlgorithm.SCF
    elif algorithm_type == "OPT3":
        algo = pygcmc.DrudeAlgorithm.OPT3
    elif algorithm_type == "FBP":
        algo = pygcmc.DrudeAlgorithm.FBP
    
    # Calculate energy
    if algorithm_type == "SCF":
        energy = pygcmc.DrudeComplete.calculateEnergy(state)
    else:
        energy = pygcmc.DrudeComplete.calculateEnergy(state, algo)
    
    # Get Drude positions for analysis
    drude_displacements = []
    for i in range(n_waters):
        parent_idx = i * 5
        drude_idx = i * 5 + 1
        dx = state.atoms[drude_idx].x - state.atoms[parent_idx].x
        dy = state.atoms[drude_idx].y - state.atoms[parent_idx].y
        dz = state.atoms[drude_idx].z - state.atoms[parent_idx].z
        dist = np.sqrt(dx*dx + dy*dy + dz*dz)
        drude_displacements.append(dist)
    
    pygcmc.DrudeComplete.clear()
    
    return {
        'algorithm': algorithm_type,
        'energy': energy,
        'energy_per_water': energy / n_waters,
        'avg_drude_displacement': np.mean(drude_displacements),
        'max_drude_displacement': np.max(drude_displacements),
        'n_waters': n_waters
    }

def main():
    """Compare energies with full Thole screening"""
    print("="*60)
    print("Drude Algorithm Energy Comparison")
    print("Full Thole Screening Test")
    print("="*60)
    
    # Test on smaller system with full Thole
    algorithms = ["SCF", "OPT3", "FBP"]
    results = []
    
    for algo in algorithms:
        try:
            result = benchmark_with_full_thole(algo, n_waters=27)
            results.append(result)
            print(f"{algo} completed: E = {result['energy']:.2f} kJ/mol")
        except Exception as e:
            print(f"Error with {algo}: {e}")
    
    # Display comparison
    if results:
        print("\n" + "="*60)
        print("ENERGY COMPARISON (27 waters, full Thole screening)")
        print("="*60)
        print(f"{'Algorithm':<10} {'Total E (kJ/mol)':<18} {'E/water':<15} {'Avg d_Drude (nm)':<18}")
        print("-"*60)
        
        for r in results:
            print(f"{r['algorithm']:<10} "
                  f"{r['energy']:>16.2f} "
                  f"{r['energy_per_water']:>14.2f} "
                  f"{r['avg_drude_displacement']:>17.4f}")
        
        # Check consistency
        energies = [r['energy'] for r in results]
        energy_range = max(energies) - min(energies)
        print(f"\nEnergy range: {energy_range:.2f} kJ/mol")
        if energy_range > 1.0:
            print("WARNING: Algorithms giving different energies!")
            print("\nPossible causes:")
            print("1. OPT3/FBP implementations may be incomplete")
            print("2. Different convergence criteria between algorithms")
            print("3. Initial Drude positions affecting results")

if __name__ == "__main__":
    main()
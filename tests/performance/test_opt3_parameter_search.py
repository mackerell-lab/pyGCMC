#!/usr/bin/env python
"""Comprehensive search for optimal OPT3 parameters"""

import sys
import numpy as np
import time
from itertools import product
sys.path.append('/home/zhaomt/gcmc/test107/pygcmc_dev/build/modules/bindings')

import pygcmc

def create_stable_water_system(n_waters):
    """Create a stable water system for testing"""
    state = pygcmc.MCState()
    
    # Use larger box to avoid convergence issues
    box_size = (n_waters * 0.035)**(1/3)  # Slightly lower density
    state.info.box = [box_size, box_size, box_size]
    state.info.cutoff = min(0.9, box_size/2 - 0.1)
    
    atoms = []
    charges = [1.71636, -1.71636, 0.55733, 0.55733, -1.11466]
    types = [0, 1, 2, 2, 3]
    
    # Grid placement
    n_per_side = int(np.ceil(n_waters**(1/3)))
    spacing = box_size / n_per_side
    
    n_placed = 0
    for ix in range(n_per_side):
        for iy in range(n_per_side):
            for iz in range(n_per_side):
                if n_placed >= n_waters:
                    break
                
                x = (ix + 0.5) * spacing
                y = (iy + 0.5) * spacing  
                z = (iz + 0.5) * spacing
                
                # Water positions
                positions = [
                    [x, y, z],                    # O
                    [x, y, z],                    # D
                    [x + 0.09572, y, z],         # H1
                    [x - 0.03, y + 0.09, z],     # H2
                    [x + 0.015, y + 0.011, z]    # M-site
                ]
                
                for j in range(5):
                    a = pygcmc.MCAtom()
                    a.x, a.y, a.z = positions[j]
                    a.charge = charges[j]
                    a.type = types[j]
                    atoms.append(a)
                
                n_placed += 1
    
    state.atoms = atoms
    state.activeAtomCount = len(atoms)
    
    # Create residues
    residues = []
    for i in range(n_waters):
        res = pygcmc.MCResidue()
        res.atomStart = 5 * i
        res.atomCount = 5
        res.active = True
        res.type = 0
        residues.append(res)
    
    state.residues = residues
    state.activeResidueCount = len(residues)
    
    # Force field
    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 4
    ff.numMovementTypes = 4
    ff.ljSigma = [0.318395] + [0.0] * 15
    ff.ljEps = [0.88257] + [0.0] * 15
    state.forcefield = ff
    
    return state

def grid_search_parameters():
    """Grid search for optimal OPT3 parameters"""
    print("=== Grid Search for Optimal OPT3 Parameters ===\n")
    
    # Create test systems
    test_systems = [
        create_stable_water_system(8),
        create_stable_water_system(16),
    ]
    
    # Generate coefficient combinations
    # c0 = 0 (fixed), c1 + c2 + c3 = 1
    steps = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    
    valid_coeffs = []
    for c1 in steps:
        for c2 in steps:
            c3 = 1.0 - c1 - c2
            if 0 <= c3 <= 1.0 and abs(c1 + c2 + c3 - 1.0) < 1e-6:
                valid_coeffs.append([0.0, c1, c2, c3])
    
    print(f"Testing {len(valid_coeffs)} coefficient combinations\n")
    
    # Test each combination
    results = []
    
    for coeffs in valid_coeffs:
        total_error = 0
        total_accuracy = 0
        n_tests = 0
        
        for i, state in enumerate(test_systems):
            n_waters = [8, 16][i]
            
            # Setup Drude force
            drude_force = pygcmc.DrudeForce()
            
            for j in range(n_waters):
                drude_force.addParticle(
                    drudeIndex=5*j + 1,
                    parentIndex=5*j,
                    aniso1Index=-1, aniso2Index=-1,
                    aniso3Index=-1, aniso4Index=-1,
                    charge=-1.71636,
                    polarizability=0.000978253,
                    aniso12=0.0, aniso34=0.0
                )
            
            # Add screened pairs
            for j in range(n_waters):
                for k in range(j+1, n_waters):
                    drude_force.addScreenedPair(j, k, 1.3)
            
            # Get SCF energy
            drude_force.setAlgorithm(pygcmc.DrudeAlgorithm.SCF)
            energy_scf = drude_force.calculateEnergySCF(state)
            
            # Get OPT3 energy
            drude_force.setAlgorithm(pygcmc.DrudeAlgorithm.OPT3)
            drude_force.setOPT3Coefficients(*coeffs)
            
            # Reset Drude positions
            for j in range(n_waters):
                drude_idx = 5*j + 1
                parent_idx = 5*j
                state.atoms[drude_idx].x = state.atoms[parent_idx].x
                state.atoms[drude_idx].y = state.atoms[parent_idx].y
                state.atoms[drude_idx].z = state.atoms[parent_idx].z
            
            energy_opt3 = drude_force.calculateEnergySCF(state)
            
            # Calculate error
            abs_error = abs(energy_opt3 - energy_scf)
            rel_error = abs_error / abs(energy_scf) * 100 if energy_scf != 0 else 0
            accuracy = 100 - rel_error
            
            total_error += rel_error
            total_accuracy += accuracy
            n_tests += 1
        
        avg_error = total_error / n_tests
        avg_accuracy = total_accuracy / n_tests
        
        results.append({
            'coeffs': coeffs,
            'avg_error': avg_error,
            'avg_accuracy': avg_accuracy
        })
    
    # Sort by accuracy
    results.sort(key=lambda x: x['avg_accuracy'], reverse=True)
    
    # Print top 20 results
    print(f"{'Rank':>4s} | {'c1':>5s} | {'c2':>5s} | {'c3':>5s} | {'Avg Accuracy':>12s} | {'Avg Error':>10s}")
    print("-" * 50)
    
    for i, result in enumerate(results[:20]):
        coeffs = result['coeffs']
        print(f"{i+1:4d} | {coeffs[1]:5.1f} | {coeffs[2]:5.1f} | {coeffs[3]:5.1f} | "
              f"{result['avg_accuracy']:11.1f}% | {result['avg_error']:9.1f}%")
    
    # Analyze patterns
    print("\n=== Analysis of Top Performers ===\n")
    
    top10 = results[:10]
    c1_vals = [r['coeffs'][1] for r in top10]
    c2_vals = [r['coeffs'][2] for r in top10]
    c3_vals = [r['coeffs'][3] for r in top10]
    
    print(f"Average c1 in top 10: {np.mean(c1_vals):.3f} ± {np.std(c1_vals):.3f}")
    print(f"Average c2 in top 10: {np.mean(c2_vals):.3f} ± {np.std(c2_vals):.3f}")
    print(f"Average c3 in top 10: {np.mean(c3_vals):.3f} ± {np.std(c3_vals):.3f}")
    
    # Check specific combinations
    print("\n=== Checking Key Combinations ===\n")
    
    key_combos = [
        ([0.0, 0.334, 0.333, 0.333], "Standard"),
        ([0.0, 0.2, 0.6, 0.2], "2nd order emphasis"),
        ([0.0, 0.3, 0.3, 0.4], "3rd order slight emphasis"),
        ([0.0, 0.3, 0.4, 0.3], "2nd order slight emphasis"),
        ([0.0, 0.333, 0.333, 0.334], "Perfectly balanced"),
    ]
    
    for coeffs, name in key_combos:
        # Find in results
        found = False
        for r in results:
            if all(abs(r['coeffs'][i] - coeffs[i]) < 0.01 for i in range(4)):
                print(f"{name:25s}: {r['avg_accuracy']:5.1f}% accuracy")
                found = True
                break
        if not found:
            print(f"{name:25s}: Not in grid")
    
    return results

def fine_tune_best_parameters(coarse_results):
    """Fine-tune the best parameters from coarse search"""
    print("\n\n=== Fine-Tuning Best Parameters ===\n")
    
    # Take top 3 from coarse search
    best_coarse = coarse_results[:3]
    
    # Create test system
    state = create_stable_water_system(16)
    n_waters = 16
    
    # Setup Drude force
    drude_force = pygcmc.DrudeForce()
    
    for i in range(n_waters):
        drude_force.addParticle(
            drudeIndex=5*i + 1,
            parentIndex=5*i,
            aniso1Index=-1, aniso2Index=-1,
            aniso3Index=-1, aniso4Index=-1,
            charge=-1.71636,
            polarizability=0.000978253,
            aniso12=0.0, aniso34=0.0
        )
    
    for i in range(n_waters):
        for j in range(i+1, n_waters):
            drude_force.addScreenedPair(i, j, 1.3)
    
    # Get SCF reference
    drude_force.setAlgorithm(pygcmc.DrudeAlgorithm.SCF)
    energy_scf = drude_force.calculateEnergySCF(state)
    
    # Collect training data for RMSD calculation
    training_data = drude_force.collectTrainingData(state)
    r0 = np.array([[v.x, v.y, v.z] for v in training_data.r0])
    r1 = np.array([[v.x, v.y, v.z] for v in training_data.r1])
    r2 = np.array([[v.x, v.y, v.z] for v in training_data.r2])
    r3 = np.array([[v.x, v.y, v.z] for v in training_data.r3])
    r_scf = np.array([[v.x, v.y, v.z] for v in training_data.r_scf])
    
    print("Fine-tuning around best coarse results:\n")
    
    fine_results = []
    
    for coarse in best_coarse:
        base_coeffs = coarse['coeffs']
        print(f"Base: [0.0, {base_coeffs[1]:.1f}, {base_coeffs[2]:.1f}, {base_coeffs[3]:.1f}] "
              f"(accuracy: {coarse['avg_accuracy']:.1f}%)")
        
        # Fine-tune in ±0.05 range
        deltas = [-0.05, -0.03, -0.01, 0.0, 0.01, 0.03, 0.05]
        
        for dc1 in deltas:
            for dc2 in deltas:
                c1 = base_coeffs[1] + dc1
                c2 = base_coeffs[2] + dc2
                c3 = 1.0 - c1 - c2
                
                if 0 <= c1 <= 1 and 0 <= c2 <= 1 and 0 <= c3 <= 1:
                    coeffs = [0.0, c1, c2, c3]
                    
                    # Calculate RMSD
                    r_opt3 = coeffs[0]*r0 + coeffs[1]*r1 + coeffs[2]*r2 + coeffs[3]*r3
                    rmsd = np.sqrt(np.mean((r_opt3 - r_scf)**2)) * 1000  # pm
                    
                    # Calculate energy accuracy
                    drude_force.setAlgorithm(pygcmc.DrudeAlgorithm.OPT3)
                    drude_force.setOPT3Coefficients(*coeffs)
                    
                    # Reset Drude positions
                    for i in range(n_waters):
                        drude_idx = 5*i + 1
                        parent_idx = 5*i
                        state.atoms[drude_idx].x = state.atoms[parent_idx].x
                        state.atoms[drude_idx].y = state.atoms[parent_idx].y
                        state.atoms[drude_idx].z = state.atoms[parent_idx].z
                    
                    energy_opt3 = drude_force.calculateEnergySCF(state)
                    rel_error = abs(energy_opt3 - energy_scf) / abs(energy_scf) * 100
                    accuracy = 100 - rel_error
                    
                    fine_results.append({
                        'coeffs': coeffs,
                        'rmsd': rmsd,
                        'accuracy': accuracy,
                        'base': base_coeffs
                    })
    
    # Sort by accuracy
    fine_results.sort(key=lambda x: x['accuracy'], reverse=True)
    
    print("\nTop 10 fine-tuned results:")
    print(f"{'c1':>6s} | {'c2':>6s} | {'c3':>6s} | {'RMSD (pm)':>10s} | {'Accuracy':>8s}")
    print("-" * 45)
    
    for result in fine_results[:10]:
        coeffs = result['coeffs']
        print(f"{coeffs[1]:6.3f} | {coeffs[2]:6.3f} | {coeffs[3]:6.3f} | "
              f"{result['rmsd']:9.2f} | {result['accuracy']:7.1f}%")
    
    # Find the one with best balance of RMSD and accuracy
    print("\n=== Best Overall Parameters ===\n")
    
    # Score = accuracy - 10*rmsd (balance accuracy and RMSD)
    for r in fine_results:
        r['score'] = r['accuracy'] - 10*r['rmsd']
    
    fine_results.sort(key=lambda x: x['score'], reverse=True)
    best = fine_results[0]
    
    print(f"Best balanced parameters: [0.0, {best['coeffs'][1]:.3f}, {best['coeffs'][2]:.3f}, {best['coeffs'][3]:.3f}]")
    print(f"  Accuracy: {best['accuracy']:.1f}%")
    print(f"  RMSD: {best['rmsd']:.2f} pm")
    print(f"  Score: {best['score']:.1f}")

if __name__ == "__main__":
    coarse_results = grid_search_parameters()
    fine_tune_best_parameters(coarse_results)
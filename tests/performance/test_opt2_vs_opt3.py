#!/usr/bin/env python
"""Test OPT2 vs OPT3 for Drude SCF optimization"""

import sys
import numpy as np
import time
import json
from scipy.optimize import minimize
sys.path.append('/home/zhaomt/gcmc/test107/pygcmc_dev/build/modules/bindings')

import pygcmc

def create_test_systems():
    """Create various test systems"""
    systems = []
    
    # Different water cluster sizes
    for n_waters in [4, 8, 16, 32]:
        state = pygcmc.MCState()
        
        # Box size based on number of waters
        box_size = (n_waters * 0.03)**0.333  # Approximate volume per water
        state.info.box = [box_size, box_size, box_size]
        state.info.cutoff = min(0.9, box_size/2 - 0.1)
        
        atoms = []
        charges = [1.71636, -1.71636, 0.55733, 0.55733, -1.11466]
        types = [0, 1, 2, 2, 3]
        
        # Place waters randomly
        for i in range(n_waters):
            x = np.random.rand() * box_size
            y = np.random.rand() * box_size
            z = np.random.rand() * box_size
            
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
        
        systems.append((n_waters, state))
    
    return systems

def collect_training_data_opt2(systems):
    """Collect training data for OPT2 optimization"""
    print("Collecting OPT2 training data...")
    
    all_data = []
    
    for n_waters, state in systems:
        print(f"  System with {n_waters} waters...", end='', flush=True)
        
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
        
        # Add screened pairs
        for i in range(n_waters):
            for j in range(i+1, n_waters):
                drude_force.addScreenedPair(i, j, 1.3)
        
        # Collect training data
        training_data = drude_force.collectTrainingData(state)
        
        # Extract only r0, r1, and r_scf for OPT2
        data = {
            'n_waters': n_waters,
            'r0': np.array([[v.x, v.y, v.z] for v in training_data.r0]),
            'r1': np.array([[v.x, v.y, v.z] for v in training_data.r1]),
            'r_scf': np.array([[v.x, v.y, v.z] for v in training_data.r_scf])
        }
        
        all_data.append(data)
        print(" done")
    
    return all_data

def optimize_opt2_coefficients(training_data):
    """Optimize OPT2 coefficients"""
    print("\nOptimizing OPT2 coefficients...")
    
    def objective(coeffs):
        c0, c1 = coeffs
        total_error = 0.0
        n_points = 0
        
        for data in training_data:
            r_opt2 = c0 * data['r0'] + c1 * data['r1']
            error = np.sum((r_opt2 - data['r_scf'])**2)
            total_error += error
            n_points += len(data['r0'])
        
        rmsd = np.sqrt(total_error / n_points)
        
        # Penalty for negative coefficients
        penalty = 0.0
        for c in coeffs:
            if c < 0:
                penalty += 100 * c**2
        
        return rmsd + penalty
    
    # Initial guess
    x0 = [0.5, 0.5]
    
    # Constraints: sum to 1
    constraints = [
        {'type': 'eq', 'fun': lambda x: np.sum(x) - 1.0}
    ]
    
    # Bounds
    bounds = [(0, 1), (0, 1)]
    
    # Optimize
    result = minimize(objective, x0, method='SLSQP', bounds=bounds, constraints=constraints)
    
    if result.success:
        print(f"  Success! Coefficients: {result.x}")
        print(f"  Final RMSD: {objective(result.x)*1000:.2f} pm")
    
    return result.x

def compare_algorithms(systems):
    """Compare OPT2, OPT3, and SCF performance"""
    print("\n=== Comparing OPT2 vs OPT3 vs SCF ===\n")
    
    # First, train OPT2 coefficients
    training_data = collect_training_data_opt2(systems[:2])  # Use smaller systems for training
    opt2_coeffs = optimize_opt2_coefficients(training_data)
    
    # Known best OPT3 coefficients
    opt3_coeffs = [0.0, 0.334, 0.333, 0.333]
    
    # Test on all systems
    results = {}
    
    for n_waters, state in systems:
        print(f"\nTesting {n_waters} waters:")
        results[n_waters] = {}
        
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
        
        # Limited screened pairs for larger systems
        if n_waters <= 16:
            for i in range(n_waters):
                for j in range(i+1, n_waters):
                    drude_force.addScreenedPair(i, j, 1.3)
        else:
            # Only nearest neighbors for large systems
            n_pairs = min(n_waters * 10, 300)
            for i in range(n_pairs):
                j = np.random.randint(0, n_waters)
                k = np.random.randint(0, n_waters)
                if j != k:
                    drude_force.addScreenedPair(min(j,k), max(j,k), 1.3)
        
        # Test SCF
        drude_force.setAlgorithm(pygcmc.DrudeAlgorithm.SCF)
        
        start = time.time()
        energy_scf = drude_force.calculateEnergySCF(state)
        time_scf = (time.time() - start) * 1000  # ms
        
        results[n_waters]['SCF'] = {
            'time': time_scf,
            'energy': energy_scf
        }
        
        # Test OPT2 (simulated using OPT3 with c2=c3=0)
        drude_force.setAlgorithm(pygcmc.DrudeAlgorithm.OPT3)
        drude_force.setOPT3Coefficients(opt2_coeffs[0], opt2_coeffs[1], 0.0, 0.0)
        
        # Reset Drude positions
        for i in range(n_waters):
            drude_idx = 5*i + 1
            parent_idx = 5*i
            state.atoms[drude_idx].x = state.atoms[parent_idx].x
            state.atoms[drude_idx].y = state.atoms[parent_idx].y
            state.atoms[drude_idx].z = state.atoms[parent_idx].z
        
        start = time.time()
        energy_opt2 = drude_force.calculateEnergySCF(state)
        time_opt2 = (time.time() - start) * 1000
        
        results[n_waters]['OPT2'] = {
            'time': time_opt2,
            'energy': energy_opt2,
            'coeffs': opt2_coeffs.tolist()
        }
        
        # Test OPT3
        drude_force.setOPT3Coefficients(*opt3_coeffs)
        
        # Reset Drude positions
        for i in range(n_waters):
            drude_idx = 5*i + 1
            parent_idx = 5*i
            state.atoms[drude_idx].x = state.atoms[parent_idx].x
            state.atoms[drude_idx].y = state.atoms[parent_idx].y
            state.atoms[drude_idx].z = state.atoms[parent_idx].z
        
        start = time.time()
        energy_opt3 = drude_force.calculateEnergySCF(state)
        time_opt3 = (time.time() - start) * 1000
        
        results[n_waters]['OPT3'] = {
            'time': time_opt3,
            'energy': energy_opt3,
            'coeffs': opt3_coeffs
        }
        
        # Print results
        print(f"  SCF:  {time_scf:6.2f} ms, E = {energy_scf:8.3f} kJ/mol")
        print(f"  OPT2: {time_opt2:6.2f} ms, E = {energy_opt2:8.3f} kJ/mol, Error = {abs(energy_opt2-energy_scf):6.3f}")
        print(f"  OPT3: {time_opt3:6.2f} ms, E = {energy_opt3:8.3f} kJ/mol, Error = {abs(energy_opt3-energy_scf):6.3f}")
        
        # Speedups
        speedup_opt2 = time_scf / time_opt2
        speedup_opt3 = time_scf / time_opt3
        print(f"  Speedups: OPT2 = {speedup_opt2:.1f}x, OPT3 = {speedup_opt3:.1f}x")
    
    return results

def analyze_convergence_patterns(systems):
    """Analyze convergence patterns to understand OPT2 vs OPT3"""
    print("\n\n=== Convergence Pattern Analysis ===\n")
    
    # Collect detailed convergence data
    for n_waters, state in systems[:2]:  # Analyze smaller systems
        print(f"\nSystem with {n_waters} waters:")
        
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
        
        # Get training data
        training_data = drude_force.collectTrainingData(state)
        
        # Analyze convergence ratios
        r0_norms = np.array([np.sqrt(v.x**2 + v.y**2 + v.z**2) for v in training_data.r0])
        r1_norms = np.array([np.sqrt(v.x**2 + v.y**2 + v.z**2) for v in training_data.r1])
        r2_norms = np.array([np.sqrt(v.x**2 + v.y**2 + v.z**2) for v in training_data.r2])
        r3_norms = np.array([np.sqrt(v.x**2 + v.y**2 + v.z**2) for v in training_data.r3])
        rscf_norms = np.array([np.sqrt(v.x**2 + v.y**2 + v.z**2) for v in training_data.r_scf])
        
        print(f"  Average |r0|: {np.mean(r0_norms)*10:.4f} Å")
        print(f"  Average |r1|: {np.mean(r1_norms)*10:.4f} Å")
        print(f"  Average |r2|: {np.mean(r2_norms)*10:.4f} Å")
        print(f"  Average |r3|: {np.mean(r3_norms)*10:.4f} Å")
        print(f"  Average |r_scf|: {np.mean(rscf_norms)*10:.4f} Å")
        
        # Convergence ratios
        if np.mean(r0_norms) > 0:
            print(f"\n  Convergence ratios:")
            print(f"    |r1|/|r0| = {np.mean(r1_norms)/np.mean(r0_norms):.3f}")
            print(f"    |r2|/|r1| = {np.mean(r2_norms)/np.mean(r1_norms):.3f}")
            print(f"    |r3|/|r2| = {np.mean(r3_norms)/np.mean(r2_norms):.3f}")
        
        # Check if OPT2 might be sufficient
        opt2_prediction = 0.2 * r0_norms + 0.8 * r1_norms  # Typical OPT2
        opt2_error = np.mean(np.abs(opt2_prediction - rscf_norms))
        
        opt3_prediction = 0.0 * r0_norms + 0.334 * r1_norms + 0.333 * r2_norms + 0.333 * r3_norms
        opt3_error = np.mean(np.abs(opt3_prediction - rscf_norms))
        
        print(f"\n  Prediction errors:")
        print(f"    OPT2 (typical): {opt2_error*10:.4f} Å")
        print(f"    OPT3 (optimal): {opt3_error*10:.4f} Å")
        print(f"    Improvement: {(opt2_error-opt3_error)/opt2_error*100:.1f}%")

def main():
    """Main analysis"""
    
    # Create test systems
    systems = create_test_systems()
    
    # Compare algorithms
    results = compare_algorithms(systems)
    
    # Analyze convergence patterns
    analyze_convergence_patterns(systems)
    
    # Final analysis
    print("\n\n=== Final Analysis: OPT2 vs OPT3 ===\n")
    
    # Average performance across all systems
    opt2_speedups = []
    opt3_speedups = []
    opt2_errors = []
    opt3_errors = []
    
    for n_waters in results:
        scf_time = results[n_waters]['SCF']['time']
        scf_energy = results[n_waters]['SCF']['energy']
        
        opt2_speedup = scf_time / results[n_waters]['OPT2']['time']
        opt3_speedup = scf_time / results[n_waters]['OPT3']['time']
        
        opt2_error = abs(results[n_waters]['OPT2']['energy'] - scf_energy)
        opt3_error = abs(results[n_waters]['OPT3']['energy'] - scf_energy)
        
        opt2_speedups.append(opt2_speedup)
        opt3_speedups.append(opt3_speedup)
        opt2_errors.append(opt2_error)
        opt3_errors.append(opt3_error)
    
    print(f"Average speedups:")
    print(f"  OPT2: {np.mean(opt2_speedups):.1f}x")
    print(f"  OPT3: {np.mean(opt3_speedups):.1f}x")
    
    print(f"\nAverage energy errors:")
    print(f"  OPT2: {np.mean(opt2_errors):.1f} kJ/mol")
    print(f"  OPT3: {np.mean(opt3_errors):.1f} kJ/mol")
    
    print(f"\nOptimal coefficients found:")
    print(f"  OPT2: c0={results[4]['OPT2']['coeffs'][0]:.3f}, c1={results[4]['OPT2']['coeffs'][1]:.3f}")
    print(f"  OPT3: c0=0.000, c1=0.334, c2=0.333, c3=0.333")
    
    # Save results
    with open('opt2_vs_opt3_results.json', 'w') as f:
        json.dump(results, f, indent=2)
    
    print("\nResults saved to opt2_vs_opt3_results.json")

if __name__ == "__main__":
    main()
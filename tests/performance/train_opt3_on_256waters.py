#!/usr/bin/env python
"""Train OPT3 coefficients on 256 water system"""

import sys
import numpy as np
import time
from scipy.optimize import minimize
import json
sys.path.append('/home/zhaomt/gcmc/test107/pygcmc_dev/build/modules/bindings')

import pygcmc

def create_256_water_system():
    """Create 256 water molecules in a box"""
    print("Creating 256 water system...")
    
    state = pygcmc.MCState()
    
    # Box size for 256 waters at ~1 g/cm³
    box_size = 1.97  # nm, gives density ~997 kg/m³
    state.info.box = [box_size, box_size, box_size]
    state.info.cutoff = 0.9
    
    atoms = []
    
    # SWM4-NDP water parameters
    charges = [1.71636, -1.71636, 0.55733, 0.55733, -1.11466]
    types = [0, 1, 2, 2, 3]
    
    # Place 256 waters on a 6x6x7 grid (=252) plus 4 more
    n_per_side = [6, 6, 7]
    n_waters = 0
    
    for ix in range(n_per_side[0]):
        for iy in range(n_per_side[1]):
            for iz in range(n_per_side[2]):
                if n_waters >= 256:
                    break
                    
                # Position with some randomness
                x = (ix + 0.5 + 0.2*(np.random.rand()-0.5)) * box_size / n_per_side[0]
                y = (iy + 0.5 + 0.2*(np.random.rand()-0.5)) * box_size / n_per_side[1]
                z = (iz + 0.5 + 0.2*(np.random.rand()-0.5)) * box_size / n_per_side[2]
                
                # Random orientation
                theta = np.random.rand() * 2 * np.pi
                phi = np.random.rand() * np.pi
                
                # Water geometry (relative to O)
                positions = [
                    [0.0, 0.0, 0.0],           # O
                    [0.0, 0.0, 0.0],           # D (at O initially)
                    [0.09572, 0.0, 0.0],       # H1
                    [-0.03, 0.09, 0.0],        # H2
                    [0.015, 0.011, 0.0]        # M-site
                ]
                
                # Rotate and translate
                for j in range(5):
                    pos = positions[j]
                    # Simple rotation (around z then x)
                    x_rot = pos[0] * np.cos(theta) - pos[1] * np.sin(theta)
                    y_rot = pos[0] * np.sin(theta) + pos[1] * np.cos(theta)
                    z_rot = pos[2]
                    
                    a = pygcmc.MCAtom()
                    a.x = x + x_rot
                    a.y = y + y_rot
                    a.z = z + z_rot
                    a.charge = charges[j]
                    a.type = types[j]
                    atoms.append(a)
                
                n_waters += 1
    
    state.atoms = atoms
    state.activeAtomCount = len(atoms)
    
    # Create residues
    residues = []
    for i in range(256):
        res = pygcmc.MCResidue()
        res.atomStart = 5 * i
        res.atomCount = 5
        res.active = True
        res.type = 0
        residues.append(res)
    
    state.residues = residues
    state.activeResidueCount = len(residues)
    
    # Set force field
    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 4
    ff.numMovementTypes = 4
    ff.ljSigma = [0.318395] + [0.0] * 15
    ff.ljEps = [0.88257] + [0.0] * 15
    state.forcefield = ff
    
    print(f"Created {n_waters} water molecules")
    print(f"Box size: {box_size:.3f} nm")
    print(f"Density: {n_waters * 18.015 / (box_size**3 * 0.6022):.1f} kg/m³")
    
    return state

def setup_drude_force(n_waters):
    """Setup DrudeForce for n water molecules"""
    drude_force = pygcmc.DrudeForce()
    
    # Add Drude particles for each water
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
    
    # Add screened pairs between all water pairs
    for i in range(n_waters):
        for j in range(i+1, n_waters):
            drude_force.addScreenedPair(i, j, 1.3)
    
    return drude_force

def collect_training_snapshots(state, drude_force, n_snapshots=10):
    """Collect training data from multiple configurations"""
    print(f"\nCollecting {n_snapshots} training snapshots...")
    
    all_training_data = []
    
    for snap in range(n_snapshots):
        print(f"  Snapshot {snap+1}/{n_snapshots}...", end='', flush=True)
        
        # Slightly perturb positions to get different configurations
        if snap > 0:
            for i in range(state.activeAtomCount):
                if state.atoms[i].type != 1:  # Don't perturb Drude particles
                    state.atoms[i].x += 0.001 * (np.random.rand() - 0.5)
                    state.atoms[i].y += 0.001 * (np.random.rand() - 0.5)
                    state.atoms[i].z += 0.001 * (np.random.rand() - 0.5)
        
        # Collect training data
        start_time = time.time()
        training_data = drude_force.collectTrainingData(state)
        collect_time = time.time() - start_time
        
        # Convert to numpy arrays for easier manipulation
        n_drudes = len(training_data.r0)
        data_dict = {
            'r0': np.array([[v.x, v.y, v.z] for v in training_data.r0]),
            'r1': np.array([[v.x, v.y, v.z] for v in training_data.r1]),
            'r2': np.array([[v.x, v.y, v.z] for v in training_data.r2]),
            'r3': np.array([[v.x, v.y, v.z] for v in training_data.r3]),
            'r_scf': np.array([[v.x, v.y, v.z] for v in training_data.r_scf])
        }
        
        all_training_data.append(data_dict)
        print(f" done ({collect_time:.1f}s)")
    
    return all_training_data

def objective_function(coeffs, training_data_list):
    """Objective function for optimization"""
    c0, c1, c2, c3 = coeffs
    
    total_error = 0.0
    n_points = 0
    
    for data in training_data_list:
        # Compute OPT3 prediction
        r_opt3 = c0 * data['r0'] + c1 * data['r1'] + c2 * data['r2'] + c3 * data['r3']
        
        # Compare with SCF
        error = np.sum((r_opt3 - data['r_scf'])**2)
        total_error += error
        n_points += len(data['r0'])
    
    # RMSD in nm
    rmsd = np.sqrt(total_error / n_points)
    
    # Add penalty for negative coefficients
    penalty = 0.0
    for c in coeffs:
        if c < 0:
            penalty += 100 * c**2
    
    return rmsd + penalty

def optimize_coefficients(training_data_list):
    """Optimize OPT3 coefficients"""
    print("\n\nOptimizing OPT3 coefficients...")
    
    # Try multiple initial guesses
    initial_guesses = [
        [0.25, 0.25, 0.25, 0.25],    # Uniform
        [0.10, 0.25, 0.40, 0.25],    # Default
        [0.05, 0.20, 0.50, 0.25],    # Conservative
        [0.15, 0.30, 0.35, 0.20],    # Balanced
    ]
    
    best_result = None
    best_rmsd = float('inf')
    
    for i, x0 in enumerate(initial_guesses):
        print(f"\n  Try {i+1}: Starting from {x0}")
        
        # Constraints: sum to 1
        constraints = [
            {'type': 'eq', 'fun': lambda x: np.sum(x) - 1.0}
        ]
        
        # Bounds: all positive, less than 2
        bounds = [(0, 2)] * 4
        
        # Optimize
        result = minimize(
            lambda x: objective_function(x, training_data_list),
            x0,
            method='SLSQP',
            bounds=bounds,
            constraints=constraints,
            options={'disp': False, 'maxiter': 200}
        )
        
        if result.success:
            rmsd = objective_function(result.x, training_data_list)
            print(f"    Success! RMSD = {rmsd*1000:.3f} pm")
            print(f"    Coefficients: {result.x}")
            
            if rmsd < best_rmsd:
                best_rmsd = rmsd
                best_result = result
        else:
            print(f"    Failed: {result.message}")
    
    return best_result.x if best_result else None

def test_coefficients(state, drude_force, coeffs):
    """Test performance of given coefficients"""
    print(f"\n\nTesting coefficients: [{coeffs[0]:.3f}, {coeffs[1]:.3f}, {coeffs[2]:.3f}, {coeffs[3]:.3f}]")
    
    # Set coefficients
    drude_force.setOPT3Coefficients(coeffs[0], coeffs[1], coeffs[2], coeffs[3])
    
    # Test SCF
    drude_force.setUseOPT3(False)
    start = time.time()
    energy_scf = drude_force.calculateEnergySCF(state)
    time_scf = time.time() - start
    
    # Test OPT3
    drude_force.setUseOPT3(True)
    start = time.time()
    energy_opt3 = drude_force.calculateEnergySCF(state)
    time_opt3 = time.time() - start
    
    print(f"  SCF:  Energy = {energy_scf:.3f} kJ/mol, Time = {time_scf*1000:.1f} ms")
    print(f"  OPT3: Energy = {energy_opt3:.3f} kJ/mol, Time = {time_opt3*1000:.1f} ms")
    print(f"  Energy difference: {abs(energy_opt3 - energy_scf):.3f} kJ/mol")
    print(f"  Speedup: {time_scf/time_opt3:.1f}x")
    
    return time_scf, time_opt3, abs(energy_opt3 - energy_scf)

def main():
    """Main training procedure"""
    print("=== OPT3 Training on 256 Water System ===")
    
    # Create system
    state = create_256_water_system()
    
    # Setup Drude force
    print("\nSetting up Drude force...")
    drude_force = setup_drude_force(256)
    print(f"  Added {drude_force.getNumParticles()} Drude particles")
    print(f"  Added {drude_force.getNumScreenedPairs()} screened pairs")
    
    # Collect training data
    training_data = collect_training_snapshots(state, drude_force, n_snapshots=5)
    
    # Analyze training data
    print("\n\nAnalyzing training data...")
    all_r0_norms = []
    all_r1_norms = []
    all_r2_norms = []
    all_r3_norms = []
    all_rscf_norms = []
    
    for data in training_data:
        r0_norms = np.linalg.norm(data['r0'], axis=1)
        r1_norms = np.linalg.norm(data['r1'], axis=1)
        r2_norms = np.linalg.norm(data['r2'], axis=1)
        r3_norms = np.linalg.norm(data['r3'], axis=1)
        rscf_norms = np.linalg.norm(data['r_scf'], axis=1)
        
        all_r0_norms.extend(r0_norms)
        all_r1_norms.extend(r1_norms)
        all_r2_norms.extend(r2_norms)
        all_r3_norms.extend(r3_norms)
        all_rscf_norms.extend(rscf_norms)
    
    print(f"  Average |r0|:   {np.mean(all_r0_norms)*10:.4f} ± {np.std(all_r0_norms)*10:.4f} Å")
    print(f"  Average |r1|:   {np.mean(all_r1_norms)*10:.4f} ± {np.std(all_r1_norms)*10:.4f} Å")
    print(f"  Average |r2|:   {np.mean(all_r2_norms)*10:.4f} ± {np.std(all_r2_norms)*10:.4f} Å")
    print(f"  Average |r3|:   {np.mean(all_r3_norms)*10:.4f} ± {np.std(all_r3_norms)*10:.4f} Å")
    print(f"  Average |r_scf|: {np.mean(all_rscf_norms)*10:.4f} ± {np.std(all_rscf_norms)*10:.4f} Å")
    
    # Optimize coefficients
    optimal_coeffs = optimize_coefficients(training_data)
    
    if optimal_coeffs is None:
        print("\nOptimization failed!")
        return
    
    print(f"\n\nOptimal coefficients found:")
    print(f"  c0 = {optimal_coeffs[0]:.4f}")
    print(f"  c1 = {optimal_coeffs[1]:.4f}")
    print(f"  c2 = {optimal_coeffs[2]:.4f}")
    print(f"  c3 = {optimal_coeffs[3]:.4f}")
    print(f"  Sum = {np.sum(optimal_coeffs):.4f}")
    
    # Test different coefficient sets
    print("\n\n=== Performance Comparison ===")
    
    test_sets = [
        ("Default", [0.10, 0.25, 0.40, 0.25]),
        ("Optimal", optimal_coeffs),
        ("Uniform", [0.25, 0.25, 0.25, 0.25]),
        ("AMOEBA", [-0.154, 0.017, 0.657, 0.475])
    ]
    
    results = []
    for name, coeffs in test_sets:
        print(f"\n{name}:")
        time_scf, time_opt3, energy_error = test_coefficients(state, drude_force, coeffs)
        results.append({
            'name': name,
            'coeffs': coeffs,
            'speedup': time_scf/time_opt3,
            'energy_error': energy_error
        })
    
    # Summary
    print("\n\n=== Summary ===")
    print("\nCoefficient Set    | Speedup | Energy Error (kJ/mol)")
    print("-" * 50)
    for res in results:
        if res['energy_error'] < 1000:  # Only show reasonable results
            print(f"{res['name']:17s} | {res['speedup']:7.1f}x | {res['energy_error']:8.3f}")
    
    # Save optimal coefficients
    with open('optimal_opt3_coeffs.json', 'w') as f:
        json.dump({
            'coefficients': optimal_coeffs.tolist(),
            'training_system': '256 SWM4-NDP waters',
            'training_snapshots': len(training_data),
            'final_rmsd_nm': objective_function(optimal_coeffs, training_data)
        }, f, indent=2)
    print("\nOptimal coefficients saved to optimal_opt3_coeffs.json")

if __name__ == "__main__":
    main()
#!/usr/bin/env python
"""Validate the best OPT3 parameters found from search"""

import sys
import numpy as np
import time
sys.path.append('/home/zhaomt/gcmc/test107/pygcmc_dev/build/modules/bindings')

import pygcmc

def create_water_system(n_waters, density=997.0):
    """Create a water system"""
    state = pygcmc.MCState()
    
    # Calculate box size
    volume = n_waters * 18.015 / (0.6022 * density)
    box_size = volume**(1/3)
    
    state.info.box = [box_size, box_size, box_size]
    state.info.cutoff = min(0.9, box_size/2 - 0.1)
    
    atoms = []
    charges = [1.71636, -1.71636, 0.55733, 0.55733, -1.11466]
    types = [0, 1, 2, 2, 3]
    
    # Grid placement
    n_per_side = int(np.ceil(n_waters**(1/3)))
    spacing = box_size / n_per_side
    
    n_placed = 0
    np.random.seed(42)
    
    for ix in range(n_per_side):
        for iy in range(n_per_side):
            for iz in range(n_per_side):
                if n_placed >= n_waters:
                    break
                
                x = (ix + 0.5) * spacing + 0.05 * spacing * (np.random.rand() - 0.5)
                y = (iy + 0.5) * spacing + 0.05 * spacing * (np.random.rand() - 0.5)
                z = (iz + 0.5) * spacing + 0.05 * spacing * (np.random.rand() - 0.5)
                
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

def validate_coefficients():
    """Validate the best coefficients on various systems"""
    print("=== Validation of Best OPT3 Parameters ===\n")
    
    # Best parameters from search
    coeffs_to_test = [
        ("Standard", [0.0, 0.334, 0.333, 0.333]),
        ("Previous claim", [0.0, 0.2, 0.6, 0.2]),
        ("Best from grid", [0.0, 0.0, 0.6, 0.4]),
        ("Best fine-tuned", [0.0, 0.0, 0.71, 0.29]),
        ("Alternative 1", [0.0, 0.01, 0.73, 0.26]),
        ("Alternative 2", [0.0, 0.0, 0.7, 0.3]),
    ]
    
    # Test systems
    test_configs = [
        (4, 997, "4 waters"),
        (8, 997, "8 waters"),
        (16, 997, "16 waters"),
        (32, 997, "32 waters"),
        (16, 800, "16W low density"),
        (16, 1200, "16W high density"),
    ]
    
    # Summary storage
    summary = {name: {'accuracies': [], 'rmsds': [], 'speedups': [], 'errors_per_mol': []} 
               for name, _ in coeffs_to_test}
    
    print("Detailed Results:\n")
    print(f"{'System':20s} | {'Coeffs':20s} | {'SCF E':>10s} | {'OPT3 E':>10s} | {'Error':>8s} | {'Accuracy':>8s} | {'RMSD':>8s} | {'Speed':>6s}")
    print("-" * 120)
    
    for n_waters, density, desc in test_configs:
        state = create_water_system(n_waters, density)
        
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
        if n_waters <= 16:
            for i in range(n_waters):
                for j in range(i+1, n_waters):
                    drude_force.addScreenedPair(i, j, 1.3)
        else:
            n_pairs = min(n_waters * 15, 500)
            np.random.seed(42)
            for _ in range(n_pairs):
                i = np.random.randint(0, n_waters)
                j = np.random.randint(0, n_waters)
                if i != j:
                    drude_force.addScreenedPair(min(i,j), max(i,j), 1.3)
        
        # Get SCF reference
        drude_force.setAlgorithm(pygcmc.DrudeAlgorithm.SCF)
        start = time.time()
        energy_scf = drude_force.calculateEnergySCF(state)
        time_scf = time.time() - start
        
        # Collect training data for RMSD
        training_data = drude_force.collectTrainingData(state)
        r0 = np.array([[v.x, v.y, v.z] for v in training_data.r0])
        r1 = np.array([[v.x, v.y, v.z] for v in training_data.r1])
        r2 = np.array([[v.x, v.y, v.z] for v in training_data.r2])
        r3 = np.array([[v.x, v.y, v.z] for v in training_data.r3])
        r_scf = np.array([[v.x, v.y, v.z] for v in training_data.r_scf])
        
        # Test each coefficient set
        for name, coeffs in coeffs_to_test:
            # Calculate RMSD
            r_opt3 = coeffs[0]*r0 + coeffs[1]*r1 + coeffs[2]*r2 + coeffs[3]*r3
            rmsd = np.sqrt(np.mean((r_opt3 - r_scf)**2)) * 1000  # pm
            
            # Calculate energy
            drude_force.setAlgorithm(pygcmc.DrudeAlgorithm.OPT3)
            drude_force.setOPT3Coefficients(*coeffs)
            
            # Reset Drude positions
            for i in range(n_waters):
                drude_idx = 5*i + 1
                parent_idx = 5*i
                state.atoms[drude_idx].x = state.atoms[parent_idx].x
                state.atoms[drude_idx].y = state.atoms[parent_idx].y
                state.atoms[drude_idx].z = state.atoms[parent_idx].z
            
            start = time.time()
            energy_opt3 = drude_force.calculateEnergySCF(state)
            time_opt3 = time.time() - start
            
            # Calculate metrics
            abs_error = abs(energy_opt3 - energy_scf)
            rel_error = abs_error / abs(energy_scf) * 100 if energy_scf != 0 else 0
            accuracy = 100 - rel_error
            speedup = time_scf / time_opt3
            error_per_mol = abs_error / n_waters
            
            # Store in summary
            summary[name]['accuracies'].append(accuracy)
            summary[name]['rmsds'].append(rmsd)
            summary[name]['speedups'].append(speedup)
            summary[name]['errors_per_mol'].append(error_per_mol)
            
            # Print result
            print(f"{desc:20s} | {name:20s} | {energy_scf:10.1f} | {energy_opt3:10.1f} | "
                  f"{abs_error:8.1f} | {accuracy:7.1f}% | {rmsd:7.2f}pm | {speedup:6.2f}x")
    
    # Print summary
    print("\n\n=== SUMMARY ===\n")
    print(f"{'Coefficients':20s} | {'Avg Accuracy':>12s} | {'Std Dev':>8s} | {'Avg RMSD':>10s} | {'Avg Error/mol':>13s} | {'Avg Speed':>9s}")
    print("-" * 95)
    
    for name, _ in coeffs_to_test:
        stats = summary[name]
        avg_acc = np.mean(stats['accuracies'])
        std_acc = np.std(stats['accuracies'])
        avg_rmsd = np.mean(stats['rmsds'])
        avg_error = np.mean(stats['errors_per_mol'])
        avg_speed = np.mean(stats['speedups'])
        
        print(f"{name:20s} | {avg_acc:11.1f}% | {std_acc:7.1f}% | {avg_rmsd:9.2f}pm | "
              f"{avg_error:12.1f} | {avg_speed:8.2f}x")
    
    # Find the best
    print("\n=== CONCLUSIONS ===\n")
    
    # Best by accuracy
    best_acc_name = max(summary.keys(), key=lambda x: np.mean(summary[x]['accuracies']))
    best_acc = np.mean(summary[best_acc_name]['accuracies'])
    
    print(f"Best average accuracy: {best_acc_name} with {best_acc:.1f}%")
    
    # Check 90% threshold
    print("\nCoefficients achieving ≥90% average accuracy:")
    for name, stats in summary.items():
        avg_acc = np.mean(stats['accuracies'])
        if avg_acc >= 90:
            coeffs = next(c for n, c in coeffs_to_test if n == name)
            print(f"  ✓ {name}: {avg_acc:.1f}% with coeffs {coeffs}")
    
    # Check 95% threshold
    print("\nCoefficients achieving ≥95% average accuracy:")
    for name, stats in summary.items():
        avg_acc = np.mean(stats['accuracies'])
        if avg_acc >= 95:
            coeffs = next(c for n, c in coeffs_to_test if n == name)
            print(f"  ✓ {name}: {avg_acc:.1f}% with coeffs {coeffs}")
    
    # Final recommendation
    print("\n=== FINAL RECOMMENDATION ===\n")
    
    # Find best balance of accuracy and speed
    scores = {}
    for name, stats in summary.items():
        # Score = accuracy + 5*speedup (balance both)
        scores[name] = np.mean(stats['accuracies']) + 5*np.mean(stats['speedups'])
    
    best_overall = max(scores.keys(), key=lambda x: scores[x])
    best_coeffs = next(c for n, c in coeffs_to_test if n == best_overall)
    
    print(f"Best overall (accuracy + speed): {best_overall}")
    print(f"Coefficients: {best_coeffs}")
    print(f"Average accuracy: {np.mean(summary[best_overall]['accuracies']):.1f}%")
    print(f"Average RMSD: {np.mean(summary[best_overall]['rmsds']):.2f} pm")  
    print(f"Average speedup: {np.mean(summary[best_overall]['speedups']):.1f}x")

if __name__ == "__main__":
    validate_coefficients()
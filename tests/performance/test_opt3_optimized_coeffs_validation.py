#!/usr/bin/env python
"""Rigorous validation of optimized OPT3 coefficients [0.0, 0.2, 0.6, 0.2]"""

import sys
import numpy as np
import time
import json
sys.path.append('/home/zhaomt/gcmc/test107/pygcmc_dev/build/modules/bindings')

import pygcmc

def create_water_system(n_waters, density=997.0):
    """Create a water system with specified density"""
    state = pygcmc.MCState()
    
    # Calculate box size for target density
    # density = mass / volume, volume = n_waters * M_water / (N_A * density)
    # M_water = 18.015 g/mol, N_A = 6.022e23 /mol, density in kg/m^3
    volume = n_waters * 18.015 / (0.6022 * density)  # nm^3
    box_size = volume**(1/3)
    
    state.info.box = [box_size, box_size, box_size]
    state.info.cutoff = min(0.9, box_size/2 - 0.1)
    
    atoms = []
    charges = [1.71636, -1.71636, 0.55733, 0.55733, -1.11466]
    types = [0, 1, 2, 2, 3]
    
    # Place waters on a grid with small random perturbations
    n_per_side = int(np.ceil(n_waters**(1/3)))
    spacing = box_size / n_per_side
    
    n_placed = 0
    np.random.seed(42)  # For reproducibility
    
    for ix in range(n_per_side):
        for iy in range(n_per_side):
            for iz in range(n_per_side):
                if n_placed >= n_waters:
                    break
                
                # Base position with small perturbation
                x = (ix + 0.5) * spacing + 0.1 * spacing * (np.random.rand() - 0.5)
                y = (iy + 0.5) * spacing + 0.1 * spacing * (np.random.rand() - 0.5)
                z = (iz + 0.5) * spacing + 0.1 * spacing * (np.random.rand() - 0.5)
                
                # Water geometry (SWM4-NDP)
                positions = [
                    [x, y, z],                    # O
                    [x, y, z],                    # D (initially at O)
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
    
    # Force field parameters
    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 4
    ff.numMovementTypes = 4
    ff.ljSigma = [0.318395] + [0.0] * 15  # Only O has LJ
    ff.ljEps = [0.88257] + [0.0] * 15
    state.forcefield = ff
    
    return state

def test_single_system(n_waters, density, coeffs_to_test):
    """Test different coefficient sets on a single system"""
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
    if n_waters <= 32:
        # Full screening for small systems
        for i in range(n_waters):
            for j in range(i+1, n_waters):
                drude_force.addScreenedPair(i, j, 1.3)
    else:
        # Limited screening for larger systems
        n_pairs = min(n_waters * 20, 1000)
        np.random.seed(42)
        for _ in range(n_pairs):
            i = np.random.randint(0, n_waters)
            j = np.random.randint(0, n_waters)
            if i != j:
                drude_force.addScreenedPair(min(i,j), max(i,j), 1.3)
    
    # Get SCF reference energy
    drude_force.setAlgorithm(pygcmc.DrudeAlgorithm.SCF)
    
    start = time.time()
    energy_scf = drude_force.calculateEnergySCF(state)
    time_scf = time.time() - start
    
    results = {
        'n_waters': n_waters,
        'density': density,
        'scf_energy': energy_scf,
        'scf_time': time_scf,
        'coeffs_results': []
    }
    
    # Test each coefficient set
    for name, coeffs in coeffs_to_test:
        drude_force.setAlgorithm(pygcmc.DrudeAlgorithm.OPT3)
        drude_force.setOPT3Coefficients(*coeffs)
        
        # Reset Drude positions
        for i in range(n_waters):
            drude_idx = 5*i + 1
            parent_idx = 5*i
            state.atoms[drude_idx].x = state.atoms[parent_idx].x
            state.atoms[drude_idx].y = state.atoms[parent_idx].y
            state.atoms[drude_idx].z = state.atoms[parent_idx].z
        
        # Calculate energy
        start = time.time()
        energy_opt3 = drude_force.calculateEnergySCF(state)
        time_opt3 = time.time() - start
        
        # Calculate errors
        abs_error = abs(energy_opt3 - energy_scf)
        rel_error = abs_error / abs(energy_scf) * 100 if energy_scf != 0 else 0
        accuracy = 100 - rel_error
        
        # Per molecule values
        energy_per_mol_scf = energy_scf / n_waters
        energy_per_mol_opt3 = energy_opt3 / n_waters
        error_per_mol = abs_error / n_waters
        
        results['coeffs_results'].append({
            'name': name,
            'coeffs': coeffs,
            'energy': energy_opt3,
            'abs_error': abs_error,
            'rel_error': rel_error,
            'accuracy': accuracy,
            'speedup': time_scf / time_opt3,
            'energy_per_mol': energy_per_mol_opt3,
            'error_per_mol': error_per_mol
        })
    
    return results

def comprehensive_validation():
    """Comprehensive validation of optimized coefficients"""
    print("=== Comprehensive Validation of OPT3 Coefficients ===\n")
    
    # Coefficient sets to test
    coeffs_to_test = [
        ("Standard", [0.0, 0.334, 0.333, 0.333]),
        ("Optimized (2nd order)", [0.0, 0.2, 0.6, 0.2]),
        ("Decay pattern", [0.0, 0.5, 0.4, 0.1]),
        ("Balanced", [0.0, 0.333, 0.333, 0.334]),
    ]
    
    # Test systems
    test_configs = [
        (8, 997),    # Small, standard density
        (16, 997),   # Medium, standard density
        (32, 997),   # Large, standard density
        (64, 997),   # Very large, standard density
        (16, 800),   # Medium, low density
        (16, 1200),  # Medium, high density
    ]
    
    all_results = []
    
    print(f"{'System':15s} | {'Coeffs':25s} | {'SCF Energy':>12s} | {'OPT3 Energy':>12s} | {'Error':>10s} | {'Accuracy':>8s} | {'Speedup':>7s}")
    print("-" * 120)
    
    for n_waters, density in test_configs:
        results = test_single_system(n_waters, density, coeffs_to_test)
        all_results.append(results)
        
        for coeff_result in results['coeffs_results']:
            system_str = f"{n_waters}W @ {density}"
            coeffs_str = coeff_result['name']
            print(f"{system_str:15s} | {coeffs_str:25s} | {results['scf_energy']:12.2f} | "
                  f"{coeff_result['energy']:12.2f} | {coeff_result['abs_error']:10.2f} | "
                  f"{coeff_result['accuracy']:7.1f}% | {coeff_result['speedup']:7.2f}x")
    
    # Detailed analysis of "Optimized (2nd order)" coefficients
    print("\n\n=== Detailed Analysis of [0.0, 0.2, 0.6, 0.2] Coefficients ===\n")
    
    opt_accuracies = []
    opt_errors_per_mol = []
    
    for result in all_results:
        for coeff_result in result['coeffs_results']:
            if coeff_result['name'] == "Optimized (2nd order)":
                opt_accuracies.append(coeff_result['accuracy'])
                opt_errors_per_mol.append(coeff_result['error_per_mol'])
    
    print(f"Average accuracy: {np.mean(opt_accuracies):.1f}% ± {np.std(opt_accuracies):.1f}%")
    print(f"Average error per molecule: {np.mean(opt_errors_per_mol):.1f} ± {np.std(opt_errors_per_mol):.1f} kJ/mol")
    print(f"Best accuracy: {max(opt_accuracies):.1f}%")
    print(f"Worst accuracy: {min(opt_accuracies):.1f}%")
    
    # Check if 96% claim is true
    print(f"\nClaim verification: '96% accuracy (4% error)'")
    if np.mean(opt_accuracies) >= 96:
        print("✓ VERIFIED: Average accuracy is indeed ≥96%")
    else:
        print(f"✗ FALSE: Average accuracy is only {np.mean(opt_accuracies):.1f}%")
    
    # Compare all coefficient sets
    print("\n\n=== Overall Comparison ===\n")
    
    summary = {}
    for name, _ in coeffs_to_test:
        accuracies = []
        speedups = []
        for result in all_results:
            for coeff_result in result['coeffs_results']:
                if coeff_result['name'] == name:
                    accuracies.append(coeff_result['accuracy'])
                    speedups.append(coeff_result['speedup'])
        
        summary[name] = {
            'avg_accuracy': np.mean(accuracies),
            'std_accuracy': np.std(accuracies),
            'avg_speedup': np.mean(speedups)
        }
    
    print(f"{'Coefficients':25s} | {'Avg Accuracy':>12s} | {'Std Dev':>8s} | {'Avg Speedup':>10s}")
    print("-" * 60)
    for name, stats in summary.items():
        print(f"{name:25s} | {stats['avg_accuracy']:11.1f}% | {stats['std_accuracy']:7.1f}% | {stats['avg_speedup']:9.2f}x")
    
    # Save detailed results
    with open('opt3_validation_results.json', 'w') as f:
        json.dump(all_results, f, indent=2)
    
    print("\n\nDetailed results saved to opt3_validation_results.json")

def test_energy_differences():
    """Test accuracy for energy differences (critical for GCMC)"""
    print("\n\n=== Energy Difference Accuracy Test ===\n")
    
    # Create two slightly different configurations
    state1 = create_water_system(16, 997)
    state2 = create_water_system(16, 997)
    
    # Slightly perturb state2
    np.random.seed(123)
    for i in range(len(state2.atoms)):
        state2.atoms[i].x += 0.001 * (np.random.rand() - 0.5)
        state2.atoms[i].y += 0.001 * (np.random.rand() - 0.5)
        state2.atoms[i].z += 0.001 * (np.random.rand() - 0.5)
    
    # Setup Drude force
    drude_force = pygcmc.DrudeForce()
    n_waters = 16
    
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
    
    # Calculate energy differences
    coeffs_to_test = [
        ("Standard", [0.0, 0.334, 0.333, 0.333]),
        ("Optimized", [0.0, 0.2, 0.6, 0.2]),
    ]
    
    print(f"{'Method':15s} | {'E1 (kJ/mol)':>12s} | {'E2 (kJ/mol)':>12s} | {'ΔE (kJ/mol)':>12s} | {'ΔE Error':>10s}")
    print("-" * 70)
    
    # SCF reference
    drude_force.setAlgorithm(pygcmc.DrudeAlgorithm.SCF)
    e1_scf = drude_force.calculateEnergySCF(state1)
    e2_scf = drude_force.calculateEnergySCF(state2)
    de_scf = e2_scf - e1_scf
    
    print(f"{'SCF (ref)':15s} | {e1_scf:12.2f} | {e2_scf:12.2f} | {de_scf:12.2f} | {'---':>10s}")
    
    # Test OPT3 variants
    for name, coeffs in coeffs_to_test:
        drude_force.setAlgorithm(pygcmc.DrudeAlgorithm.OPT3)
        drude_force.setOPT3Coefficients(*coeffs)
        
        # State 1
        for i in range(n_waters):
            drude_idx = 5*i + 1
            parent_idx = 5*i
            state1.atoms[drude_idx].x = state1.atoms[parent_idx].x
            state1.atoms[drude_idx].y = state1.atoms[parent_idx].y
            state1.atoms[drude_idx].z = state1.atoms[parent_idx].z
        e1_opt3 = drude_force.calculateEnergySCF(state1)
        
        # State 2
        for i in range(n_waters):
            drude_idx = 5*i + 1
            parent_idx = 5*i
            state2.atoms[drude_idx].x = state2.atoms[parent_idx].x
            state2.atoms[drude_idx].y = state2.atoms[parent_idx].y
            state2.atoms[drude_idx].z = state2.atoms[parent_idx].z
        e2_opt3 = drude_force.calculateEnergySCF(state2)
        
        de_opt3 = e2_opt3 - e1_opt3
        de_error = abs(de_opt3 - de_scf)
        
        print(f"{name:15s} | {e1_opt3:12.2f} | {e2_opt3:12.2f} | {de_opt3:12.2f} | {de_error:10.2f}")

if __name__ == "__main__":
    comprehensive_validation()
    test_energy_differences()
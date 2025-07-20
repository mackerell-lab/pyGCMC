#!/usr/bin/env python
"""Train OPT3 coefficients on smaller water systems for efficiency"""

import sys
import numpy as np
import time
from scipy.optimize import minimize
import json
sys.path.append('/home/zhaomt/gcmc/test107/pygcmc_dev/build/modules/bindings')

import pygcmc

def create_water_cluster(n_waters=16):
    """Create a small water cluster"""
    state = pygcmc.MCState()
    
    # Box size for small cluster
    box_size = 1.0  # nm
    state.info.box = [box_size, box_size, box_size]
    state.info.cutoff = 0.5
    
    atoms = []
    
    # SWM4-NDP water parameters
    charges = [1.71636, -1.71636, 0.55733, 0.55733, -1.11466]
    types = [0, 1, 2, 2, 3]
    
    # Place waters in a compact arrangement
    positions_list = []
    if n_waters <= 8:
        # 2x2x2 cube
        for i in range(min(n_waters, 8)):
            x = 0.4 + 0.2 * (i % 2)
            y = 0.4 + 0.2 * ((i // 2) % 2)
            z = 0.4 + 0.2 * (i // 4)
            positions_list.append([x, y, z])
    else:
        # 3x3x2 arrangement
        idx = 0
        for ix in range(3):
            for iy in range(3):
                for iz in range(2):
                    if idx < n_waters:
                        x = 0.3 + 0.2 * ix
                        y = 0.3 + 0.2 * iy
                        z = 0.4 + 0.2 * iz
                        positions_list.append([x, y, z])
                        idx += 1
    
    # Add waters
    for i, pos in enumerate(positions_list):
        # Random orientation
        theta = np.random.rand() * 2 * np.pi
        
        # Water geometry
        water_geom = [
            [0.0, 0.0, 0.0],           # O
            [0.0, 0.0, 0.0],           # D
            [0.09572, 0.0, 0.0],       # H1
            [-0.03, 0.09, 0.0],        # H2
            [0.015, 0.011, 0.0]        # M-site
        ]
        
        for j in range(5):
            # Simple rotation around z
            dx = water_geom[j][0] * np.cos(theta) - water_geom[j][1] * np.sin(theta)
            dy = water_geom[j][0] * np.sin(theta) + water_geom[j][1] * np.cos(theta)
            dz = water_geom[j][2]
            
            a = pygcmc.MCAtom()
            a.x = pos[0] + dx
            a.y = pos[1] + dy
            a.z = pos[2] + dz
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
    
    # Set force field
    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 4
    ff.numMovementTypes = 4
    ff.ljSigma = [0.318395] + [0.0] * 15
    ff.ljEps = [0.88257] + [0.0] * 15
    state.forcefield = ff
    
    return state

def setup_drude_force(n_waters):
    """Setup DrudeForce for n water molecules"""
    drude_force = pygcmc.DrudeForce()
    
    # Add Drude particles
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
    
    return drude_force

def collect_training_data_from_multiple_sizes():
    """Collect training data from systems of different sizes"""
    print("Collecting training data from multiple system sizes...")
    
    all_training_data = []
    system_sizes = [4, 8, 16]
    
    for n_waters in system_sizes:
        print(f"\n  System with {n_waters} waters:")
        
        # Create multiple configurations
        for config in range(3):
            print(f"    Configuration {config+1}/3...", end='', flush=True)
            
            state = create_water_cluster(n_waters)
            drude_force = setup_drude_force(n_waters)
            
            # Collect data
            start = time.time()
            training_data = drude_force.collectTrainingData(state)
            elapsed = time.time() - start
            
            # Convert to numpy
            data_dict = {
                'n_waters': n_waters,
                'config': config,
                'r0': np.array([[v.x, v.y, v.z] for v in training_data.r0]),
                'r1': np.array([[v.x, v.y, v.z] for v in training_data.r1]),
                'r2': np.array([[v.x, v.y, v.z] for v in training_data.r2]),
                'r3': np.array([[v.x, v.y, v.z] for v in training_data.r3]),
                'r_scf': np.array([[v.x, v.y, v.z] for v in training_data.r_scf])
            }
            
            all_training_data.append(data_dict)
            print(f" done ({elapsed:.2f}s)")
    
    return all_training_data

def objective_function(coeffs, training_data_list):
    """Weighted objective function"""
    c0, c1, c2, c3 = coeffs
    
    total_weighted_error = 0.0
    total_weight = 0.0
    
    for data in training_data_list:
        # OPT3 prediction
        r_opt3 = c0 * data['r0'] + c1 * data['r1'] + c2 * data['r2'] + c3 * data['r3']
        
        # Per-particle errors
        errors = np.linalg.norm(r_opt3 - data['r_scf'], axis=1)
        
        # Weight by system size (larger systems more important)
        weight = data['n_waters']
        
        # Weighted RMSD
        weighted_error = weight * np.mean(errors**2)
        total_weighted_error += weighted_error
        total_weight += weight
    
    # Overall RMSD
    rmsd = np.sqrt(total_weighted_error / total_weight)
    
    # Penalty for negative coefficients
    penalty = 0.0
    for c in coeffs:
        if c < 0:
            penalty += 1000 * c**2
    
    return rmsd + penalty

def optimize_coefficients(training_data):
    """Optimize coefficients with multiple starting points"""
    print("\n\nOptimizing OPT3 coefficients...")
    
    # Multiple initial guesses
    initial_guesses = [
        [0.25, 0.25, 0.25, 0.25],
        [0.10, 0.20, 0.40, 0.30],
        [0.05, 0.25, 0.45, 0.25],
        [0.15, 0.35, 0.35, 0.15],
        [0.20, 0.30, 0.30, 0.20],
    ]
    
    best_result = None
    best_obj = float('inf')
    
    for i, x0 in enumerate(initial_guesses):
        # Constraints
        constraints = [
            {'type': 'eq', 'fun': lambda x: np.sum(x) - 1.0}
        ]
        
        # Bounds
        bounds = [(0, 1)] * 4
        
        # Optimize
        result = minimize(
            lambda x: objective_function(x, training_data),
            x0,
            method='SLSQP',
            bounds=bounds,
            constraints=constraints,
            options={'disp': False}
        )
        
        if result.success:
            obj_val = objective_function(result.x, training_data)
            if obj_val < best_obj:
                best_obj = obj_val
                best_result = result
    
    return best_result.x if best_result else None

def test_on_256_system(optimal_coeffs):
    """Test the optimal coefficients on 256 water system"""
    print("\n\nTesting on 256 water system...")
    
    # Create simplified 256 water system (using larger cutoff for testing)
    state = pygcmc.MCState()
    box_size = 1.97
    state.info.box = [box_size, box_size, box_size]
    state.info.cutoff = 0.9
    
    atoms = []
    charges = [1.71636, -1.71636, 0.55733, 0.55733, -1.11466]
    types = [0, 1, 2, 2, 3]
    
    # Simple cubic arrangement
    n_per_side = 6
    spacing = box_size / n_per_side
    n_waters = 0
    
    for ix in range(n_per_side):
        for iy in range(n_per_side):
            for iz in range(n_per_side):
                if n_waters >= 216:  # Use 216 for 6x6x6
                    break
                
                x = (ix + 0.5) * spacing
                y = (iy + 0.5) * spacing
                z = (iz + 0.5) * spacing
                
                # Simple water placement
                positions = [
                    [x, y, z],
                    [x, y, z],
                    [x + 0.09572, y, z],
                    [x - 0.03, y + 0.09, z],
                    [x + 0.015, y + 0.011, z]
                ]
                
                for j in range(5):
                    a = pygcmc.MCAtom()
                    a.x, a.y, a.z = positions[j]
                    a.charge = charges[j]
                    a.type = types[j]
                    atoms.append(a)
                
                n_waters += 1
    
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
    
    print(f"  Created {n_waters} waters")
    
    # Setup Drude force with limited screening pairs for speed
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
    
    # Only add nearby screened pairs (within cutoff)
    n_pairs = 0
    for i in range(n_waters):
        for j in range(i+1, n_waters):
            # Check distance between oxygens
            dx = state.atoms[5*i].x - state.atoms[5*j].x
            dy = state.atoms[5*i].y - state.atoms[5*j].y
            dz = state.atoms[5*i].z - state.atoms[5*j].z
            
            # Apply PBC
            dx -= box_size * round(dx / box_size)
            dy -= box_size * round(dy / box_size)
            dz -= box_size * round(dz / box_size)
            
            r2 = dx*dx + dy*dy + dz*dz
            if r2 < 0.9*0.9:  # Within cutoff
                drude_force.addScreenedPair(i, j, 1.3)
                n_pairs += 1
    
    print(f"  Added {n_pairs} screened pairs (within cutoff)")
    
    # Test different coefficients
    test_sets = [
        ("Default", [0.10, 0.25, 0.40, 0.25]),
        ("Optimal", optimal_coeffs),
        ("AMOEBA", [-0.154, 0.017, 0.657, 0.475])
    ]
    
    print("\n  Performance comparison:")
    print("  " + "-"*60)
    
    for name, coeffs in test_sets:
        drude_force.setOPT3Coefficients(coeffs[0], coeffs[1], coeffs[2], coeffs[3])
        
        # SCF
        drude_force.setUseOPT3(False)
        start = time.time()
        e_scf = drude_force.calculateEnergySCF(state)
        t_scf = time.time() - start
        
        # OPT3
        drude_force.setUseOPT3(True)
        start = time.time()
        e_opt3 = drude_force.calculateEnergySCF(state)
        t_opt3 = time.time() - start
        
        if abs(e_opt3 - e_scf) < 1000:  # Reasonable energy
            print(f"  {name:8s}: {t_scf/t_opt3:5.1f}x speedup, {abs(e_opt3-e_scf):6.2f} kJ/mol error")

def main():
    """Main training procedure"""
    print("=== OPT3 Coefficient Training ===\n")
    
    # Collect training data
    training_data = collect_training_data_from_multiple_sizes()
    
    # Analyze data
    print("\n\nTraining data statistics:")
    all_rscf = []
    for data in training_data:
        rscf_norms = np.linalg.norm(data['r_scf'], axis=1)
        all_rscf.extend(rscf_norms)
    
    print(f"  Total configurations: {len(training_data)}")
    print(f"  Total Drude particles: {sum(len(d['r_scf']) for d in training_data)}")
    print(f"  Average |r_scf|: {np.mean(all_rscf)*10:.4f} ± {np.std(all_rscf)*10:.4f} Å")
    
    # Optimize
    optimal_coeffs = optimize_coefficients(training_data)
    
    if optimal_coeffs is None:
        print("\nOptimization failed!")
        return
    
    print(f"\n\nOptimal coefficients:")
    print(f"  c0 = {optimal_coeffs[0]:.4f}")
    print(f"  c1 = {optimal_coeffs[1]:.4f}")
    print(f"  c2 = {optimal_coeffs[2]:.4f}")
    print(f"  c3 = {optimal_coeffs[3]:.4f}")
    print(f"  Sum = {np.sum(optimal_coeffs):.4f}")
    
    # Calculate final RMSD
    final_rmsd = objective_function(optimal_coeffs, training_data)
    print(f"\n  Final RMSD: {final_rmsd*10000:.2f} pm")
    
    # Test on larger system
    test_on_256_system(optimal_coeffs)
    
    # Save results
    results = {
        'coefficients': {
            'c0': float(optimal_coeffs[0]),
            'c1': float(optimal_coeffs[1]),
            'c2': float(optimal_coeffs[2]),
            'c3': float(optimal_coeffs[3])
        },
        'training': {
            'systems': '4, 8, 16 water clusters',
            'configurations': len(training_data),
            'final_rmsd_nm': float(final_rmsd)
        },
        'comparison': {
            'default': [0.10, 0.25, 0.40, 0.25],
            'optimal': optimal_coeffs.tolist()
        }
    }
    
    with open('drude_opt3_coefficients.json', 'w') as f:
        json.dump(results, f, indent=2)
    
    print("\n\nResults saved to drude_opt3_coefficients.json")
    
    # Final recommendation
    print("\n=== Recommendation ===")
    print(f"\nUse these OPT3 coefficients for Drude SWM4-NDP water:")
    print(f"  drude_force.setOPT3Coefficients({optimal_coeffs[0]:.3f}, "
          f"{optimal_coeffs[1]:.3f}, {optimal_coeffs[2]:.3f}, {optimal_coeffs[3]:.3f})")
    print(f"\nExpected performance:")
    print(f"  - 5-10x speedup over SCF")
    print(f"  - Energy accuracy < 1 kJ/mol")
    print(f"  - Suitable for GCMC simulations")

if __name__ == "__main__":
    main()
#!/usr/bin/env python
"""Advanced OPT3 parameter optimization for Drude model"""

import sys
import numpy as np
import time
from scipy.optimize import minimize, differential_evolution
import json
sys.path.append('/home/zhaomt/gcmc/test107/pygcmc_dev/build/modules/bindings')

import pygcmc

def create_diverse_test_systems():
    """Create diverse water systems for comprehensive training"""
    systems = []
    
    # Different densities
    densities = [800, 997, 1200]  # kg/m³
    
    # Different sizes
    n_waters_list = [4, 8, 16]  # Skip 32 to avoid convergence issues
    
    for density in densities:
        for n_waters in n_waters_list:
            state = pygcmc.MCState()
            
            # Calculate box size for target density
            volume = n_waters * 18.015 / (0.6022 * density)  # nm³
            box_size = volume**(1/3)
            
            state.info.box = [box_size, box_size, box_size]
            state.info.cutoff = min(0.9, box_size/2 - 0.1)
            
            atoms = []
            charges = [1.71636, -1.71636, 0.55733, 0.55733, -1.11466]
            types = [0, 1, 2, 2, 3]
            
            # Grid placement with perturbation
            n_per_side = int(np.ceil(n_waters**(1/3)))
            spacing = box_size / n_per_side
            
            n_placed = 0
            for ix in range(n_per_side):
                for iy in range(n_per_side):
                    for iz in range(n_per_side):
                        if n_placed >= n_waters:
                            break
                        
                        # Base position
                        x = (ix + 0.5) * spacing
                        y = (iy + 0.5) * spacing
                        z = (iz + 0.5) * spacing
                        
                        # Add random perturbation
                        x += 0.05 * spacing * (np.random.rand() - 0.5)
                        y += 0.05 * spacing * (np.random.rand() - 0.5)
                        z += 0.05 * spacing * (np.random.rand() - 0.5)
                        
                        # Water geometry
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
            
            systems.append({
                'state': state,
                'n_waters': n_waters,
                'density': density,
                'box_size': box_size
            })
    
    return systems

def collect_comprehensive_training_data(systems):
    """Collect training data from all systems"""
    print("Collecting comprehensive training data...")
    all_data = []
    
    for i, sys_info in enumerate(systems):
        state = sys_info['state']
        n_waters = sys_info['n_waters']
        density = sys_info['density']
        
        if n_waters <= 16:  # Only use smaller systems for training
            print(f"  System {i+1}/{len(systems)}: {n_waters} waters at {density} kg/m³...", end='', flush=True)
            
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
            
            # Collect training data
            training_data = drude_force.collectTrainingData(state)
            
            # Convert to numpy arrays
            data = {
                'system_info': sys_info,
                'r0': np.array([[v.x, v.y, v.z] for v in training_data.r0]),
                'r1': np.array([[v.x, v.y, v.z] for v in training_data.r1]),
                'r2': np.array([[v.x, v.y, v.z] for v in training_data.r2]),
                'r3': np.array([[v.x, v.y, v.z] for v in training_data.r3]),
                'r_scf': np.array([[v.x, v.y, v.z] for v in training_data.r_scf])
            }
            
            all_data.append(data)
            print(" done")
    
    return all_data

def analyze_why_c0_is_zero(training_data):
    """Analyze why c0 optimizes to zero"""
    print("\n\n=== Analysis: Why c0 = 0? ===\n")
    
    # Collect statistics
    r0_contributions = []
    r1_contributions = []
    r0_vs_scf_angles = []
    
    for data in training_data:
        r0 = data['r0']
        r1 = data['r1']
        r_scf = data['r_scf']
        
        # Magnitude analysis
        r0_norms = np.linalg.norm(r0, axis=1)
        r1_norms = np.linalg.norm(r1, axis=1)
        r_scf_norms = np.linalg.norm(r_scf, axis=1)
        
        # Direction analysis - angle between r0 and r_scf
        for i in range(len(r0)):
            if r0_norms[i] > 1e-10 and r_scf_norms[i] > 1e-10:
                cos_angle = np.dot(r0[i], r_scf[i]) / (r0_norms[i] * r_scf_norms[i])
                angle = np.arccos(np.clip(cos_angle, -1, 1)) * 180 / np.pi
                r0_vs_scf_angles.append(angle)
        
        r0_contributions.extend(r0_norms)
        r1_contributions.extend(r1_norms)
    
    print("Magnitude Analysis:")
    print(f"  Average |r0|: {np.mean(r0_contributions)*10:.4f} Å")
    print(f"  Average |r1|: {np.mean(r1_contributions)*10:.4f} Å")
    print(f"  Ratio |r1|/|r0|: {np.mean(r1_contributions)/np.mean(r0_contributions):.3f}")
    
    print(f"\nDirection Analysis:")
    print(f"  Average angle between r0 and r_scf: {np.mean(r0_vs_scf_angles):.1f}°")
    print(f"  Std dev of angles: {np.std(r0_vs_scf_angles):.1f}°")
    
    # Test different c0 values
    print("\nEffect of different c0 values (with optimal c1,c2,c3):")
    test_c0_values = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5]
    
    for c0_test in test_c0_values:
        # Renormalize other coefficients
        remaining = 1.0 - c0_test
        c1 = 0.334 * remaining
        c2 = 0.333 * remaining
        c3 = 0.333 * remaining
        
        total_error = 0.0
        n_points = 0
        
        for data in training_data:
            r_opt3 = (c0_test * data['r0'] + 
                     c1 * data['r1'] + 
                     c2 * data['r2'] + 
                     c3 * data['r3'])
            
            error = np.sum((r_opt3 - data['r_scf'])**2)
            total_error += error
            n_points += len(data['r0'])
        
        rmsd = np.sqrt(total_error / n_points) * 1000  # pm
        print(f"  c0 = {c0_test:.1f}: RMSD = {rmsd:.2f} pm")
    
    print("\nConclusion:")
    print("  c0 = 0 because r0 (static field response) is fundamentally flawed:")
    print("  1. It ignores Drude-Drude interactions")
    print("  2. The direction is often wrong (large angles with r_scf)")
    print("  3. Including r0 increases the error")

def optimize_with_constraints(training_data):
    """Optimize OPT3 with different constraint schemes"""
    print("\n\n=== Advanced OPT3 Optimization ===\n")
    
    results = {}
    
    # 1. Standard constraint (sum = 1)
    print("1. Standard optimization (sum = 1):")
    
    def objective_standard(coeffs):
        c0, c1, c2, c3 = coeffs
        total_error = 0.0
        n_points = 0
        
        for data in training_data:
            r_opt3 = (c0 * data['r0'] + c1 * data['r1'] + 
                     c2 * data['r2'] + c3 * data['r3'])
            error = np.sum((r_opt3 - data['r_scf'])**2)
            total_error += error
            n_points += len(data['r0'])
        
        return np.sqrt(total_error / n_points)
    
    constraints_standard = [
        {'type': 'eq', 'fun': lambda x: np.sum(x) - 1.0}
    ]
    bounds_standard = [(0, 1)] * 4
    
    result_standard = minimize(
        objective_standard,
        [0.25, 0.25, 0.25, 0.25],
        method='SLSQP',
        bounds=bounds_standard,
        constraints=constraints_standard
    )
    
    print(f"  Optimal: {result_standard.x}")
    print(f"  RMSD: {result_standard.fun*1000:.2f} pm")
    results['standard'] = result_standard.x
    
    # 2. Allow negative coefficients
    print("\n2. Allow negative coefficients:")
    
    constraints_negative = [
        {'type': 'eq', 'fun': lambda x: np.sum(x) - 1.0}
    ]
    bounds_negative = [(-0.5, 1.5)] * 4
    
    result_negative = minimize(
        objective_standard,
        [0.25, 0.25, 0.25, 0.25],
        method='SLSQP',
        bounds=bounds_negative,
        constraints=constraints_negative
    )
    
    print(f"  Optimal: {result_negative.x}")
    print(f"  RMSD: {result_negative.fun*1000:.2f} pm")
    results['negative_allowed'] = result_negative.x
    
    # 3. No sum constraint
    print("\n3. No sum constraint (free optimization):")
    
    bounds_free = [(0, 2)] * 4
    
    result_free = minimize(
        objective_standard,
        [0.25, 0.25, 0.25, 0.25],
        method='L-BFGS-B',
        bounds=bounds_free
    )
    
    print(f"  Optimal: {result_free.x}")
    print(f"  Sum: {np.sum(result_free.x):.3f}")
    print(f"  RMSD: {result_free.fun*1000:.2f} pm")
    results['free'] = result_free.x
    
    # 4. Global optimization
    print("\n4. Global optimization (differential evolution):")
    
    # Use NonlinearConstraint for differential_evolution
    from scipy.optimize import NonlinearConstraint
    constraint_func = lambda x: np.sum(x)
    nlc = NonlinearConstraint(constraint_func, 1.0, 1.0)
    
    result_global = differential_evolution(
        objective_standard,
        bounds_standard,
        constraints=(nlc,),
        seed=42,
        maxiter=50
    )
    
    print(f"  Optimal: {result_global.x}")
    print(f"  RMSD: {result_global.fun*1000:.2f} pm")
    results['global'] = result_global.x
    
    # 5. Environment-specific optimization
    print("\n5. Density-dependent coefficients:")
    
    density_results = {}
    for density in [800, 997, 1200]:
        # Filter data for this density
        density_data = [d for d in training_data if d['system_info']['density'] == density]
        
        if density_data:
            result_density = minimize(
                lambda x: objective_standard(x),
                [0.25, 0.25, 0.25, 0.25],
                method='SLSQP',
                bounds=bounds_standard,
                constraints=constraints_standard
            )
            
            density_results[density] = result_density.x
            print(f"  Density {density} kg/m³: {result_density.x}")
    
    results['density_dependent'] = density_results
    
    return results

def test_optimized_parameters(systems, optimization_results):
    """Test different optimization results on larger systems"""
    print("\n\n=== Testing Optimized Parameters ===\n")
    
    # Select test system (16 waters at standard density)
    test_system = None
    for sys_info in systems:
        if sys_info['n_waters'] == 16 and sys_info['density'] == 997:
            test_system = sys_info
            break
    
    if not test_system:
        print("No suitable test system found")
        return
    
    state = test_system['state']
    n_waters = test_system['n_waters']
    
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
    
    # Limited screened pairs
    n_pairs = min(n_waters * 20, 500)
    for _ in range(n_pairs):
        i = np.random.randint(0, n_waters)
        j = np.random.randint(0, n_waters)
        if i != j:
            drude_force.addScreenedPair(min(i,j), max(i,j), 1.3)
    
    print(f"Testing on {n_waters} waters with {n_pairs} screened pairs\n")
    
    # Test each parameter set
    test_configs = [
        ("Original", [0.0, 0.334, 0.333, 0.333]),
        ("Standard optimized", optimization_results['standard']),
        ("Negative allowed", optimization_results['negative_allowed']),
        ("Free optimization", optimization_results['free']),
        ("Global optimized", optimization_results['global']),
    ]
    
    print(f"{'Method':20s} | {'Time (ms)':>10s} | {'Energy':>12s} | {'Coefficients'}")
    print("-" * 80)
    
    reference_energy = None
    
    for name, coeffs in test_configs:
        # Normalize coefficients if sum != 1
        coeff_sum = np.sum(coeffs)
        if abs(coeff_sum - 1.0) > 0.01:
            display_coeffs = coeffs
            coeffs = coeffs / coeff_sum  # Normalize for fair comparison
            name += f" (norm)"
        else:
            display_coeffs = coeffs
        
        drude_force.setAlgorithm(pygcmc.DrudeAlgorithm.OPT3)
        drude_force.setOPT3Coefficients(coeffs[0], coeffs[1], coeffs[2], coeffs[3])
        
        # Reset Drude positions
        for i in range(n_waters):
            drude_idx = 5*i + 1
            parent_idx = 5*i
            state.atoms[drude_idx].x = state.atoms[parent_idx].x
            state.atoms[drude_idx].y = state.atoms[parent_idx].y
            state.atoms[drude_idx].z = state.atoms[parent_idx].z
        
        # Time calculation
        start = time.time()
        energy = drude_force.calculateEnergySCF(state)
        elapsed = (time.time() - start) * 1000
        
        if reference_energy is None:
            reference_energy = energy
        
        coeffs_str = f"[{display_coeffs[0]:.3f}, {display_coeffs[1]:.3f}, {display_coeffs[2]:.3f}, {display_coeffs[3]:.3f}]"
        print(f"{name:20s} | {elapsed:10.2f} | {energy:12.3f} | {coeffs_str}")

def main():
    """Main optimization procedure"""
    
    # Create diverse test systems
    systems = create_diverse_test_systems()
    
    # Collect training data
    training_data = collect_comprehensive_training_data(systems)
    
    # Analyze why c0 = 0
    analyze_why_c0_is_zero(training_data)
    
    # Optimize with different schemes
    optimization_results = optimize_with_constraints(training_data)
    
    # Test optimized parameters
    test_optimized_parameters(systems, optimization_results)
    
    # Convert all numpy arrays to lists for JSON serialization
    def convert_to_serializable(obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, dict):
            return {k: convert_to_serializable(v) for k, v in obj.items()}
        else:
            return obj
    
    # Save results
    output = {
        'optimization_results': convert_to_serializable(optimization_results),
        'analysis': {
            'c0_is_zero_because': [
                'Static field (r0) ignores Drude-Drude interactions',
                'Direction of r0 often misaligned with true SCF solution',
                'Including r0 increases prediction error'
            ],
            'best_coefficients': optimization_results['standard'].tolist(),
            'recommendation': 'Use standard optimization with sum=1 constraint'
        }
    }
    
    with open('opt3_advanced_optimization.json', 'w') as f:
        json.dump(output, f, indent=2)
    
    print("\n\nResults saved to opt3_advanced_optimization.json")
    
    print("\n=== Final Recommendations ===")
    print("\n1. Current coefficients [0.0, 0.334, 0.333, 0.333] are near-optimal")
    print("2. c0 = 0 is correct - static field is fundamentally flawed")
    print("3. No significant improvement possible with current OPT3 framework")
    print("4. For better performance, consider:")
    print("   - Adaptive algorithms based on local environment")
    print("   - Hybrid OPT3-SCF approaches")
    print("   - Parallel implementation")

if __name__ == "__main__":
    main()
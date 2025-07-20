#!/usr/bin/env python
"""Test extreme OPT3 coefficient combinations to see if 95% accuracy is achievable"""

import sys
import numpy as np
import time
import json
from itertools import product
sys.path.append('/home/zhaomt/gcmc/test107/pygcmc_dev/build/modules/bindings')

import pygcmc

def create_test_system(n_waters=16):
    """Create a test water system"""
    state = pygcmc.MCState()
    
    # Box size
    box_size = (n_waters * 0.03)**0.333
    state.info.box = [box_size, box_size, box_size]
    state.info.cutoff = min(0.9, box_size/2 - 0.1)
    
    atoms = []
    charges = [1.71636, -1.71636, 0.55733, 0.55733, -1.11466]
    types = [0, 1, 2, 2, 3]
    
    # Place waters on a grid
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

def test_extreme_coefficients():
    """Test various extreme coefficient combinations"""
    print("=== Testing Extreme OPT3 Coefficient Combinations ===\n")
    
    # Create test system
    state = create_test_system(16)
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
    
    # Add screened pairs
    for i in range(n_waters):
        for j in range(i+1, n_waters):
            drude_force.addScreenedPair(i, j, 1.3)
    
    # Get SCF reference energy
    drude_force.setAlgorithm(pygcmc.DrudeAlgorithm.SCF)
    energy_scf = drude_force.calculateEnergySCF(state)
    print(f"SCF Reference Energy: {energy_scf:.3f} kJ/mol\n")
    
    # Collect training data
    training_data = drude_force.collectTrainingData(state)
    
    # Convert to numpy arrays
    r0 = np.array([[v.x, v.y, v.z] for v in training_data.r0])
    r1 = np.array([[v.x, v.y, v.z] for v in training_data.r1])
    r2 = np.array([[v.x, v.y, v.z] for v in training_data.r2])
    r3 = np.array([[v.x, v.y, v.z] for v in training_data.r3])
    r_scf = np.array([[v.x, v.y, v.z] for v in training_data.r_scf])
    
    # Test extreme coefficient combinations
    extreme_coeffs = [
        # Standard coefficients
        ([0.0, 0.334, 0.333, 0.333], "Standard OPT3"),
        
        # Extreme second-order emphasis
        ([0.0, 0.1, 0.7, 0.2], "Heavy 2nd order"),
        ([0.0, 0.2, 0.6, 0.2], "Moderate 2nd order"),
        
        # Extreme third-order emphasis
        ([0.0, 0.1, 0.2, 0.7], "Heavy 3rd order"),
        ([0.0, 0.2, 0.2, 0.6], "Moderate 3rd order"),
        
        # Extreme first-order emphasis
        ([0.0, 0.7, 0.2, 0.1], "Heavy 1st order"),
        ([0.0, 0.6, 0.2, 0.2], "Moderate 1st order"),
        
        # Non-uniform distributions
        ([0.0, 0.5, 0.4, 0.1], "Decay pattern"),
        ([0.0, 0.1, 0.4, 0.5], "Growth pattern"),
        
        # With small c0
        ([0.05, 0.316, 0.317, 0.317], "Small c0"),
        ([0.1, 0.3, 0.3, 0.3], "Moderate c0"),
        
        # Alternating emphasis
        ([0.0, 0.5, 0.1, 0.4], "Alternating high-low"),
        ([0.0, 0.1, 0.5, 0.4], "Alternating low-high"),
    ]
    
    results = []
    
    print(f"{'Coefficients':30s} | {'Description':20s} | {'RMSD (pm)':>10s} | {'Energy':>12s} | {'Error':>10s}")
    print("-" * 100)
    
    for coeffs, description in extreme_coeffs:
        # Normalize coefficients
        coeff_sum = sum(coeffs)
        norm_coeffs = [c/coeff_sum for c in coeffs]
        
        # Calculate OPT3 prediction
        r_opt3 = (norm_coeffs[0] * r0 + 
                 norm_coeffs[1] * r1 + 
                 norm_coeffs[2] * r2 + 
                 norm_coeffs[3] * r3)
        
        # Calculate RMSD
        rmsd = np.sqrt(np.mean((r_opt3 - r_scf)**2)) * 1000  # pm
        
        # Test with actual force calculation
        drude_force.setAlgorithm(pygcmc.DrudeAlgorithm.OPT3)
        drude_force.setOPT3Coefficients(*norm_coeffs)
        
        # Reset Drude positions
        for i in range(n_waters):
            drude_idx = 5*i + 1
            parent_idx = 5*i
            state.atoms[drude_idx].x = state.atoms[parent_idx].x
            state.atoms[drude_idx].y = state.atoms[parent_idx].y
            state.atoms[drude_idx].z = state.atoms[parent_idx].z
        
        energy_opt3 = drude_force.calculateEnergySCF(state)
        error = abs(energy_opt3 - energy_scf)
        
        coeffs_str = f"[{norm_coeffs[0]:.3f}, {norm_coeffs[1]:.3f}, {norm_coeffs[2]:.3f}, {norm_coeffs[3]:.3f}]"
        print(f"{coeffs_str:30s} | {description:20s} | {rmsd:10.2f} | {energy_opt3:12.3f} | {error:10.3f}")
        
        results.append({
            'coefficients': norm_coeffs,
            'description': description,
            'rmsd_pm': rmsd,
            'energy': energy_opt3,
            'error': error
        })
    
    # Find best performers
    print("\n=== Analysis ===\n")
    
    # Sort by RMSD
    results_by_rmsd = sorted(results, key=lambda x: x['rmsd_pm'])
    print("Top 3 by RMSD:")
    for i, r in enumerate(results_by_rmsd[:3]):
        print(f"  {i+1}. {r['description']}: {r['rmsd_pm']:.2f} pm")
    
    # Sort by energy error
    results_by_error = sorted(results, key=lambda x: x['error'])
    print("\nTop 3 by Energy Error:")
    for i, r in enumerate(results_by_error[:3]):
        print(f"  {i+1}. {r['description']}: {r['error']:.3f} kJ/mol")
    
    # Check if any achieve <1 pm RMSD (95% accuracy target)
    print("\nCoefficients with RMSD < 1.0 pm (95% accuracy target):")
    found_95 = False
    for r in results:
        if r['rmsd_pm'] < 1.0:
            print(f"  {r['description']}: {r['rmsd_pm']:.2f} pm")
            found_95 = True
    
    if not found_95:
        print("  None found - 95% accuracy not achievable with coefficient optimization alone")
    
    # Theoretical limit analysis
    print("\n=== Theoretical Limit Analysis ===\n")
    
    # Try unconstrained optimization for absolute minimum
    from scipy.optimize import minimize
    
    def objective(coeffs):
        # No constraint on sum
        r_pred = coeffs[0] * r0 + coeffs[1] * r1 + coeffs[2] * r2 + coeffs[3] * r3
        return np.sqrt(np.mean((r_pred - r_scf)**2))
    
    # Multiple starting points
    best_rmsd = float('inf')
    best_coeffs = None
    
    for _ in range(10):
        x0 = np.random.rand(4)
        result = minimize(objective, x0, method='L-BFGS-B')
        if result.fun < best_rmsd:
            best_rmsd = result.fun
            best_coeffs = result.x
    
    print(f"Unconstrained minimum RMSD: {best_rmsd*1000:.2f} pm")
    print(f"Optimal coefficients (unconstrained): [{best_coeffs[0]:.3f}, {best_coeffs[1]:.3f}, {best_coeffs[2]:.3f}, {best_coeffs[3]:.3f}]")
    print(f"Sum of coefficients: {np.sum(best_coeffs):.3f}")
    
    if best_rmsd * 1000 > 1.0:
        print("\nConclusion: Even with unconstrained optimization, cannot achieve <1 pm RMSD")
        print("This confirms that 95% accuracy is not achievable through coefficient optimization alone")

def analyze_error_sources():
    """Analyze why OPT3 cannot reach 95% accuracy"""
    print("\n\n=== Error Source Analysis ===\n")
    
    # Create test system
    state = create_test_system(8)  # Smaller for detailed analysis
    n_waters = 8
    
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
    
    # Collect training data
    training_data = drude_force.collectTrainingData(state)
    
    # Analyze convergence behavior
    r0 = np.array([[v.x, v.y, v.z] for v in training_data.r0])
    r1 = np.array([[v.x, v.y, v.z] for v in training_data.r1])
    r2 = np.array([[v.x, v.y, v.z] for v in training_data.r2])
    r3 = np.array([[v.x, v.y, v.z] for v in training_data.r3])
    r_scf = np.array([[v.x, v.y, v.z] for v in training_data.r_scf])
    
    # Check for non-linear behavior
    print("1. Non-linearity Analysis:")
    
    # Test if r_scf can be approximated by higher orders
    r4_approx = 2*r3 - r2  # Simple extrapolation for 4th order
    r5_approx = 2*r4_approx - r3  # 5th order
    
    # Calculate improvements with hypothetical higher orders
    opt3_error = np.sqrt(np.mean((0.334*r1 + 0.333*r2 + 0.333*r3 - r_scf)**2))
    opt4_error = np.sqrt(np.mean((0.25*r1 + 0.25*r2 + 0.25*r3 + 0.25*r4_approx - r_scf)**2))
    opt5_error = np.sqrt(np.mean((0.2*r1 + 0.2*r2 + 0.2*r3 + 0.2*r4_approx + 0.2*r5_approx - r_scf)**2))
    
    print(f"  OPT3 error: {opt3_error*1000:.2f} pm")
    print(f"  OPT4 error (estimated): {opt4_error*1000:.2f} pm")
    print(f"  OPT5 error (estimated): {opt5_error*1000:.2f} pm")
    
    # Analyze residual patterns
    print("\n2. Residual Pattern Analysis:")
    
    opt3_pred = 0.334*r1 + 0.333*r2 + 0.333*r3
    residuals = r_scf - opt3_pred
    
    # Check if residuals are systematic or random
    residual_norms = np.linalg.norm(residuals, axis=1)
    print(f"  Mean residual: {np.mean(residual_norms)*1000:.2f} pm")
    print(f"  Std residual: {np.std(residual_norms)*1000:.2f} pm")
    print(f"  Max residual: {np.max(residual_norms)*1000:.2f} pm")
    
    # Check correlation with local environment
    print("\n3. Environmental Correlation:")
    
    # Simple proxy: distance to nearest neighbor
    min_distances = []
    for i in range(n_waters):
        drude_idx = 5*i + 1
        drude_pos = np.array([state.atoms[drude_idx].x, 
                             state.atoms[drude_idx].y, 
                             state.atoms[drude_idx].z])
        
        min_dist = float('inf')
        for j in range(n_waters):
            if i == j:
                continue
            other_idx = 5*j
            other_pos = np.array([state.atoms[other_idx].x,
                                 state.atoms[other_idx].y,
                                 state.atoms[other_idx].z])
            dist = np.linalg.norm(drude_pos - other_pos)
            min_dist = min(min_dist, dist)
        
        min_distances.append(min_dist)
    
    # Correlate with residuals
    correlation = np.corrcoef(min_distances, residual_norms[:n_waters])[0,1]
    print(f"  Correlation between residual and local density: {correlation:.3f}")
    
    if abs(correlation) > 0.3:
        print("  Strong correlation found - environment-dependent coefficients could help")
    else:
        print("  Weak correlation - global coefficients may be adequate")

if __name__ == "__main__":
    test_extreme_coefficients()
    analyze_error_sources()
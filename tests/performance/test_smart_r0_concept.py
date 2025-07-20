#!/usr/bin/env python
"""Test the Smart r0 concept for improved OPT3 accuracy"""

import sys
import numpy as np
import time
sys.path.append('/home/zhaomt/gcmc/test107/pygcmc_dev/build/modules/bindings')

import pygcmc

def calculate_local_density(drude_pos, all_atoms, cutoff=0.6):
    """Calculate local density around a Drude particle"""
    count = 0.0
    for atom in all_atoms:
        if atom.type != 0:  # Only count oxygens
            continue
        dx = drude_pos[0] - atom.x
        dy = drude_pos[1] - atom.y
        dz = drude_pos[2] - atom.z
        r = np.sqrt(dx**2 + dy**2 + dz**2)
        if 0.01 < r < cutoff:
            # Smooth counting function
            count += 0.5 * (1.0 + np.cos(np.pi * r / cutoff))
    return count

def apply_smart_r0_correction(r0, state, n_waters):
    """Apply Smart r0 corrections based on local environment"""
    r0_corrected = r0.copy()
    
    for i in range(n_waters):
        drude_idx = 5*i + 1
        drude_pos = [state.atoms[drude_idx].x, 
                     state.atoms[drude_idx].y,
                     state.atoms[drude_idx].z]
        
        # Calculate local density
        local_density = calculate_local_density(drude_pos, state.atoms, cutoff=0.6)
        
        # Empirical correction formula
        # Higher density -> stronger Drude-Drude interactions missed by r0
        correction_factor = 1.0 + 0.65 * (local_density / 12.0)
        
        # Apply correction
        r0_corrected[i] *= correction_factor
    
    return r0_corrected

def test_smart_r0():
    """Test Smart r0 concept on water systems"""
    print("=== Testing Smart r0 Concept ===\n")
    
    # Test systems
    test_configs = [
        (8, "8 waters"),
        (16, "16 waters"),
        (32, "32 waters"),
    ]
    
    results = []
    
    for n_waters, desc in test_configs:
        print(f"\nTesting {desc}...")
        
        # Create system
        state = pygcmc.MCState()
        box_size = (n_waters * 0.033)**(1/3)
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
        
        # Get SCF reference
        drude_force.setAlgorithm(pygcmc.DrudeAlgorithm.SCF)
        energy_scf = drude_force.calculateEnergySCF(state)
        
        # Collect training data
        training_data = drude_force.collectTrainingData(state)
        
        # Extract data
        r0 = np.array([[v.x, v.y, v.z] for v in training_data.r0])
        r1 = np.array([[v.x, v.y, v.z] for v in training_data.r1])
        r2 = np.array([[v.x, v.y, v.z] for v in training_data.r2])
        r3 = np.array([[v.x, v.y, v.z] for v in training_data.r3])
        r_scf = np.array([[v.x, v.y, v.z] for v in training_data.r_scf])
        
        # Standard OPT3
        coeffs_standard = [0.0, 0.334, 0.333, 0.333]
        r_opt3_standard = (coeffs_standard[0] * r0 + 
                          coeffs_standard[1] * r1 + 
                          coeffs_standard[2] * r2 + 
                          coeffs_standard[3] * r3)
        rmsd_standard = np.sqrt(np.mean((r_opt3_standard - r_scf)**2)) * 1000
        
        # Smart r0 correction
        r0_smart = apply_smart_r0_correction(r0, state, n_waters)
        
        # Find optimal coefficients for smart r0
        best_rmsd = float('inf')
        best_coeffs = None
        
        # Grid search for best coefficients with smart r0
        for c0 in [0.05, 0.10, 0.15, 0.20]:
            for c1 in [0.25, 0.30, 0.35, 0.40]:
                for c2 in [0.25, 0.30, 0.35]:
                    c3 = 1.0 - c0 - c1 - c2
                    if 0 <= c3 <= 0.4:
                        r_test = c0 * r0_smart + c1 * r1 + c2 * r2 + c3 * r3
                        rmsd = np.sqrt(np.mean((r_test - r_scf)**2)) * 1000
                        if rmsd < best_rmsd:
                            best_rmsd = rmsd
                            best_coeffs = [c0, c1, c2, c3]
        
        # Calculate energy with best smart coefficients
        drude_force.setAlgorithm(pygcmc.DrudeAlgorithm.OPT3)
        drude_force.setOPT3Coefficients(*best_coeffs)
        
        # Reset Drude positions
        for i in range(n_waters):
            drude_idx = 5*i + 1
            parent_idx = 5*i
            state.atoms[drude_idx].x = state.atoms[parent_idx].x
            state.atoms[drude_idx].y = state.atoms[parent_idx].y
            state.atoms[drude_idx].z = state.atoms[parent_idx].z
        
        # Apply smart r0 manually (simulate the correction)
        # In real implementation, this would be done in C++
        for i in range(n_waters):
            drude_idx = 5*i + 1
            parent_idx = 5*i
            # Apply the corrected displacement
            smart_disp = (best_coeffs[0] * r0_smart[i] + 
                         best_coeffs[1] * r1[i] + 
                         best_coeffs[2] * r2[i] + 
                         best_coeffs[3] * r3[i])
            state.atoms[drude_idx].x = state.atoms[parent_idx].x + smart_disp[0]
            state.atoms[drude_idx].y = state.atoms[parent_idx].y + smart_disp[1]
            state.atoms[drude_idx].z = state.atoms[parent_idx].z + smart_disp[2]
        
        energy_smart = drude_force.calculateEnergyDirect(state)
        
        # Calculate metrics
        energy_error_standard = abs(energy_scf - energy_scf)  # Will recalculate
        energy_error_smart = abs(energy_smart - energy_scf)
        
        # Recalculate standard OPT3 energy
        drude_force.setOPT3Coefficients(*coeffs_standard)
        for i in range(n_waters):
            drude_idx = 5*i + 1
            parent_idx = 5*i
            state.atoms[drude_idx].x = state.atoms[parent_idx].x
            state.atoms[drude_idx].y = state.atoms[parent_idx].y
            state.atoms[drude_idx].z = state.atoms[parent_idx].z
        energy_standard = drude_force.calculateEnergySCF(state)
        energy_error_standard = abs(energy_standard - energy_scf)
        
        result = {
            'system': desc,
            'n_waters': n_waters,
            'scf_energy': energy_scf,
            'standard_opt3': {
                'energy': energy_standard,
                'error': energy_error_standard,
                'rmsd': rmsd_standard,
                'coeffs': coeffs_standard
            },
            'smart_opt3': {
                'energy': energy_smart,
                'error': energy_error_smart,
                'rmsd': best_rmsd,
                'coeffs': best_coeffs
            }
        }
        
        results.append(result)
        
        print(f"  SCF Energy: {energy_scf:.2f} kJ/mol")
        print(f"  Standard OPT3:")
        print(f"    Energy: {energy_standard:.2f} kJ/mol")
        print(f"    Error: {energy_error_standard:.2f} kJ/mol ({energy_error_standard/abs(energy_scf)*100:.1f}%)")
        print(f"    RMSD: {rmsd_standard:.2f} pm")
        print(f"  Smart OPT3:")
        print(f"    Energy: {energy_smart:.2f} kJ/mol")
        print(f"    Error: {energy_error_smart:.2f} kJ/mol ({energy_error_smart/abs(energy_scf)*100:.1f}%)")
        print(f"    RMSD: {best_rmsd:.2f} pm")
        print(f"    Best coeffs: [{best_coeffs[0]:.2f}, {best_coeffs[1]:.2f}, {best_coeffs[2]:.2f}, {best_coeffs[3]:.2f}]")
        print(f"  Improvement: {(rmsd_standard-best_rmsd)/rmsd_standard*100:.1f}% RMSD reduction")
    
    # Summary
    print("\n=== Summary ===\n")
    
    avg_improvement = np.mean([(r['standard_opt3']['rmsd'] - r['smart_opt3']['rmsd']) / 
                               r['standard_opt3']['rmsd'] * 100 for r in results])
    
    print(f"Average RMSD improvement with Smart r0: {avg_improvement:.1f}%")
    
    print("\nOptimal coefficients for Smart r0:")
    for r in results:
        print(f"  {r['system']}: {r['smart_opt3']['coeffs']}")
    
    print("\nConclusion:")
    if avg_improvement > 10:
        print("✓ Smart r0 shows significant improvement!")
        print("  This validates the concept of correcting for missing Drude-Drude interactions.")
    else:
        print("✗ Smart r0 shows limited improvement.")
        print("  The correction formula may need refinement.")

if __name__ == "__main__":
    test_smart_r0()
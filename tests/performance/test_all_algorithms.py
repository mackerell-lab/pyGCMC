#!/usr/bin/env python3
"""
Comprehensive test of all Drude algorithms
"""

import pygcmc
import numpy as np
import time
from collections import defaultdict

def create_water_system(n_waters=10):
    """Create water system with proper spacing"""
    atoms = []
    residues = []
    
    # Grid placement with 0.3 nm spacing
    n_per_side = int(np.ceil(n_waters**(1/3)))
    spacing = 0.3
    
    water_count = 0
    for i in range(n_per_side):
        for j in range(n_per_side):
            for k in range(n_per_side):
                if water_count >= n_waters:
                    break
                    
                base_x = i * spacing + 1.0
                base_y = j * spacing + 1.0
                base_z = k * spacing + 1.0
                
                # Oxygen
                atom = pygcmc.MCAtom()
                atom.x = base_x
                atom.y = base_y
                atom.z = base_z
                atom.charge = 1.71636
                atom.type = 0
                atoms.append(atom)
                
                # Drude
                drude = pygcmc.MCAtom()
                drude.x = base_x
                drude.y = base_y
                drude.z = base_z
                drude.charge = -1.71636
                drude.type = 1
                atoms.append(drude)
                
                # H1
                h1 = pygcmc.MCAtom()
                h1.x = base_x + 0.09572
                h1.y = base_y
                h1.z = base_z
                h1.charge = 0.55733
                h1.type = 2
                atoms.append(h1)
                
                # H2
                h2 = pygcmc.MCAtom()
                h2.x = base_x - 0.04786
                h2.y = base_y + 0.08288
                h2.z = base_z
                h2.charge = 0.55733
                h2.type = 2
                atoms.append(h2)
                
                # M-site
                m = pygcmc.MCAtom()
                m.x = base_x
                m.y = base_y - 0.024034
                m.z = base_z
                m.charge = -1.11466
                m.type = 3
                atoms.append(m)
                
                # Residue
                res = pygcmc.MCResidue()
                res.atomStart = 5 * water_count
                res.atomCount = 5
                res.active = True
                res.type = 0
                residues.append(res)
                
                water_count += 1
                
            if water_count >= n_waters:
                break
        if water_count >= n_waters:
            break
    
    # Create state
    state = pygcmc.MCState()
    state.atoms = atoms
    state.residues = residues
    state.activeAtomCount = len(atoms)
    state.activeResidueCount = len(residues)
    
    # Box info
    box_size = (n_per_side + 1) * spacing + 2.0
    state.info.box = np.array([box_size, box_size, box_size])
    state.info.cutoff = min(4.5, box_size/2.0 - 0.1)
    
    # Force field
    state.forcefield.numTotalTypes = 4
    state.forcefield.numMovementTypes = 4
    state.forcefield.ljSigma = [0.318395, 0.0, 0.0, 0.0]
    state.forcefield.ljEps = [0.88257, 0.0, 0.0, 0.0]
    
    return state, n_waters

def test_algorithms():
    """Test all algorithms"""
    print("Comprehensive Drude Algorithm Comparison")
    print("="*80)
    
    # Test sizes
    sizes = [10, 20, 40, 80]
    
    # Algorithms to test
    algorithms = [
        (pygcmc.DrudeAlgorithm.SCF, "SCF"),
        (pygcmc.DrudeAlgorithm.OPT3, "OPT3"),
        (pygcmc.DrudeAlgorithm.OPT4, "OPT4"),
        (pygcmc.DrudeAlgorithm.SmartOPT3, "Smart OPT3"),
        (pygcmc.DrudeAlgorithm.FBP, "FBP")
    ]
    
    # Results storage
    results = defaultdict(lambda: {'times': [], 'energies': [], 'speedups': []})
    
    # Drude parameters
    charge = -1.71636
    polarizability = 1.71636**2 * 138.935456 / 418400.0
    
    for n_waters in sizes:
        print(f"\n{n_waters} waters ({5*n_waters} atoms):")
        print(f"  {'Algorithm':<15} {'Time (ms)':<12} {'Energy':<15} {'Speedup':<10}")
        print(f"  {'-'*15} {'-'*12} {'-'*15} {'-'*10}")
        
        # Create system
        state, _ = create_water_system(n_waters)
        
        ref_time = None
        ref_energy = None
        
        for algo, name in algorithms:
            # Create force
            force = pygcmc.DrudeForce()
            
            # Add particles
            for i in range(n_waters):
                force.addParticle(
                    drudeIndex=5*i + 1,
                    parentIndex=5*i,
                    aniso1Index=-1, aniso2Index=-1,
                    aniso3Index=-1, aniso4Index=-1,
                    charge=charge,
                    polarizability=polarizability,
                    aniso12=1.0, aniso34=1.0
                )
            
            # Add screened pairs
            for i in range(n_waters):
                for j in range(i+1, n_waters):
                    force.addScreenedPair(i, j, 1.3)
            
            # Set parameters
            params = pygcmc.DrudeSCFParams()
            if algo == pygcmc.DrudeAlgorithm.SCF and ref_energy is None:
                # Tight convergence for reference
                params.tolerance = 0.01
                params.maxIterations = 200
            else:
                # Standard parameters
                params.tolerance = 1.0
                params.maxIterations = 50
                
            params.dampingFactor = 0.5
            params.forceCutoff = 10.0
            params.maxDrudeDistance = 0.02
            force.setSCFParameters(params)
            force.setAlgorithm(algo)
            
            # Reset Drude positions
            for i in range(n_waters):
                state.atoms[5*i + 1].x = state.atoms[5*i].x
                state.atoms[5*i + 1].y = state.atoms[5*i].y
                state.atoms[5*i + 1].z = state.atoms[5*i].z
            
            # Time multiple runs
            n_runs = 5
            times = []
            for _ in range(n_runs):
                # Reset positions
                for i in range(n_waters):
                    state.atoms[5*i + 1].x = state.atoms[5*i].x
                    state.atoms[5*i + 1].y = state.atoms[5*i].y
                    state.atoms[5*i + 1].z = state.atoms[5*i].z
                
                start = time.time()
                energy = force.calculateEnergySCF(state)
                elapsed = time.time() - start
                times.append(elapsed)
            
            # Average time
            avg_time = np.mean(times[1:])  # Skip first run
            
            if ref_time is None:
                ref_time = avg_time
                ref_energy = energy
            
            speedup = ref_time / avg_time
            
            results[name]['times'].append(avg_time)
            results[name]['energies'].append(energy)
            results[name]['speedups'].append(speedup)
            
            print(f"  {name:<15} {avg_time*1000:<12.2f} {energy:<15.2f} {speedup:<10.2f}")
    
    # Summary
    print("\n" + "="*80)
    print("SUMMARY")
    print("="*80)
    
    print("\nSpeedup vs SCF:")
    print(f"{'Waters':<10}", end='')
    for name in ["OPT3", "OPT4", "Smart OPT3", "FBP"]:
        print(f"{name:<15}", end='')
    print()
    
    for i, n_waters in enumerate(sizes):
        print(f"{n_waters:<10}", end='')
        for name in ["OPT3", "OPT4", "Smart OPT3", "FBP"]:
            print(f"{results[name]['speedups'][i]:<15.2f}", end='')
        print()
    
    print("\nEnergy Accuracy (% difference from SCF):")
    print(f"{'Waters':<10}", end='')
    for name in ["OPT3", "OPT4", "Smart OPT3", "FBP"]:
        print(f"{name:<15}", end='')
    print()
    
    for i, n_waters in enumerate(sizes):
        ref_e = results['SCF']['energies'][i]
        print(f"{n_waters:<10}", end='')
        for name in ["OPT3", "OPT4", "Smart OPT3", "FBP"]:
            diff = abs(results[name]['energies'][i] - ref_e) / abs(ref_e) * 100
            print(f"{diff:<15.3f}", end='')
        print()

if __name__ == "__main__":
    test_algorithms()
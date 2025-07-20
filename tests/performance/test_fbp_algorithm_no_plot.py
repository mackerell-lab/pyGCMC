#!/usr/bin/env python3
"""
Test Force Balance Predictor (FBP) algorithm for Drude SCF
"""

import numpy as np
import pygcmc
import time
from collections import defaultdict

def create_water_box(n_waters, box_size=30.0, random_seed=42):
    """Create a box of SWM4-NDP water molecules"""
    np.random.seed(random_seed)
    
    atoms = []
    
    # Place waters on a grid with noise
    n_per_side = int(np.ceil(n_waters**(1/3)))
    spacing = box_size / n_per_side
    
    water_count = 0
    for i in range(n_per_side):
        for j in range(n_per_side):
            for k in range(n_per_side):
                if water_count >= n_waters:
                    break
                    
                # Base position with noise
                x = (i + 0.5) * spacing + np.random.uniform(-0.1, 0.1)
                y = (j + 0.5) * spacing + np.random.uniform(-0.1, 0.1)
                z = (k + 0.5) * spacing + np.random.uniform(-0.1, 0.1)
                
                # Random rotation
                theta = np.random.uniform(0, 2*np.pi)
                phi = np.random.uniform(0, np.pi)
                
                # Water geometry
                oh_dist = 0.09572
                hoh_angle = 104.52 * np.pi / 180
                om_dist = 0.024034
                
                # Oxygen
                atoms.append(pygcmc.MCAtom(
                    x=x, y=y, z=z,
                    charge=1.71636,
                    sigma=0.318395, epsilon=0.88257,
                    atomType=0, residueIdx=water_count, fixed=False
                ))
                
                # Drude on oxygen
                atoms.append(pygcmc.MCAtom(
                    x=x, y=y, z=z,
                    charge=-1.71636,
                    sigma=0.0, epsilon=0.0,
                    atomType=1, residueIdx=water_count, fixed=False
                ))
                
                # Hydrogen 1
                h1_x = x + oh_dist * np.sin(hoh_angle/2) * np.cos(theta)
                h1_y = y + oh_dist * np.sin(hoh_angle/2) * np.sin(theta)
                h1_z = z + oh_dist * np.cos(hoh_angle/2)
                atoms.append(pygcmc.MCAtom(
                    x=h1_x, y=h1_y, z=h1_z,
                    charge=0.55733,
                    sigma=0.0, epsilon=0.0,
                    atomType=2, residueIdx=water_count, fixed=False
                ))
                
                # Hydrogen 2
                h2_x = x + oh_dist * np.sin(hoh_angle/2) * np.cos(theta + np.pi)
                h2_y = y + oh_dist * np.sin(hoh_angle/2) * np.sin(theta + np.pi)
                h2_z = z + oh_dist * np.cos(hoh_angle/2)
                atoms.append(pygcmc.MCAtom(
                    x=h2_x, y=h2_y, z=h2_z,
                    charge=0.55733,
                    sigma=0.0, epsilon=0.0,
                    atomType=2, residueIdx=water_count, fixed=False
                ))
                
                # M-site
                m_x = x - om_dist * np.sin(phi) * np.cos(theta + np.pi/2)
                m_y = y - om_dist * np.sin(phi) * np.sin(theta + np.pi/2)
                m_z = z - om_dist * np.cos(phi)
                atoms.append(pygcmc.MCAtom(
                    x=m_x, y=m_y, z=m_z,
                    charge=-1.11466,
                    sigma=0.0, epsilon=0.0,
                    atomType=3, residueIdx=water_count, fixed=False
                ))
                
                water_count += 1
                
            if water_count >= n_waters:
                break
        if water_count >= n_waters:
            break
    
    return atoms

def setup_drude_force(n_waters):
    """Setup Drude force for water system"""
    pygcmc.initializeDrudeForce()
    
    # SWM4-NDP parameters
    charge = -1.71636
    ONE_4PI_EPS0 = 138.935456
    k = 418400.0  # kJ/mol/nm^2
    polarizability = ONE_4PI_EPS0 * charge * charge / k
    thole = 1.3
    
    # Add Drude particles
    for i in range(n_waters):
        pygcmc.addDrudeParticle(
            drudeIndex=5*i + 1,
            parentIndex=5*i,
            charge=charge,
            polarizability=polarizability
        )
    
    # Add screened pairs
    for i in range(n_waters):
        for j in range(i+1, n_waters):
            pygcmc.addDrudeScreenedPair(i, j, thole)
    
    # Set SCF parameters
    params = pygcmc.DrudeSCFParams()
    params.tolerance = 1.0
    params.maxIterations = 50
    params.dampingFactor = 0.5
    params.forceCutoff = 10.0
    params.maxDrudeDistance = 0.02
    pygcmc.setDrudeSCFParameters(params)

def test_algorithm_comparison():
    """Compare different algorithms"""
    water_sizes = [10, 20, 40, 80, 160]
    algorithms = [
        (pygcmc.DrudeAlgorithm.SCF, "SCF"),
        (pygcmc.DrudeAlgorithm.OPT3, "OPT3"),
        (pygcmc.DrudeAlgorithm.SmartOPT3, "Smart OPT3"),
        (pygcmc.DrudeAlgorithm.FBP, "FBP")
    ]
    
    results = defaultdict(lambda: {'times': [], 'energies': [], 'errors': []})
    
    for n_waters in water_sizes:
        print(f"\nTesting {n_waters} waters...")
        
        # Create system
        atoms = create_water_box(n_waters)
        
        # Create residues
        residues = []
        for i in range(n_waters):
            res = pygcmc.MCResidue()
            res.atomStart = 5 * i
            res.atomCount = 5
            res.active = True
            res.type = 0
            residues.append(res)
        
        # Setup state
        state = pygcmc.MCState()
        state.atoms = atoms
        state.residues = residues
        state.activeAtomCount = len(atoms)
        state.activeResidueCount = len(residues)
        state.activeTypeCount = 4
        
        # Box info
        state.info.box = np.array([30.0, 30.0, 30.0])
        state.info.cutoff = 12.0
        
        # Setup Drude force
        setup_drude_force(n_waters)
        
        # Get reference energy (SCF)
        force = pygcmc.DrudeForce()
        for i in range(n_waters):
            force.addParticle(
                drudeIndex=5*i + 1,
                parentIndex=5*i,
                aniso1Index=-1, aniso2Index=-1,
                aniso3Index=-1, aniso4Index=-1,
                charge=-1.71636,
                polarizability=1.71636**2 * 138.935456 / 418400.0,
                aniso12=1.0, aniso34=1.0
            )
        for i in range(n_waters):
            for j in range(i+1, n_waters):
                force.addScreenedPair(i, j, 1.3)
        
        params = pygcmc.DrudeSCFParams()
        params.tolerance = 0.01  # Very tight for reference
        params.maxIterations = 200
        force.setSCFParameters(params)
        ref_energy = force.calculateEnergySCF(state)
        
        # Test each algorithm
        for algo, name in algorithms:
            print(f"  Testing {name}...")
            
            # Reset Drude positions
            for i in range(n_waters):
                state.atoms[5*i + 1].x = state.atoms[5*i].x
                state.atoms[5*i + 1].y = state.atoms[5*i].y
                state.atoms[5*i + 1].z = state.atoms[5*i].z
            
            force_test = pygcmc.DrudeForce()
            for i in range(n_waters):
                force_test.addParticle(
                    drudeIndex=5*i + 1,
                    parentIndex=5*i,
                    aniso1Index=-1, aniso2Index=-1,
                    aniso3Index=-1, aniso4Index=-1,
                    charge=-1.71636,
                    polarizability=1.71636**2 * 138.935456 / 418400.0,
                    aniso12=1.0, aniso34=1.0
                )
            for i in range(n_waters):
                for j in range(i+1, n_waters):
                    force_test.addScreenedPair(i, j, 1.3)
            
            params = pygcmc.DrudeSCFParams()
            params.tolerance = 1.0
            params.maxIterations = 50
            force_test.setSCFParameters(params)
            force_test.setAlgorithm(algo)
            
            # Time the calculation
            start = time.time()
            energy = force_test.calculateEnergySCF(state)
            elapsed = time.time() - start
            
            # Calculate error
            error = abs(energy - ref_energy) / abs(ref_energy) * 100.0
            
            results[name]['times'].append(elapsed)
            results[name]['energies'].append(energy)
            results[name]['errors'].append(error)
            
            print(f"    Time: {elapsed:.4f}s, Energy: {energy:.2f}, Error: {error:.2%}")
    
    # Summary table
    print("\n" + "="*80)
    print("SUMMARY: Force Balance Predictor Performance")
    print("="*80)
    
    for n_waters in water_sizes:
        idx = water_sizes.index(n_waters)
        print(f"\n{n_waters} waters:")
        print(f"  {'Algorithm':<15} {'Time (ms)':<12} {'Speedup':<10} {'Error (%)':<10}")
        print(f"  {'-'*15} {'-'*12} {'-'*10} {'-'*10}")
        
        scf_time = results['SCF']['times'][idx]
        for name in ["SCF", "OPT3", "Smart OPT3", "FBP"]:
            time_ms = results[name]['times'][idx] * 1000
            speedup = scf_time / results[name]['times'][idx]
            error = results[name]['errors'][idx]
            print(f"  {name:<15} {time_ms:<12.2f} {speedup:<10.2f} {error:<10.4f}")

if __name__ == "__main__":
    test_algorithm_comparison()
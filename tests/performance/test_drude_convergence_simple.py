#!/usr/bin/env python
"""Simple test to analyze Drude SCF convergence issues"""

import sys
import time
import numpy as np
sys.path.append('/home/zhaomt/gcmc/test107/pygcmc_dev/build/modules/bindings')

import pygcmc

def create_water_box_with_perturbation(n_waters, perturbation=0.0):
    """Create water box with optional random perturbation"""
    n_dim = int(n_waters**(1/3) + 0.5)
    box_size = n_dim * 0.31
    
    state = pygcmc.MCState()
    state.info.box = [box_size, box_size, box_size]
    state.info.cutoff = min(box_size/2 - 0.01, 1.2)
    state.info.setTemperature(300.0)
    
    atoms = []
    residues = []
    
    # PSF parameters
    qO_core = 1.66260
    qD = -1.76260
    qM = -0.957100
    qH = 0.528550
    alpha = 0.0013  # nm^3
    
    mol_id = 0
    for ix in range(n_dim):
        for iy in range(n_dim):
            for iz in range(n_dim):
                if mol_id >= n_waters:
                    break
                    
                x = (ix + 0.5) * 0.31
                y = (iy + 0.5) * 0.31
                z = (iz + 0.5) * 0.31
                
                # Add random perturbation
                if perturbation > 0:
                    x += (np.random.random() - 0.5) * perturbation
                    y += (np.random.random() - 0.5) * perturbation
                    z += (np.random.random() - 0.5) * perturbation
                
                # Oxygen
                o = pygcmc.MCAtom()
                o.x, o.y, o.z = x, y, z
                o.charge = qO_core
                o.type = 0
                atoms.append(o)
                
                # Drude
                d = pygcmc.MCAtom()
                d.x, d.y, d.z = x, y, z
                d.charge = qD
                d.type = 1
                atoms.append(d)
                
                # H1
                h1 = pygcmc.MCAtom()
                h1.x = x + 0.09572
                h1.y = y
                h1.z = z
                h1.charge = qH
                h1.type = 2
                atoms.append(h1)
                
                # H2
                angle = 104.52 * np.pi / 180
                h2 = pygcmc.MCAtom()
                h2.x = x + 0.09572 * np.cos(angle)
                h2.y = y + 0.09572 * np.sin(angle)
                h2.z = z
                h2.charge = qH
                h2.type = 2
                atoms.append(h2)
                
                # M-site
                weights = {'O': 0.786646558, 'H1': 0.106676721, 'H2': 0.106676721}
                m = pygcmc.MCAtom()
                m.x = weights['O'] * o.x + weights['H1'] * h1.x + weights['H2'] * h2.x
                m.y = weights['O'] * o.y + weights['H1'] * h1.y + weights['H2'] * h2.y
                m.z = weights['O'] * o.z + weights['H1'] * h1.z + weights['H2'] * h2.z
                m.charge = qM
                m.type = 3
                atoms.append(m)
                
                # Residue
                res = pygcmc.MCResidue()
                res.atomStart = mol_id * 5
                res.atomCount = 5
                res.active = True
                res.type = 0
                residues.append(res)
                
                mol_id += 1
    
    state.atoms = atoms
    state.activeAtomCount = len(atoms)
    state.residues = residues
    state.activeResidueCount = mol_id
    
    # Force field
    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 4
    ff.numMovementTypes = 4
    
    ljSigma = [0.0] * 16
    ljEps = [0.0] * 16
    ljSigma[0] = 0.318395
    ljEps[0] = 0.88257
    
    ff.ljSigma = ljSigma
    ff.ljEps = ljEps
    state.forcefield = ff
    
    # Drude force
    drude_force = pygcmc.DrudeForce()
    
    for i in range(mol_id):
        drude_idx = i * 5 + 1
        parent_idx = i * 5
        
        drude_force.addParticle(
            drudeIndex=drude_idx,
            parentIndex=parent_idx,
            aniso1Index=-1,
            aniso2Index=-1,
            aniso3Index=-1,
            aniso4Index=-1,
            charge=qD,
            polarizability=alpha,
            aniso12=1.0,
            aniso34=1.0
        )
    
    return state, drude_force

def test_convergence_issues():
    """Test what causes convergence issues"""
    print("=" * 80)
    print("DRUDE SCF CONVERGENCE ISSUE ANALYSIS")
    print("=" * 80)
    
    # Test 1: Effect of tolerance
    print("\n1. EFFECT OF TOLERANCE ON CONVERGENCE")
    print("-" * 50)
    
    state, drude_force = create_water_box_with_perturbation(27)  # 3x3x3
    
    tolerances = [10.0, 1.0, 0.1, 0.01, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8]
    
    print(f"{'Tolerance':>12} | {'Energy':>15} | {'Time (ms)':>10}")
    print("-" * 40)
    
    for tol in tolerances:
        params = pygcmc.DrudeSCFParams()
        params.tolerance = tol
        params.maxIterations = 50
        params.maxDrudeDistance = 0.02
        drude_force.setSCFParameters(params)
        
        # Measure time
        t0 = time.time()
        energy = drude_force.calculateEnergySCF(state)
        t1 = time.time()
        
        print(f"{tol:12.0e} | {energy:15.6f} | {(t1-t0)*1000:10.3f}")
    
    # Test 2: Effect of initial configuration
    print("\n\n2. EFFECT OF INITIAL CONFIGURATION")
    print("-" * 50)
    print("Testing with tolerance = 1e-6")
    
    params = pygcmc.DrudeSCFParams()
    params.tolerance = 1e-6
    params.maxIterations = 50
    params.maxDrudeDistance = 0.02
    
    perturbations = [0.0, 0.001, 0.01, 0.1]
    
    print(f"{'Perturbation':>12} | {'Energy':>15} | {'Std Dev':>10}")
    print("-" * 40)
    
    for pert in perturbations:
        energies = []
        for _ in range(5):
            state, drude_force = create_water_box_with_perturbation(8, pert)
            drude_force.setSCFParameters(params)
            energy = drude_force.calculateEnergySCF(state)
            energies.append(energy)
        
        mean_e = np.mean(energies)
        std_e = np.std(energies)
        print(f"{pert:12.3f} | {mean_e:15.6f} | {std_e:10.6f}")
    
    # Test 3: Effect of max iterations
    print("\n\n3. EFFECT OF MAX ITERATIONS")
    print("-" * 50)
    
    state, drude_force = create_water_box_with_perturbation(27)
    
    max_iters = [5, 10, 20, 50, 100, 200]
    
    print(f"{'Max Iter':>10} | {'Energy':>15} | {'Time (ms)':>10}")
    print("-" * 40)
    
    for max_iter in max_iters:
        params = pygcmc.DrudeSCFParams()
        params.tolerance = 1e-6
        params.maxIterations = max_iter
        params.maxDrudeDistance = 0.02
        drude_force.setSCFParameters(params)
        
        t0 = time.time()
        energy = drude_force.calculateEnergySCF(state)
        t1 = time.time()
        
        print(f"{max_iter:10d} | {energy:15.6f} | {(t1-t0)*1000:10.3f}")
    
    # Test 4: Practical recommendations
    print("\n\n4. RECOMMENDED PARAMETERS vs DEFAULT")
    print("-" * 50)
    
    configs = [
        ("Default (1e-6)", 1e-6, 50),
        ("Loose (0.1)", 0.1, 50),
        ("Moderate (1e-3)", 1e-3, 50),
        ("OpenMM-like (1.0)", 1.0, 100),
        ("CHARMM-like (1e-5)", 1e-5, 100)
    ]
    
    print(f"{'Config':>20} | {'8 waters':>10} | {'27 waters':>10} | {'64 waters':>10}")
    print("-" * 55)
    
    for name, tol, max_iter in configs:
        times = []
        for n_waters in [8, 27, 64]:
            state, drude_force = create_water_box_with_perturbation(n_waters)
            
            params = pygcmc.DrudeSCFParams()
            params.tolerance = tol
            params.maxIterations = max_iter
            params.maxDrudeDistance = 0.02
            drude_force.setSCFParameters(params)
            
            # Time 10 calculations
            t0 = time.time()
            for _ in range(10):
                energy = drude_force.calculateEnergySCF(state)
            t1 = time.time()
            
            avg_time = (t1 - t0) / 10 * 1000  # ms
            times.append(avg_time)
        
        print(f"{name:>20} | {times[0]:10.2f} | {times[1]:10.2f} | {times[2]:10.2f}")

if __name__ == "__main__":
    test_convergence_issues()
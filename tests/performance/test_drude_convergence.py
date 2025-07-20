#!/usr/bin/env python
"""Test Drude SCF convergence with different parameters"""

import sys
import time
import numpy as np
sys.path.append('/home/zhaomt/gcmc/test107/pygcmc_dev/build/modules/bindings')

import pygcmc

def create_water_dimer():
    """Create two SWM4-NDP water molecules"""
    state = pygcmc.MCState()
    state.info.box = [3.0, 3.0, 3.0]
    state.info.cutoff = 1.4
    state.info.setTemperature(300.0)
    
    atoms = []
    residues = []
    
    # PSF parameters
    qO_core = 1.66260
    qD = -1.76260
    qM = -0.957100
    qH = 0.528550
    alpha = 0.0013  # nm^3
    
    # Water 1 at origin
    positions1 = {
        'O': [0.0, 0.0, 0.0],
        'D': [0.0, 0.0, 0.0],
        'H1': [0.09572, 0.0, 0.0],
        'H2': [-0.02399, 0.09277, 0.0],
        'M': [0.0, 0.0127, 0.0]
    }
    
    # Water 2 at 0.3 nm distance
    positions2 = {
        'O': [0.3, 0.0, 0.0],
        'D': [0.3, 0.0, 0.0],
        'H1': [0.39572, 0.0, 0.0],
        'H2': [0.27601, 0.09277, 0.0],
        'M': [0.3, 0.0127, 0.0]
    }
    
    # Add atoms
    for i, (pos1, pos2) in enumerate([(positions1, positions2)]):
        for positions in [pos1, pos2]:
            # O
            o = pygcmc.MCAtom()
            o.x, o.y, o.z = positions['O']
            o.charge = qO_core
            o.type = 0
            atoms.append(o)
            
            # D
            d = pygcmc.MCAtom()
            d.x, d.y, d.z = positions['D']
            d.charge = qD
            d.type = 1
            atoms.append(d)
            
            # H1
            h1 = pygcmc.MCAtom()
            h1.x, h1.y, h1.z = positions['H1']
            h1.charge = qH
            h1.type = 2
            atoms.append(h1)
            
            # H2
            h2 = pygcmc.MCAtom()
            h2.x, h2.y, h2.z = positions['H2']
            h2.charge = qH
            h2.type = 2
            atoms.append(h2)
            
            # M
            m = pygcmc.MCAtom()
            m.x, m.y, m.z = positions['M']
            m.charge = qM
            m.type = 3
            atoms.append(m)
    
    # Residues
    for i in range(2):
        res = pygcmc.MCResidue()
        res.atomStart = i * 5
        res.atomCount = 5
        res.active = True
        res.type = 0
        residues.append(res)
    
    state.atoms = atoms
    state.activeAtomCount = len(atoms)
    state.residues = residues
    state.activeResidueCount = len(residues)
    
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
    drude_force.addParticle(1, 0, -1, -1, -1, -1, qD, alpha, 1.0, 1.0)
    drude_force.addParticle(6, 5, -1, -1, -1, -1, qD, alpha, 1.0, 1.0)
    
    return state, drude_force

def test_convergence_parameters():
    """Test different SCF convergence parameters"""
    print("=" * 80)
    print("DRUDE SCF CONVERGENCE ANALYSIS")
    print("=" * 80)
    
    state, drude_force = create_water_dimer()
    
    # Test different tolerances
    tolerances = [1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9]
    max_iters = [10, 20, 50, 100, 200]
    
    print("\n1. TOLERANCE SWEEP (max_iter=50)")
    print("-" * 50)
    print(f"{'Tolerance':>10} | {'Converged':>10} | {'Iterations':>12} | {'Final RMS':>12} | {'Energy':>12}")
    print("-" * 50)
    
    for tol in tolerances:
        params = pygcmc.DrudeSCFParams()
        params.tolerance = tol
        params.maxIterations = 50
        params.maxDrudeDistance = 0.02
        drude_force.setSCFParameters(params)
        
        # Run SCF
        energy = drude_force.calculateEnergySCF(state)
        stats = drude_force.getSCFStatistics()
        
        print(f"{tol:10.0e} | {stats['converged']:>10} | {stats['iterations']:12d} | {stats['finalRMS']:12.2e} | {energy:12.6f}")
    
    print("\n\n2. MAX ITERATIONS SWEEP (tolerance=1e-6)")
    print("-" * 50)
    print(f"{'Max Iter':>10} | {'Converged':>10} | {'Iterations':>12} | {'Final RMS':>12} | {'Energy':>12}")
    print("-" * 50)
    
    for max_iter in max_iters:
        params = pygcmc.DrudeSCFParams()
        params.tolerance = 1e-6
        params.maxIterations = max_iter
        params.maxDrudeDistance = 0.02
        drude_force.setSCFParameters(params)
        
        energy = drude_force.calculateEnergySCF(state)
        stats = drude_force.getSCFStatistics()
        
        print(f"{max_iter:10d} | {stats['converged']:>10} | {stats['iterations']:12d} | {stats['finalRMS']:12.2e} | {energy:12.6f}")
    
    # Test convergence history
    print("\n\n3. CONVERGENCE HISTORY (tolerance=1e-6, max_iter=100)")
    print("-" * 50)
    
    params = pygcmc.DrudeSCFParams()
    params.tolerance = 1e-6
    params.maxIterations = 100
    params.maxDrudeDistance = 0.02
    drude_force.setSCFParameters(params)
    
    # Move water slightly to create different initial conditions
    for test in range(5):
        # Perturb positions
        dx = 0.001 * test
        state.atoms[5].x = 0.3 + dx
        state.atoms[6].x = 0.3 + dx
        state.atoms[7].x = 0.39572 + dx
        state.atoms[8].x = 0.27601 + dx
        state.atoms[9].x = 0.3 + dx
        
        energy = drude_force.calculateEnergySCF(state)
        stats = drude_force.getSCFStatistics()
        
        print(f"\nTest {test+1}: dx={dx:.3f}")
        print(f"  Converged: {stats['converged']}")
        print(f"  Iterations: {stats['iterations']}")
        print(f"  Final RMS: {stats['finalRMS']:.2e}")
        print(f"  Energy: {energy:.6f}")

def test_system_size_convergence():
    """Test how system size affects convergence"""
    print("\n\n4. SYSTEM SIZE EFFECT ON CONVERGENCE")
    print("-" * 50)
    print(f"{'N waters':>10} | {'Converged':>10} | {'Avg Iter':>10} | {'Max Iter':>10} | {'Avg RMS':>12}")
    print("-" * 50)
    
    for n_dim in [2, 3, 4, 5]:
        n_waters = n_dim ** 3
        
        # Create water box
        state = pygcmc.MCState()
        box_size = n_dim * 0.31
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
        alpha = 0.0013
        
        mol_id = 0
        for ix in range(n_dim):
            for iy in range(n_dim):
                for iz in range(n_dim):
                    x = (ix + 0.5) * 0.31
                    y = (iy + 0.5) * 0.31
                    z = (iz + 0.5) * 0.31
                    
                    # Add 5 atoms per water
                    for i, (charge, dx, dy, dz, atype) in enumerate([
                        (qO_core, 0, 0, 0, 0),
                        (qD, 0, 0, 0, 1),
                        (qH, 0.09572, 0, 0, 2),
                        (qH, -0.02399, 0.09277, 0, 2),
                        (qM, 0, 0.0127, 0, 3)
                    ]):
                        atom = pygcmc.MCAtom()
                        atom.x = x + dx
                        atom.y = y + dy
                        atom.z = z + dz
                        atom.charge = charge
                        atom.type = atype
                        atoms.append(atom)
                    
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
        state.activeResidueCount = len(residues)
        
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
        for i in range(n_waters):
            drude_idx = i * 5 + 1
            parent_idx = i * 5
            drude_force.addParticle(drude_idx, parent_idx, -1, -1, -1, -1, qD, alpha, 1.0, 1.0)
        
        params = pygcmc.DrudeSCFParams()
        params.tolerance = 1e-6
        params.maxIterations = 50
        params.maxDrudeDistance = 0.02
        drude_force.setSCFParameters(params)
        
        # Test 10 different configurations
        iterations = []
        rms_values = []
        converged_count = 0
        
        for test in range(10):
            # Randomly move one water
            mol = np.random.randint(0, n_waters)
            start_atom = mol * 5
            dx = (np.random.random() - 0.5) * 0.01
            dy = (np.random.random() - 0.5) * 0.01
            dz = (np.random.random() - 0.5) * 0.01
            
            for i in [0, 2, 3, 4]:  # Move all except Drude
                state.atoms[start_atom + i].x += dx
                state.atoms[start_atom + i].y += dy
                state.atoms[start_atom + i].z += dz
            
            energy = drude_force.calculateEnergySCF(state)
            stats = drude_force.getSCFStatistics()
            
            iterations.append(stats['iterations'])
            rms_values.append(stats['finalRMS'])
            if stats['converged']:
                converged_count += 1
        
        avg_iter = np.mean(iterations)
        max_iter = np.max(iterations)
        avg_rms = np.mean(rms_values)
        
        print(f"{n_waters:10d} | {converged_count:>4d}/{10:<5d} | {avg_iter:10.1f} | {max_iter:10d} | {avg_rms:12.2e}")

def test_practical_parameters():
    """Test practical parameter recommendations"""
    print("\n\n5. PRACTICAL PARAMETER RECOMMENDATIONS")
    print("-" * 50)
    
    state, drude_force = create_water_dimer()
    
    # Compare with OpenMM-like parameters
    param_sets = [
        ("PyGCMC default", 1e-6, 50),
        ("OpenMM-like", 1.0, 100),
        ("Tight", 1e-8, 200),
        ("Loose", 0.01, 20),
        ("Balanced", 1e-4, 50)
    ]
    
    print(f"{'Config':>15} | {'Tol':>10} | {'MaxIter':>8} | {'Conv':>5} | {'Iter':>5} | {'Time(ms)':>10}")
    print("-" * 70)
    
    for name, tol, max_iter in param_sets:
        params = pygcmc.DrudeSCFParams()
        params.tolerance = tol
        params.maxIterations = max_iter
        params.maxDrudeDistance = 0.02
        drude_force.setSCFParameters(params)
        
        # Time the calculation
        times = []
        converged_all = True
        total_iters = 0
        
        for _ in range(100):
            t0 = time.time()
            energy = drude_force.calculateEnergySCF(state)
            t1 = time.time()
            times.append(t1 - t0)
            
            stats = drude_force.getSCFStatistics()
            if not stats['converged']:
                converged_all = False
            total_iters += stats['iterations']
        
        avg_time = np.mean(times) * 1000
        avg_iters = total_iters / 100
        
        print(f"{name:>15} | {tol:10.0e} | {max_iter:8d} | {'Yes' if converged_all else 'No':>5} | {avg_iters:5.1f} | {avg_time:10.3f}")

if __name__ == "__main__":
    test_convergence_parameters()
    test_system_size_convergence()
    test_practical_parameters()
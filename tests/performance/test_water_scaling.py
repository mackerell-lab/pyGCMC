#!/usr/bin/env python
"""Test actual performance scaling with different water box sizes"""

import sys
import time
import numpy as np
sys.path.append('/home/zhaomt/gcmc/test107/pygcmc_dev/build/modules/bindings')

import pygcmc

def create_water_box(n_per_dim):
    """Create a cubic water box with n_per_dim waters per dimension"""
    n_waters = n_per_dim ** 3
    box_size = n_per_dim * 0.31  # 0.31 nm spacing
    
    state = pygcmc.MCState()
    state.info.box = [box_size, box_size, box_size]
    state.info.cutoff = min(box_size/2 - 0.01, 1.2)
    state.info.setTemperature(300.0)
    
    atoms = []
    residues = []
    
    # Water geometry
    angle = 104.52 * np.pi / 180
    
    mol_id = 0
    for ix in range(n_per_dim):
        for iy in range(n_per_dim):
            for iz in range(n_per_dim):
                x = (ix + 0.5) * 0.31
                y = (iy + 0.5) * 0.31
                z = (iz + 0.5) * 0.31
                
                # Oxygen
                o = pygcmc.MCAtom()
                o.x, o.y, o.z = x, y, z
                o.charge = -0.834
                o.type = 0
                atoms.append(o)
                
                # H1
                h1 = pygcmc.MCAtom()
                h1.x = x + 0.09572
                h1.y = y
                h1.z = z
                h1.charge = 0.417
                h1.type = 1
                atoms.append(h1)
                
                # H2
                h2 = pygcmc.MCAtom()
                h2.x = x + 0.09572 * np.cos(angle)
                h2.y = y + 0.09572 * np.sin(angle)
                h2.z = z
                h2.charge = 0.417
                h2.type = 1
                atoms.append(h2)
                
                # Residue
                res = pygcmc.MCResidue()
                res.atomStart = mol_id * 3
                res.atomCount = 3
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
    ff.numTotalTypes = 2
    ff.numMovementTypes = 2
    ff.ljSigma = [0.315, 0.0, 0.0, 0.0]
    ff.ljEps = [0.636, 0.0, 0.0, 0.0]
    state.forcefield = ff
    
    return state, n_waters, box_size

def test_performance(state, n_waters, n_steps=1000):
    """Test performance with movement and energy calculations"""
    target_mol = n_waters // 2
    start_atom = target_mol * 3
    
    # Warmup
    for _ in range(10):
        pygcmc.computeSystemEnergyCutoff(state)
    
    # Test
    times = []
    for step in range(n_steps):
        # Move water slightly
        dx = (np.random.random() - 0.5) * 0.001
        dy = (np.random.random() - 0.5) * 0.001
        dz = (np.random.random() - 0.5) * 0.001
        
        for i in range(3):
            state.atoms[start_atom + i].x += dx
            state.atoms[start_atom + i].y += dy
            state.atoms[start_atom + i].z += dz
        
        # Time energy calculation
        t0 = time.time()
        pygcmc.computeSystemEnergyCutoff(state)
        t1 = time.time()
        times.append(t1 - t0)
    
    return np.array(times)

def main():
    print("=" * 70)
    print("Water Box Scaling Test - Actual Performance Measurements")
    print("=" * 70)
    
    # Test configurations: (n_per_dim, expected_n_waters)
    configs = [
        (2, 8),      # 2x2x2 = 8
        (3, 27),     # 3x3x3 = 27
        (4, 64),     # 4x4x4 = 64
        (5, 125),    # 5x5x5 = 125
        (6, 216),    # 6x6x6 = 216
        (7, 343),    # 7x7x7 = 343
        (8, 512),    # 8x8x8 = 512
        (10, 1000),  # 10x10x10 = 1000
    ]
    
    results = []
    
    print("\nCreating and testing water boxes...\n")
    print(f"{'N waters':>8} | {'Box (nm)':>8} | {'Cutoff':>7} | {'ms/step':>8} | {'steps/s':>8} | {'vs 8-water':>10}")
    print("-" * 70)
    
    base_time = None
    
    for n_per_dim, expected_n in configs:
        # Skip very large systems if they take too long
        if n_per_dim > 6 and base_time is not None and base_time * (expected_n/8)**2 > 0.1:
            print(f"{expected_n:8d} | {'skipped - would take too long':>45}")
            continue
            
        # Create system
        state, n_waters, box_size = create_water_box(n_per_dim)
        assert n_waters == expected_n
        
        # Test performance
        n_test_steps = min(1000, max(100, int(1000 * 8 / n_waters)))  # Fewer steps for larger systems
        times = test_performance(state, n_waters, n_test_steps)
        
        avg_time = np.mean(times) * 1000  # Convert to ms
        std_time = np.std(times) * 1000
        steps_per_sec = 1000 / avg_time
        
        if base_time is None:
            base_time = avg_time
            scaling = 1.0
        else:
            scaling = avg_time / base_time
        
        results.append({
            'n_waters': n_waters,
            'box_size': box_size,
            'cutoff': state.info.cutoff,
            'avg_time_ms': avg_time,
            'std_time_ms': std_time,
            'steps_per_sec': steps_per_sec,
            'scaling': scaling
        })
        
        print(f"{n_waters:8d} | {box_size:8.2f} | {state.info.cutoff:7.2f} | "
              f"{avg_time:8.3f} | {steps_per_sec:8.1f} | {scaling:10.1f}x")
    
    # Analysis
    print("\n" + "=" * 70)
    print("SCALING ANALYSIS")
    print("=" * 70)
    
    if len(results) > 1:
        # Check scaling law
        n_vals = np.array([r['n_waters'] for r in results])
        t_vals = np.array([r['avg_time_ms'] for r in results])
        
        # Fit power law: t = a * n^b
        log_n = np.log(n_vals)
        log_t = np.log(t_vals)
        
        # Linear regression in log space
        A = np.vstack([log_n, np.ones(len(log_n))]).T
        b, log_a = np.linalg.lstsq(A, log_t, rcond=None)[0]
        a = np.exp(log_a)
        
        print(f"\nFitted scaling law: time = {a:.3e} * N^{b:.2f}")
        print(f"Expected O(N²) scaling would give exponent = 2.00")
        print(f"Actual scaling exponent: {b:.2f}")
        
        # Predictions vs actual
        print("\nPredicted vs Actual times:")
        print(f"{'N waters':>8} | {'Actual (ms)':>12} | {'Predicted':>12} | {'Error %':>8}")
        print("-" * 50)
        
        for r in results:
            n = r['n_waters']
            actual = r['avg_time_ms']
            predicted = a * (n ** b)
            error = 100 * (predicted - actual) / actual
            print(f"{n:8d} | {actual:12.3f} | {predicted:12.3f} | {error:8.1f}%")
    
    # Time estimates for 20,000 steps
    print("\n" + "=" * 70)
    print("TIME ESTIMATES FOR 20,000 STEPS")
    print("=" * 70)
    
    print(f"{'N waters':>8} | {'Time (seconds)':>15} | {'Time (minutes)':>15}")
    print("-" * 45)
    
    for r in results:
        time_20k = r['avg_time_ms'] * 20000 / 1000  # seconds
        print(f"{r['n_waters']:8d} | {time_20k:15.1f} | {time_20k/60:15.1f}")

if __name__ == "__main__":
    main()
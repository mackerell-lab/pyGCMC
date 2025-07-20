#!/usr/bin/env python
"""Actually test ALL water box sizes - no predictions, only measurements"""

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

def test_performance_fixed_time(state, n_waters, max_time=5.0):
    """Test performance for a fixed amount of time"""
    target_mol = n_waters // 2
    start_atom = target_mol * 3
    
    # Warmup
    for _ in range(10):
        pygcmc.computeSystemEnergyCutoff(state)
    
    # Test for fixed time
    times = []
    start_time = time.time()
    steps = 0
    
    while True:
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
        steps += 1
        
        # Check if we've run long enough
        if t1 - start_time > max_time:
            break
        
        # For very fast systems, make sure we get enough samples
        if steps >= 10000:
            break
    
    return np.array(times), steps

def main():
    print("=" * 80)
    print("COMPLETE Water Box Performance Test - ALL SIZES MEASURED")
    print("=" * 80)
    
    # All configurations to test
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
    
    print("\nMeasuring actual performance for each system size...")
    print("(Each test runs for up to 5 seconds)\n")
    print(f"{'N waters':>8} | {'Box (nm)':>8} | {'Steps':>8} | {'ms/step':>10} | {'steps/s':>10} | {'20k steps':>12}")
    print("-" * 80)
    
    all_results = []
    
    for n_per_dim, expected_n in configs:
        # Create system
        print(f"\rCreating {expected_n:4d} water system...", end='', flush=True)
        state, n_waters, box_size = create_water_box(n_per_dim)
        assert n_waters == expected_n
        
        # Test performance
        print(f"\rTesting  {expected_n:4d} water system...", end='', flush=True)
        times, steps = test_performance_fixed_time(state, n_waters, max_time=5.0)
        
        avg_time = np.mean(times) * 1000  # Convert to ms
        std_time = np.std(times) * 1000
        steps_per_sec = 1000 / avg_time
        time_20k = avg_time * 20000 / 1000  # seconds
        
        # Format time nicely
        if time_20k < 60:
            time_str = f"{time_20k:.1f} s"
        elif time_20k < 3600:
            time_str = f"{time_20k/60:.1f} min"
        else:
            time_str = f"{time_20k/3600:.1f} hr"
        
        print(f"\r{n_waters:8d} | {box_size:8.2f} | {steps:8d} | {avg_time:10.3f} | {steps_per_sec:10.1f} | {time_str:>12}")
        
        all_results.append({
            'n_waters': n_waters,
            'box_size': box_size,
            'steps': steps,
            'avg_time_ms': avg_time,
            'std_time_ms': std_time,
            'steps_per_sec': steps_per_sec,
            'time_20k_sec': time_20k
        })
    
    # Summary
    print("\n" + "=" * 80)
    print("SUMMARY - All Actual Measurements (No Predictions)")
    print("=" * 80)
    
    print("\nActual scaling observed:")
    print(f"{'N waters':>8} | {'Actual ms/step':>15} | {'Scaling vs 8-water':>20}")
    print("-" * 50)
    
    base_time = all_results[0]['avg_time_ms']
    for r in all_results:
        scaling = r['avg_time_ms'] / base_time
        print(f"{r['n_waters']:8d} | {r['avg_time_ms']:15.3f} | {scaling:20.1f}x")
    
    # Check if it follows N^2
    print("\nScaling analysis:")
    n_vals = np.array([r['n_waters'] for r in all_results])
    t_vals = np.array([r['avg_time_ms'] for r in all_results])
    
    # Expected N^2 scaling
    expected_scaling = (n_vals / n_vals[0]) ** 2
    actual_scaling = t_vals / t_vals[0]
    
    print(f"{'N waters':>8} | {'Expected (N²)':>15} | {'Actual scaling':>15} | {'Ratio':>10}")
    print("-" * 60)
    for i, n in enumerate(n_vals):
        ratio = actual_scaling[i] / expected_scaling[i]
        print(f"{n:8d} | {expected_scaling[i]:15.1f}x | {actual_scaling[i]:15.1f}x | {ratio:10.2f}")
    
    # Save detailed results
    with open('water_performance_all_measured.txt', 'w') as f:
        f.write("Complete Water Box Performance - All Measured\n")
        f.write("=" * 60 + "\n\n")
        for r in all_results:
            f.write(f"{r['n_waters']} waters:\n")
            f.write(f"  Box size: {r['box_size']:.2f} nm\n")
            f.write(f"  Steps tested: {r['steps']}\n")
            f.write(f"  Time per step: {r['avg_time_ms']:.3f} ± {r['std_time_ms']:.3f} ms\n")
            f.write(f"  Steps per second: {r['steps_per_sec']:.1f}\n")
            f.write(f"  Time for 20k steps: {r['time_20k_sec']:.1f} seconds\n")
            f.write("\n")
    
    print(f"\nDetailed results saved to water_performance_all_measured.txt")

if __name__ == "__main__":
    main()
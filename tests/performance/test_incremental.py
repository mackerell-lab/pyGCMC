#!/usr/bin/env python
"""Incremental performance test that can be run in parts to avoid timeout"""

import sys
import time
import json
import os
import numpy as np
sys.path.append('/home/zhaomt/gcmc/test107/pygcmc_dev/build/modules/bindings')

import pygcmc

# State file to save progress
STATE_FILE = 'performance_test_state.json'
RESULTS_FILE = 'performance_results.json'

def load_state():
    """Load test state from file"""
    if os.path.exists(STATE_FILE):
        with open(STATE_FILE, 'r') as f:
            return json.load(f)
    else:
        # Initial state
        return {
            'completed': [],
            'current': None,
            'all_systems': [
                {'n_dim': 2, 'n_waters': 8},
                {'n_dim': 3, 'n_waters': 27},
                {'n_dim': 4, 'n_waters': 64},
                {'n_dim': 5, 'n_waters': 125},
                {'n_dim': 6, 'n_waters': 216},
                {'n_dim': 7, 'n_waters': 343},
                {'n_dim': 8, 'n_waters': 512},
                {'n_dim': 10, 'n_waters': 1000},
            ]
        }

def save_state(state):
    """Save test state to file"""
    with open(STATE_FILE, 'w') as f:
        json.dump(state, f, indent=2)

def load_results():
    """Load existing results"""
    if os.path.exists(RESULTS_FILE):
        with open(RESULTS_FILE, 'r') as f:
            return json.load(f)
    else:
        return []

def save_results(results):
    """Save results to file"""
    with open(RESULTS_FILE, 'w') as f:
        json.dump(results, f, indent=2)

def create_water_box(n_per_dim):
    """Create water box (same as before)"""
    n_waters = n_per_dim ** 3
    box_size = n_per_dim * 0.31
    
    state = pygcmc.MCState()
    state.info.box = [box_size, box_size, box_size]
    state.info.cutoff = min(box_size/2 - 0.01, 1.2)
    state.info.setTemperature(300.0)
    
    atoms = []
    residues = []
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

def test_system(n_per_dim, n_waters, max_time=90):
    """Test a single system with time limit"""
    print(f"\nTesting {n_waters} water system ({n_per_dim}x{n_per_dim}x{n_per_dim})...")
    
    # Create system
    state, actual_n, box_size = create_water_box(n_per_dim)
    
    # Quick test with fewer steps for large systems
    if n_waters <= 125:
        n_steps = 1000
    elif n_waters <= 343:
        n_steps = 500
    else:
        n_steps = 200
    
    print(f"  Running {n_steps} steps...")
    
    target_mol = n_waters // 2
    start_atom = target_mol * 3
    
    # Warmup
    for _ in range(10):
        pygcmc.computeSystemEnergyCutoff(state)
    
    # Test
    times = []
    start_time = time.time()
    
    for step in range(n_steps):
        # Move water
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
        
        # Check time limit
        if t1 - start_time > max_time:
            print(f"  Time limit reached after {step+1} steps")
            break
    
    avg_time_ms = np.mean(times) * 1000
    std_time_ms = np.std(times) * 1000
    
    result = {
        'n_waters': n_waters,
        'n_dim': n_per_dim,
        'box_size': box_size,
        'steps_tested': len(times),
        'avg_time_ms': avg_time_ms,
        'std_time_ms': std_time_ms,
        'steps_per_sec': 1000 / avg_time_ms,
        'time_20k_sec': avg_time_ms * 20000 / 1000
    }
    
    print(f"  Average: {avg_time_ms:.3f} ms/step")
    print(f"  Steps/sec: {result['steps_per_sec']:.1f}")
    
    return result

def main():
    print("=" * 70)
    print("INCREMENTAL PERFORMANCE TEST")
    print("=" * 70)
    
    # Load state
    state = load_state()
    results = load_results()
    
    # Show progress
    print(f"\nCompleted systems: {[s['n_waters'] for s in state['completed']]}")
    remaining = [s for s in state['all_systems'] if s not in state['completed']]
    print(f"Remaining systems: {[s['n_waters'] for s in remaining]}")
    
    if not remaining:
        print("\nAll tests completed! Showing final results:")
        print_results(results)
        return
    
    # Run next system
    start_time = time.time()
    max_runtime = 100  # seconds, leaving buffer
    
    for system in remaining:
        if time.time() - start_time > max_runtime:
            print("\nTime limit approaching, saving progress...")
            break
            
        result = test_system(system['n_dim'], system['n_waters'])
        results.append(result)
        state['completed'].append(system)
        
        # Save progress after each system
        save_state(state)
        save_results(results)
        
        print(f"  Saved. Progress: {len(state['completed'])}/{len(state['all_systems'])}")
    
    # Show current results
    print("\n" + "=" * 70)
    print("RESULTS SO FAR:")
    print_results(results)
    
    if len(state['completed']) < len(state['all_systems']):
        print(f"\nRun the script again to continue testing remaining systems.")
    else:
        print(f"\nAll systems tested! Total time across runs: check timestamps")
        # Clean up state file
        if os.path.exists(STATE_FILE):
            os.remove(STATE_FILE)

def print_results(results):
    """Print formatted results"""
    print(f"\n{'N waters':>8} | {'ms/step':>10} | {'steps/sec':>10} | {'20k steps':>12}")
    print("-" * 50)
    
    for r in sorted(results, key=lambda x: x['n_waters']):
        time_20k = r['time_20k_sec']
        if time_20k < 60:
            time_str = f"{time_20k:.1f} s"
        elif time_20k < 3600:
            time_str = f"{time_20k/60:.1f} min"
        else:
            time_str = f"{time_20k/3600:.1f} hr"
        
        print(f"{r['n_waters']:8d} | {r['avg_time_ms']:10.3f} | {r['steps_per_sec']:10.1f} | {time_str:>12}")

if __name__ == "__main__":
    main()
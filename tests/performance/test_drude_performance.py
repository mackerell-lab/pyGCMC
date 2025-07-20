#!/usr/bin/env python
"""Test Drude water model performance"""

import sys
import time
import numpy as np
sys.path.append('/home/zhaomt/gcmc/test107/pygcmc_dev/build/modules/bindings')

import pygcmc

def create_drude_water_box(n_per_dim):
    """Create a box of SWM4-NDP water molecules with Drude particles"""
    n_waters = n_per_dim ** 3
    box_size = n_per_dim * 0.31
    
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
    alpha = 0.0013  # nm^3 (1.3 Å^3)
    
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
                o.charge = qO_core
                o.type = 0
                atoms.append(o)
                
                # Drude on oxygen
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
                
                # M-site (virtual)
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
    state.activeResidueCount = len(residues)
    
    # Force field
    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 4  # O, D, H, M
    ff.numMovementTypes = 4
    
    # LJ only on oxygen
    ljSigma = [0.0] * 16
    ljEps = [0.0] * 16
    ljSigma[0] = 0.318395  # O-O
    ljEps[0] = 0.88257     # O-O
    
    ff.ljSigma = ljSigma
    ff.ljEps = ljEps
    state.forcefield = ff
    
    # Create Drude force
    drude_force = pygcmc.DrudeForce()
    
    # Add Drude particles
    for i in range(n_waters):
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
    
    # Use default SCF parameters (tolerance = 1.0)
    # No need to set parameters explicitly - defaults are good
    
    return state, drude_force, n_waters, box_size

def test_drude_performance():
    """Test Drude model performance for different system sizes"""
    print("=" * 70)
    print("DRUDE WATER MODEL PERFORMANCE TEST")
    print("=" * 70)
    
    # Test smaller systems due to expected slower performance
    test_sizes = [
        (2, 8),    # Very small
        (3, 27),   # Small
        (4, 64),   # Medium
        (5, 125),  # Practical limit?
    ]
    
    results = []
    
    print(f"\n{'N waters':>8} | {'ms/step':>10} | {'steps/sec':>10} | {'vs non-Drude':>12} | {'20k steps':>12}")
    print("-" * 70)
    
    # Non-Drude reference times (from previous measurements)
    ref_times = {8: 0.003, 27: 0.031, 64: 0.159, 125: 0.643}
    
    for n_dim, n_waters in test_sizes:
        print(f"\rTesting {n_waters:4d} waters...", end='', flush=True)
        
        # Create system
        state, drude_force, actual_n, box_size = create_drude_water_box(n_dim)
        
        # Test with fewer steps
        n_steps = min(100, 800 // n_waters)  # Fewer steps for larger systems
        target_mol = 0  # First water
        start_atom = 0
        
        # Warmup
        for _ in range(5):
            drude_force.calculateEnergySCF(state)
        
        # Test
        times = []
        for step in range(n_steps):
            # Move water
            dx = (np.random.random() - 0.5) * 0.001
            dy = (np.random.random() - 0.5) * 0.001
            dz = (np.random.random() - 0.5) * 0.001
            
            # Move all atoms of first water (except Drude)
            for i in [0, 2, 3, 4]:  # O, H1, H2, M
                state.atoms[start_atom + i].x += dx
                state.atoms[start_atom + i].y += dy
                state.atoms[start_atom + i].z += dz
            
            # Time energy + SCF
            t0 = time.time()
            energy = drude_force.calculateEnergySCF(state)
            t1 = time.time()
            times.append(t1 - t0)
        
        avg_time_ms = np.mean(times) * 1000
        steps_per_sec = 1000 / avg_time_ms
        slowdown = avg_time_ms / ref_times[n_waters]
        time_20k = avg_time_ms * 20000 / 1000
        
        if time_20k < 60:
            time_str = f"{time_20k:.1f} s"
        elif time_20k < 3600:
            time_str = f"{time_20k/60:.1f} min"
        else:
            time_str = f"{time_20k/3600:.1f} hr"
        
        print(f"\r{n_waters:8d} | {avg_time_ms:10.3f} | {steps_per_sec:10.1f} | {slowdown:12.1f}x | {time_str:>12}")
        
        results.append({
            'n_waters': n_waters,
            'avg_time_ms': avg_time_ms,
            'slowdown': slowdown,
            'time_20k': time_str
        })
        
        # Stop if getting too slow
        if avg_time_ms > 50:  # 50ms per step is quite slow
            print("\n\nStopping test - larger systems would take too long")
            break
    
    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY:")
    print("-" * 70)
    print(f"Drude is approximately {np.mean([r['slowdown'] for r in results]):.0f}x slower than non-Drude")
    print("\nPractical system sizes for Drude:")
    for r in results:
        if r['n_waters'] <= 125:
            print(f"  {r['n_waters']} waters: {r['time_20k']} for 20k steps")

if __name__ == "__main__":
    test_drude_performance()
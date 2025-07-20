#!/usr/bin/env python
"""Debug why large systems are slower and still have warnings"""

import sys
import time
import numpy as np
sys.path.append('/home/zhaomt/gcmc/test107/pygcmc_dev/build/modules/bindings')

import pygcmc

def create_drude_water_box(n_per_dim, offset=0.001):
    """Create a box of SWM4-NDP water molecules"""
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
    alpha = 0.0013
    
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
                
                # Drude - with specified offset
                d = pygcmc.MCAtom()
                d.x = x + offset * (np.random.random() - 0.5)
                d.y = y + offset * (np.random.random() - 0.5)
                d.z = z + offset * (np.random.random() - 0.5)
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
    
    return state, drude_force, n_waters

print("=" * 80)
print("DEBUGGING LARGE SYSTEM CONVERGENCE ISSUES")
print("=" * 80)

# Test 1: Compare different system sizes with same tolerance
print("\n1. SYSTEM SIZE EFFECT (tolerance = 1.0)")
print("-" * 60)

sizes = [(2, 8), (3, 27), (4, 64), (5, 125)]

for n_dim, n_waters in sizes:
    state, drude_force, _ = create_drude_water_box(n_dim)
    
    # First calculation - should take longer
    t0 = time.time()
    energy1 = drude_force.calculateEnergySCF(state)
    t1 = time.time()
    
    # Second calculation - should be faster (already converged)
    t2 = time.time()
    energy2 = drude_force.calculateEnergySCF(state)
    t3 = time.time()
    
    print(f"\n{n_waters} waters:")
    print(f"  First calc:  {(t1-t0)*1000:.1f} ms, E = {energy1:.2f} kJ/mol")
    print(f"  Second calc: {(t3-t2)*1000:.1f} ms, E = {energy2:.2f} kJ/mol")
    print(f"  Speedup: {(t1-t0)/(t3-t2):.1f}x")

# Test 2: Effect of initial offset
print("\n\n2. INITIAL OFFSET EFFECT (125 waters)")
print("-" * 60)

offsets = [0.0001, 0.001, 0.01, 0.1]

for offset in offsets:
    state, drude_force, _ = create_drude_water_box(5, offset)
    
    times = []
    for _ in range(5):
        t0 = time.time()
        energy = drude_force.calculateEnergySCF(state)
        t1 = time.time()
        times.append((t1-t0)*1000)
    
    avg_time = np.mean(times)
    print(f"\nOffset = {offset} nm:")
    print(f"  Avg time: {avg_time:.1f} ms")
    print(f"  Energy: {energy:.2f} kJ/mol")

# Test 3: Different tolerances
print("\n\n3. TOLERANCE EFFECT (125 waters)")
print("-" * 60)

state, drude_force, _ = create_drude_water_box(5)

tolerances = [100.0, 10.0, 1.0, 0.1]

for tol in tolerances:
    params = pygcmc.DrudeSCFParams()
    params.tolerance = tol
    params.maxIterations = 50
    params.maxDrudeDistance = 0.02
    drude_force.setSCFParameters(params)
    
    times = []
    for _ in range(5):
        # Perturb atoms slightly
        state.atoms[0].x += 0.0001
        
        t0 = time.time()
        energy = drude_force.calculateEnergySCF(state)
        t1 = time.time()
        times.append((t1-t0)*1000)
    
    avg_time = np.mean(times)
    print(f"\nTolerance = {tol}:")
    print(f"  Avg time: {avg_time:.1f} ms")
    print(f"  Energy: {energy:.2f} kJ/mol")

# Test 4: Track specific warning case
print("\n\n4. TRACKING SPECIFIC WARNING CASE")
print("-" * 60)

# Create fresh system
state, drude_force, _ = create_drude_water_box(5)

# Move one water significantly
print("\nMoving water 0 by 0.1 nm...")
for i in [0, 2, 3, 4]:  # Move O, H1, H2, M (not Drude)
    state.atoms[i].x += 0.1

t0 = time.time()
energy = drude_force.calculateEnergySCF(state)
t1 = time.time()

print(f"Time: {(t1-t0)*1000:.1f} ms")
print(f"Energy: {energy:.2f} kJ/mol")
print("\nIf warning appeared above, it's due to large perturbation requiring many iterations")
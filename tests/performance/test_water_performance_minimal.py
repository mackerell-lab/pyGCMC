#!/usr/bin/env python
"""Minimal water box performance test"""

import sys
import time
import numpy as np
sys.path.append('/home/zhaomt/gcmc/test107/pygcmc_dev/build/modules/bindings')

import pygcmc

# Create small water box - 2x2x2 = 8 waters
state = pygcmc.MCState()
state.info.box = [1.0, 1.0, 1.0]  # 1 nm box
state.info.cutoff = 0.4
state.info.setTemperature(300.0)

# Create 8 water molecules
atoms = []
residues = []
water_pos = [
    (0.25, 0.25, 0.25), (0.75, 0.25, 0.25),
    (0.25, 0.75, 0.25), (0.75, 0.75, 0.25),
    (0.25, 0.25, 0.75), (0.75, 0.25, 0.75),
    (0.25, 0.75, 0.75), (0.75, 0.75, 0.75)
]

for i, (x, y, z) in enumerate(water_pos):
    # Oxygen
    o = pygcmc.MCAtom()
    o.x, o.y, o.z = x, y, z
    o.charge = -0.834
    o.type = 0
    atoms.append(o)
    
    # H1
    h1 = pygcmc.MCAtom()
    h1.x, h1.y, h1.z = x + 0.09572, y, z
    h1.charge = 0.417
    h1.type = 1
    atoms.append(h1)
    
    # H2
    h2 = pygcmc.MCAtom()
    h2.x, h2.y, h2.z = x + 0.078, y + 0.055, z
    h2.charge = 0.417
    h2.type = 1
    atoms.append(h2)
    
    # Residue
    res = pygcmc.MCResidue()
    res.atomStart = i * 3
    res.atomCount = 3
    res.active = True
    res.type = 0
    residues.append(res)

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

print("Small water box test:")
print(f"  Waters: 8")
print(f"  Atoms: 24")
print(f"  Box: 1x1x1 nm")
print(f"  Cutoff: 0.4 nm")

# Initial energy
pygcmc.computeSystemEnergyCutoff(state)
energy = sum(r.energy_elec + r.energy_vdw for r in state.residues)
print(f"\nInitial energy: {energy:.3f} kJ/mol")

# Time 10000 calculations
print("\nTiming 10000 energy calculations...")
start = time.time()

for _ in range(10000):
    pygcmc.computeSystemEnergyCutoff(state)

end = time.time()
total_time = end - start

print(f"Time for 10000 calculations: {total_time:.3f} s")
print(f"Time per calculation: {total_time/10000*1000:.3f} ms")
print(f"Calculations per second: {10000/total_time:.1f}")

# Now do a simple movement test
print("\n\nMovement test (5000 steps)...")
target_mol = 0  # First water
start_atom = 0

times = []
for step in range(5000):
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

avg_time = np.mean(times) * 1000  # ms
print(f"Average time per move+energy: {avg_time:.3f} ms")
print(f"Steps per second: {1000/avg_time:.1f}")

# Estimate for larger systems
print("\n\nEstimated performance for larger systems:")
print("(Assuming O(N²) scaling for cutoff calculations)")

for n_waters in [27, 64, 125, 216, 343, 512, 1000]:
    scale_factor = (n_waters / 8) ** 2
    est_time = avg_time * scale_factor
    print(f"  {n_waters:4d} waters: ~{est_time:8.1f} ms/step, {1000/est_time:6.1f} steps/s")
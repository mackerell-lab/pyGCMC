#!/usr/bin/env python
"""Test the effect of water spacing on energy"""

import sys
import numpy as np
sys.path.append('/home/zhaomt/gcmc/test107/pygcmc_dev/build/modules/bindings')

import pygcmc

def create_water_box_with_spacing(n_waters, spacing):
    """Create water box with specified spacing"""
    n_dim = int(n_waters**(1/3) + 0.5)
    box_size = n_dim * spacing
    
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
    for i in range(n_waters):
        x = (i % n_dim + 0.5) * spacing
        y = ((i // n_dim) % n_dim + 0.5) * spacing
        z = (i // (n_dim * n_dim) + 0.5) * spacing
        
        # Add small random perturbation to break symmetry
        x += (np.random.random() - 0.5) * 0.02
        y += (np.random.random() - 0.5) * 0.02
        z += (np.random.random() - 0.5) * 0.02
        
        # Oxygen
        o = pygcmc.MCAtom()
        o.x, o.y, o.z = x, y, z
        o.charge = qO_core
        o.type = 0
        atoms.append(o)
        
        # Drude - with proper offset
        d = pygcmc.MCAtom()
        d.x = x + 0.001 * (np.random.random() - 0.5)
        d.y = y + 0.001 * (np.random.random() - 0.5)
        d.z = z + 0.001 * (np.random.random() - 0.5)
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

print("=" * 70)
print("WATER SPACING EFFECT ON ENERGY")
print("=" * 70)

# Test different spacings
spacings = [0.31, 0.35, 0.40, 0.45, 0.50]  # nm
n_waters = 27  # 3x3x3

print(f"\nTesting {n_waters} water molecules with different spacings:")
print("-" * 60)
print(f"{'Spacing (nm)':>12} | {'Energy (kJ/mol)':>15} | {'Per molecule':>12} | {'Density (g/cm³)':>15}")
print("-" * 60)

for spacing in spacings:
    state, drude_force = create_water_box_with_spacing(n_waters, spacing)
    
    # Calculate energy
    energy = drude_force.calculateEnergySCF(state)
    
    # Calculate density
    box_volume = (spacing * 3)**3  # nm³
    mass = n_waters * 18.015  # g/mol
    density = mass / (box_volume * 0.6022)  # g/cm³
    
    print(f"{spacing:12.2f} | {energy:15.2f} | {energy/n_waters:12.2f} | {density:15.3f}")

# Also test the minimum distance between atoms
print("\n\nMinimum distances in 0.31 nm spacing:")
print("-" * 40)

state, drude_force = create_water_box_with_spacing(27, 0.31)

# Check O-O distances
min_oo_dist = float('inf')
for i in range(27):
    o1_idx = i * 5
    for j in range(i+1, 27):
        o2_idx = j * 5
        dx = state.atoms[o1_idx].x - state.atoms[o2_idx].x
        dy = state.atoms[o1_idx].y - state.atoms[o2_idx].y
        dz = state.atoms[o1_idx].z - state.atoms[o2_idx].z
        dist = np.sqrt(dx*dx + dy*dy + dz*dz)
        min_oo_dist = min(min_oo_dist, dist)

print(f"Minimum O-O distance: {min_oo_dist:.3f} nm")
print(f"O-O LJ sigma: {0.318395:.3f} nm")
print(f"Ratio: {min_oo_dist/0.318395:.2f}")

if min_oo_dist < 0.318395:
    print("WARNING: O-O distance is less than LJ sigma - strong repulsion!")

# Check for any very close atoms
min_dist = float('inf')
for i in range(len(state.atoms)):
    for j in range(i+1, len(state.atoms)):
        # Skip intramolecular
        if i // 5 == j // 5:
            continue
        dx = state.atoms[i].x - state.atoms[j].x
        dy = state.atoms[i].y - state.atoms[j].y
        dz = state.atoms[i].z - state.atoms[j].z
        dist = np.sqrt(dx*dx + dy*dy + dz*dz)
        min_dist = min(min_dist, dist)

print(f"\nMinimum intermolecular distance: {min_dist:.3f} nm")
if min_dist < 0.15:
    print("WARNING: Very close atoms detected - high energy expected!")

# Test with looser tolerance for large system
print("\n\nTesting 125 waters with different tolerances:")
print("-" * 60)

state, drude_force = create_water_box_with_spacing(125, 0.40)  # Use larger spacing

tolerances = [100.0, 10.0, 1.0]

for tol in tolerances:
    params = pygcmc.DrudeSCFParams()
    params.tolerance = tol
    params.maxIterations = 100  # More iterations
    params.maxDrudeDistance = 0.02
    drude_force.setSCFParameters(params)
    
    import time
    t0 = time.time()
    energy = drude_force.calculateEnergySCF(state)
    t1 = time.time()
    
    print(f"Tolerance = {tol:5.1f}: Energy = {energy:10.2f} kJ/mol, Time = {(t1-t0)*1000:6.1f} ms")
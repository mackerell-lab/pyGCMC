#!/usr/bin/env python
"""Test with pre-optimized water configuration"""

import sys
import time
import numpy as np
sys.path.append('/home/zhaomt/gcmc/test107/pygcmc_dev/build/modules/bindings')

import pygcmc

def create_optimized_water_box(n_waters):
    """Create water box with reasonable spacing and random orientations"""
    # Calculate box size for ~1 g/cm³ density
    # 1 g/cm³ = 1000 kg/m³
    # For n water molecules: V = n * M / (ρ * NA)
    # M = 18.015 g/mol, NA = 6.022e23, ρ = 1 g/cm³
    volume_per_water = 18.015 / (1.0 * 602.2)  # nm³
    total_volume = n_waters * volume_per_water
    box_size = total_volume ** (1/3)
    
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
    
    # Place waters randomly with minimum distance constraint
    positions = []
    min_dist = 0.35  # nm - reasonable minimum O-O distance
    
    for i in range(n_waters):
        placed = False
        attempts = 0
        while not placed and attempts < 100:
            x = np.random.random() * box_size
            y = np.random.random() * box_size
            z = np.random.random() * box_size
            
            # Check distance to existing waters
            too_close = False
            for pos in positions:
                dx = x - pos[0]
                dy = y - pos[1]
                dz = z - pos[2]
                # Apply PBC
                dx -= box_size * round(dx / box_size)
                dy -= box_size * round(dy / box_size)
                dz -= box_size * round(dz / box_size)
                dist = np.sqrt(dx*dx + dy*dy + dz*dz)
                if dist < min_dist:
                    too_close = True
                    break
            
            if not too_close:
                positions.append([x, y, z])
                placed = True
            attempts += 1
        
        if not placed:
            # Fallback to grid placement
            print(f"Warning: Could not place water {i} randomly, using grid")
            n_dim = int(n_waters**(1/3) + 0.5)
            spacing = box_size / n_dim
            ix = i % n_dim
            iy = (i // n_dim) % n_dim
            iz = i // (n_dim * n_dim)
            x = (ix + 0.5) * spacing
            y = (iy + 0.5) * spacing
            z = (iz + 0.5) * spacing
            positions.append([x, y, z])
        
        # Random rotation
        theta = np.random.random() * 2 * np.pi
        phi = np.random.random() * np.pi
        psi = np.random.random() * 2 * np.pi
        
        # Rotation matrices
        Rz = np.array([[np.cos(theta), -np.sin(theta), 0],
                       [np.sin(theta), np.cos(theta), 0],
                       [0, 0, 1]])
        Ry = np.array([[np.cos(phi), 0, np.sin(phi)],
                       [0, 1, 0],
                       [-np.sin(phi), 0, np.cos(phi)]])
        Rx = np.array([[1, 0, 0],
                       [0, np.cos(psi), -np.sin(psi)],
                       [0, np.sin(psi), np.cos(psi)]])
        R = Rz @ Ry @ Rx
        
        # Water geometry (before rotation)
        water_coords = np.array([
            [0, 0, 0],  # O
            [0.09572, 0, 0],  # H1
            [-0.02399, 0.09277, 0]  # H2
        ])
        
        # Apply rotation
        rotated = water_coords @ R.T
        
        # Oxygen
        o = pygcmc.MCAtom()
        o.x = x + rotated[0, 0]
        o.y = y + rotated[0, 1]
        o.z = z + rotated[0, 2]
        o.charge = qO_core
        o.type = 0
        atoms.append(o)
        
        # Drude - with small offset
        d = pygcmc.MCAtom()
        d.x = o.x + 0.001 * (np.random.random() - 0.5)
        d.y = o.y + 0.001 * (np.random.random() - 0.5)
        d.z = o.z + 0.001 * (np.random.random() - 0.5)
        d.charge = qD
        d.type = 1
        atoms.append(d)
        
        # H1
        h1 = pygcmc.MCAtom()
        h1.x = x + rotated[1, 0]
        h1.y = y + rotated[1, 1]
        h1.z = z + rotated[1, 2]
        h1.charge = qH
        h1.type = 2
        atoms.append(h1)
        
        # H2
        h2 = pygcmc.MCAtom()
        h2.x = x + rotated[2, 0]
        h2.y = y + rotated[2, 1]
        h2.z = z + rotated[2, 2]
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
        res.atomStart = i * 5
        res.atomCount = 5
        res.active = True
        res.type = 0
        residues.append(res)
    
    state.atoms = atoms
    state.activeAtomCount = len(atoms)
    state.residues = residues
    state.activeResidueCount = n_waters
    
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
    
    return state, drude_force, box_size

print("=" * 70)
print("OPTIMIZED WATER BOX PERFORMANCE TEST")
print("=" * 70)

# Test different system sizes with proper density
sizes = [8, 27, 64, 125]

print("\n1. ENERGY COMPARISON")
print("-" * 60)
print(f"{'N waters':>8} | {'Box (nm)':>8} | {'Energy (kJ/mol)':>15} | {'Per molecule':>12} | {'vs Old':>10}")
print("-" * 60)

# Old energies for comparison
old_energies = {8: -6846, 27: 35803, 64: 21573, 125: 735304}

for n_waters in sizes:
    state, drude_force, box_size = create_optimized_water_box(n_waters)
    
    # Use reasonable tolerance
    params = pygcmc.DrudeSCFParams()
    params.tolerance = 10.0  # Looser for faster convergence
    params.maxIterations = 100
    params.maxDrudeDistance = 0.02
    drude_force.setSCFParameters(params)
    
    energy = drude_force.calculateEnergySCF(state)
    old_energy = old_energies.get(n_waters, 0)
    
    print(f"{n_waters:8d} | {box_size:8.2f} | {energy:15.2f} | {energy/n_waters:12.2f} | {old_energy/n_waters:10.0f}")

# Performance test with optimized configuration
print("\n\n2. PERFORMANCE TEST WITH OPTIMIZED CONFIGURATION")
print("-" * 60)

test_sizes = [(2, 8), (3, 27), (4, 64), (5, 125)]
ref_times = {8: 0.003, 27: 0.031, 64: 0.159, 125: 0.643}

print(f"{'N waters':>8} | {'ms/step':>10} | {'steps/sec':>10} | {'vs Non-Drude':>12}")
print("-" * 60)

for _, n_waters in test_sizes[:3]:  # Test smaller systems
    state, drude_force, _ = create_optimized_water_box(n_waters)
    
    # Use default tolerance (1.0)
    # No need to set parameters - defaults are used
    
    # Warmup
    for _ in range(5):
        drude_force.calculateEnergySCF(state)
    
    # Test
    n_steps = min(50, 400 // n_waters)
    times = []
    
    for _ in range(n_steps):
        # Move one water slightly
        mol = np.random.randint(0, n_waters)
        start_atom = mol * 5
        dx = (np.random.random() - 0.5) * 0.001
        dy = (np.random.random() - 0.5) * 0.001
        dz = (np.random.random() - 0.5) * 0.001
        
        for i in [0, 2, 3, 4]:  # Move O, H1, H2, M (not Drude)
            state.atoms[start_atom + i].x += dx
            state.atoms[start_atom + i].y += dy
            state.atoms[start_atom + i].z += dz
        
        t0 = time.time()
        energy = drude_force.calculateEnergySCF(state)
        t1 = time.time()
        times.append(t1 - t0)
    
    avg_time = np.mean(times) * 1000
    steps_per_sec = 1000 / avg_time
    vs_nondrude = avg_time / ref_times[n_waters]
    
    print(f"{n_waters:8d} | {avg_time:10.3f} | {steps_per_sec:10.1f} | {vs_nondrude:12.1f}x")
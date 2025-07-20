#!/usr/bin/env python
"""Check if energies are reasonable"""

import sys
sys.path.append('/home/zhaomt/gcmc/test107/pygcmc_dev/build/modules/bindings')

import pygcmc

# Expected energy per water molecule (rough estimate)
# Single water: 0 kJ/mol (intramolecular excluded)
# Water-water interaction: ~-5 to -10 kJ/mol per pair at 0.3 nm

def check_energy(n_waters, energy):
    """Check if energy is reasonable"""
    # Very rough estimate: each water interacts with ~6 neighbors
    # Each interaction is about -5 kJ/mol
    expected_range = (-30 * n_waters, 10 * n_waters)
    
    print(f"\n{n_waters} waters:")
    print(f"  Total energy: {energy:.2f} kJ/mol")
    print(f"  Per molecule: {energy/n_waters:.2f} kJ/mol")
    print(f"  Expected range: {expected_range[0]:.0f} to {expected_range[1]:.0f} kJ/mol")
    
    if energy < expected_range[0] or energy > expected_range[1]:
        print("  WARNING: Energy outside expected range!")
        print("  This suggests:")
        if energy > expected_range[1]:
            print("  - Possible overlapping atoms or bad initial configuration")
            print("  - Drude particles might be too far from parents")
            print("  - Force field parameters might be incorrect")
        else:
            print("  - Possible double counting of interactions")
    else:
        print("  ✓ Energy seems reasonable")

# Need to recreate the function since import won't work
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
    alpha = 0.0013  # nm^3
    
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
                import numpy as np
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
    ff.numTotalTypes = 4
    ff.numMovementTypes = 4
    
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
    
    return state, drude_force, n_waters

print("=" * 70)
print("ENERGY SANITY CHECK")
print("=" * 70)

sizes = [(2, 8), (3, 27), (4, 64), (5, 125)]

for n_dim, n_waters in sizes:
    state, drude_force, _ = create_drude_water_box(n_dim)
    
    # Calculate energy
    energy = drude_force.calculateEnergySCF(state)
    check_energy(n_waters, energy)

# Also check the standard non-Drude energy
print("\n\nFor comparison, non-Drude water energy:")
print("-" * 50)

# Create simple TIP3P-like water box
state = pygcmc.MCState()
state.info.box = [1.55, 1.55, 1.55]
state.info.cutoff = 0.7
state.info.setTemperature(300.0)

atoms = []
residues = []

# 5x5x5 = 125 waters
for i in range(125):
    x = (i % 5 + 0.5) * 0.31
    y = ((i // 5) % 5 + 0.5) * 0.31
    z = (i // 25 + 0.5) * 0.31
    
    # Simple water (O, H1, H2)
    o = pygcmc.MCAtom()
    o.x, o.y, o.z = x, y, z
    o.charge = -0.834
    o.type = 0
    atoms.append(o)
    
    h1 = pygcmc.MCAtom()
    h1.x = x + 0.09572
    h1.y = y
    h1.z = z
    h1.charge = 0.417
    h1.type = 1
    atoms.append(h1)
    
    h2 = pygcmc.MCAtom()
    h2.x = x - 0.02399
    h2.y = y + 0.09277
    h2.z = z
    h2.charge = 0.417
    h2.type = 1
    atoms.append(h2)
    
    res = pygcmc.MCResidue()
    res.atomStart = i * 3
    res.atomCount = 3
    res.active = True
    res.type = 0
    residues.append(res)

state.atoms = atoms
state.activeAtomCount = len(atoms)
state.residues = residues
state.activeResidueCount = 125

ff = pygcmc.MCForceField()
ff.numTotalTypes = 2
ff.numMovementTypes = 2
ff.ljSigma = [0.315, 0.0, 0.0, 0.0]
ff.ljEps = [0.636, 0.0, 0.0, 0.0]
state.forcefield = ff

# Calculate non-Drude energy
energy_coul = pygcmc.computeIntermolecularCoulombEnergyCutoff(state, 0)
energy_lj = pygcmc.computeIntermolecularLJEnergyCutoff(state, 0)
total_energy = energy_coul + energy_lj

print(f"\n125 TIP3P waters (non-Drude):")
print(f"  Coulomb energy: {energy_coul:.2f} kJ/mol")
print(f"  LJ energy: {energy_lj:.2f} kJ/mol")
print(f"  Total energy: {total_energy:.2f} kJ/mol")
print(f"  Per molecule: {total_energy/125:.2f} kJ/mol")
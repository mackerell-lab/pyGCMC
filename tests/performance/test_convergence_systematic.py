#!/usr/bin/env python
"""Systematic test of convergence issues"""

import sys
import numpy as np
sys.path.append('/home/zhaomt/gcmc/test107/pygcmc_dev/build/modules/bindings')

import pygcmc

def create_single_water():
    """Create a single water molecule"""
    state = pygcmc.MCState()
    state.info.box = [3.0, 3.0, 3.0]
    state.info.cutoff = 1.4
    state.info.setTemperature(300.0)
    
    atoms = []
    
    # PSF parameters
    qO_core = 1.66260
    qD = -1.76260
    qM = -0.957100
    qH = 0.528550
    alpha = 0.0013
    
    # One water at center
    positions = [
        (qO_core, 1.5, 1.5, 1.5, 0),         # O
        (qD, 1.501, 1.5, 1.5, 1),            # D (offset 0.001)
        (qH, 1.59572, 1.5, 1.5, 2),          # H1
        (qH, 1.47601, 1.59277, 1.5, 2),      # H2
        (qM, 1.5, 1.5127, 1.5, 3)            # M
    ]
    
    for charge, x, y, z, atype in positions:
        atom = pygcmc.MCAtom()
        atom.x, atom.y, atom.z = x, y, z
        atom.charge = charge
        atom.type = atype
        atoms.append(atom)
    
    res = pygcmc.MCResidue()
    res.atomStart = 0
    res.atomCount = 5
    res.active = True
    res.type = 0
    
    state.atoms = atoms
    state.activeAtomCount = 5
    state.residues = [res]
    state.activeResidueCount = 1
    
    # Force field
    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 4
    ff.numMovementTypes = 4
    ff.ljSigma = [0.318395, 0, 0, 0] * 4
    ff.ljEps = [0.88257, 0, 0, 0] * 4
    state.forcefield = ff
    
    # Drude force
    drude_force = pygcmc.DrudeForce()
    drude_force.addParticle(1, 0, -1, -1, -1, -1, qD, alpha, 1.0, 1.0)
    
    return state, drude_force

def create_water_dimer():
    """Create two water molecules"""
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
    alpha = 0.0013
    
    # Two waters 0.3 nm apart
    water_positions = [[1.0, 1.5, 1.5], [1.3, 1.5, 1.5]]
    
    mol_id = 0
    for pos in water_positions:
        x, y, z = pos
        
        positions = [
            (qO_core, x, y, z, 0),
            (qD, x+0.001, y, z, 1),
            (qH, x+0.09572, y, z, 2),
            (qH, x-0.02399, y+0.09277, z, 2),
            (qM, x, y+0.0127, z, 3)
        ]
        
        for charge, px, py, pz, atype in positions:
            atom = pygcmc.MCAtom()
            atom.x, atom.y, atom.z = px, py, pz
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
    state.activeResidueCount = 2
    
    # Force field
    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 4
    ff.numMovementTypes = 4
    ff.ljSigma = [0.318395, 0, 0, 0] * 4
    ff.ljEps = [0.88257, 0, 0, 0] * 4
    state.forcefield = ff
    
    # Drude force
    drude_force = pygcmc.DrudeForce()
    drude_force.addParticle(1, 0, -1, -1, -1, -1, qD, alpha, 1.0, 1.0)
    drude_force.addParticle(6, 5, -1, -1, -1, -1, qD, alpha, 1.0, 1.0)
    
    return state, drude_force

print("=" * 70)
print("SYSTEMATIC CONVERGENCE TEST")
print("=" * 70)

# Test 1: Single water molecule
print("\n1. SINGLE WATER MOLECULE")
print("-" * 40)
state, drude_force = create_single_water()

for i in range(5):
    print(f"\nRun {i+1}:")
    energy = drude_force.calculateEnergySCF(state)
    print(f"  Energy = {energy:.6f} kJ/mol")
    
    # Perturb slightly
    state.atoms[0].x += 0.0001

# Test 2: Water dimer
print("\n\n2. WATER DIMER")
print("-" * 40)
state, drude_force = create_water_dimer()

for i in range(5):
    print(f"\nRun {i+1}:")
    energy = drude_force.calculateEnergySCF(state)
    print(f"  Energy = {energy:.6f} kJ/mol")
    
    # Perturb first water
    state.atoms[0].x += 0.0001

# Test 3: Effect of Drude initial position
print("\n\n3. DRUDE INITIAL POSITION EFFECT")
print("-" * 40)
state, drude_force = create_single_water()

offsets = [0.0, 0.0001, 0.001, 0.01, 0.1]

for offset in offsets:
    # Reset Drude position relative to parent
    state.atoms[1].x = state.atoms[0].x + offset
    state.atoms[1].y = state.atoms[0].y
    state.atoms[1].z = state.atoms[0].z
    
    print(f"\nDrude offset = {offset} nm:")
    energy = drude_force.calculateEnergySCF(state)
    print(f"  Energy = {energy:.6f} kJ/mol")

# Test 4: Different tolerances
print("\n\n4. TOLERANCE EFFECT")
print("-" * 40)
state, drude_force = create_water_dimer()

tolerances = [100.0, 10.0, 1.0, 0.1, 0.01]

for tol in tolerances:
    params = pygcmc.DrudeSCFParams()
    params.tolerance = tol
    params.maxIterations = 50
    params.maxDrudeDistance = 0.02
    drude_force.setSCFParameters(params)
    
    print(f"\nTolerance = {tol}:")
    energy = drude_force.calculateEnergySCF(state)
    print(f"  Energy = {energy:.6f} kJ/mol")
#!/usr/bin/env python3
"""
Minimal test of Force Balance Predictor (FBP) algorithm
"""

import pygcmc
import time

# Create minimal water system
atoms = []
residues = []

# Create 10 water molecules
for i in range(10):
    base_x = i * 0.4
    
    # Oxygen
    atom = pygcmc.MCAtom()
    atom.x = base_x
    atom.y = 0.0
    atom.z = 0.0
    atom.charge = 1.71636
    atom.type = 0
    atoms.append(atom)
    
    # Drude
    drude = pygcmc.MCAtom()
    drude.x = base_x
    drude.y = 0.0
    drude.z = 0.0
    drude.charge = -1.71636
    drude.type = 1
    atoms.append(drude)
    
    # H1
    h1 = pygcmc.MCAtom()
    h1.x = base_x + 0.09572
    h1.y = 0.0
    h1.z = 0.0
    h1.charge = 0.55733
    h1.type = 2
    atoms.append(h1)
    
    # H2
    h2 = pygcmc.MCAtom()
    h2.x = base_x - 0.04786
    h2.y = 0.08288
    h2.z = 0.0
    h2.charge = 0.55733
    h2.type = 2
    atoms.append(h2)
    
    # M-site
    m = pygcmc.MCAtom()
    m.x = base_x
    m.y = -0.024034
    m.z = 0.0
    m.charge = -1.11466
    m.type = 3
    atoms.append(m)
    
    # Residue
    res = pygcmc.MCResidue()
    res.atomStart = 5 * i
    res.atomCount = 5
    res.active = True
    res.type = 0
    residues.append(res)

# Create state
state = pygcmc.MCState()
state.atoms = atoms
state.residues = residues
state.activeAtomCount = len(atoms)
state.activeResidueCount = len(residues)

# Box info
import numpy as np
state.info.box = np.array([10.0, 10.0, 10.0])
state.info.cutoff = 4.5

# Force field
state.forcefield.numTotalTypes = 4
state.forcefield.numMovementTypes = 4
state.forcefield.ljSigma = [0.318395, 0.0, 0.0, 0.0]
state.forcefield.ljEps = [0.88257, 0.0, 0.0, 0.0]

print("Testing FBP Algorithm")
print("=" * 40)

# Test algorithms
algorithms = [
    (pygcmc.DrudeAlgorithm.SCF, "SCF"),
    (pygcmc.DrudeAlgorithm.FBP, "FBP")
]

# Drude parameters
charge = -1.71636
polarizability = 1.71636**2 * 138.935456 / 418400.0

for algo, name in algorithms:
    print(f"\nTesting {name}...")
    
    # Create force
    force = pygcmc.DrudeForce()
    
    # Add particles
    for i in range(10):
        force.addParticle(
            drudeIndex=5*i + 1,
            parentIndex=5*i,
            aniso1Index=-1, aniso2Index=-1,
            aniso3Index=-1, aniso4Index=-1,
            charge=charge,
            polarizability=polarizability,
            aniso12=1.0, aniso34=1.0
        )
    
    # Add screened pairs
    for i in range(10):
        for j in range(i+1, 10):
            force.addScreenedPair(i, j, 1.3)
    
    # Set parameters
    params = pygcmc.DrudeSCFParams()
    params.tolerance = 1.0
    params.maxIterations = 50
    force.setSCFParameters(params)
    force.setAlgorithm(algo)
    
    # Reset Drude positions
    for i in range(10):
        state.atoms[5*i + 1].x = state.atoms[5*i].x
        state.atoms[5*i + 1].y = state.atoms[5*i].y
        state.atoms[5*i + 1].z = state.atoms[5*i].z
    
    # Calculate
    start = time.time()
    energy = force.calculateEnergySCF(state)
    elapsed = time.time() - start
    
    print(f"  Time: {elapsed*1000:.2f} ms")
    print(f"  Energy: {energy:.2f} kJ/mol")
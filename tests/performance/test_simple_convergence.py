#!/usr/bin/env python
"""Simple test to understand convergence"""

import sys
sys.path.append('/home/zhaomt/gcmc/test107/pygcmc_dev/build/modules/bindings')

import pygcmc

# Create minimal system - just one water
state = pygcmc.MCState()
state.info.box = [3.0, 3.0, 3.0]
state.info.cutoff = 1.4
state.info.setTemperature(300.0)

atoms = []

# One water molecule
positions = [
    (1.66260, 1.5, 1.5, 1.5, 0),      # O
    (-1.76260, 1.501, 1.5, 1.5, 1),   # D (offset 0.001)
    (0.528550, 1.59572, 1.5, 1.5, 2), # H1
    (0.528550, 1.47601, 1.59277, 1.5, 2), # H2
    (-0.957100, 1.5, 1.5127, 1.5, 3)  # M
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
drude_force.addParticle(1, 0, -1, -1, -1, -1, -1.76260, 0.0013, 1.0, 1.0)

print("Testing single water molecule SCF convergence")
print("=" * 50)

# Test 1: Default parameters
print("\nTest 1: Default (tolerance=1.0)")
energy = drude_force.calculateEnergySCF(state)
print(f"Energy = {energy:.6f} kJ/mol")

# Test 2: Looser tolerance
print("\nTest 2: Looser tolerance=10.0")
params = pygcmc.DrudeSCFParams()
params.tolerance = 10.0
params.maxIterations = 50
params.maxDrudeDistance = 0.02
drude_force.setSCFParameters(params)
energy = drude_force.calculateEnergySCF(state)
print(f"Energy = {energy:.6f} kJ/mol")

# Test 3: Very loose tolerance
print("\nTest 3: Very loose tolerance=100.0")
params.tolerance = 100.0
drude_force.setSCFParameters(params)
energy = drude_force.calculateEnergySCF(state)
print(f"Energy = {energy:.6f} kJ/mol")
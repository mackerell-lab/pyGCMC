#!/usr/bin/env python
"""Check what tolerance is actually being used"""

import sys
sys.path.append('/home/zhaomt/gcmc/test107/pygcmc_dev/build/modules/bindings')

import pygcmc

# Create a simple Drude force
drude_force = pygcmc.DrudeForce()

# Check default parameters
params = pygcmc.DrudeSCFParams()
print(f"Default tolerance: {params.tolerance}")
print(f"Default maxIterations: {params.maxIterations}")
print(f"Default maxDrudeDistance: {params.maxDrudeDistance}")

# Also create a water system and check convergence
def create_simple_water():
    state = pygcmc.MCState()
    state.info.box = [3.0, 3.0, 3.0]
    state.info.cutoff = 1.4
    state.info.setTemperature(300.0)
    
    atoms = []
    residues = []
    
    # One water molecule
    positions = [
        (1.66260, 0, 0, 0, 0),      # O
        (-1.76260, 0, 0, 0, 1),     # D
        (0.528550, 0.09572, 0, 0, 2), # H1
        (0.528550, -0.02399, 0.09277, 0, 2), # H2
        (-0.957100, 0, 0.0127, 0, 3)  # M
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
    residues.append(res)
    
    state.atoms = atoms
    state.activeAtomCount = len(atoms)
    state.residues = residues
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
    
    return state, drude_force

print("\nTesting single water molecule...")
state, drude_force = create_simple_water()

# Test with different tolerances
for tol in [10.0, 1.0, 0.1, 0.01, 1e-3, 1e-6]:
    params = pygcmc.DrudeSCFParams()
    params.tolerance = tol
    params.maxIterations = 50
    params.maxDrudeDistance = 0.02
    drude_force.setSCFParameters(params)
    
    print(f"\nTolerance = {tol}:")
    energy = drude_force.calculateEnergySCF(state)
    print(f"  Energy = {energy:.6f} kJ/mol")
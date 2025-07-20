#!/usr/bin/env python
"""Simple test to see debug output"""

import sys
sys.path.append('/home/zhaomt/gcmc/test107/pygcmc_dev/build/modules/bindings')

import pygcmc

def create_drude_water_box(n_per_dim):
    n_waters = n_per_dim ** 3
    box_size = n_per_dim * 0.31
    
    state = pygcmc.MCState()
    state.info.box = [box_size, box_size, box_size]
    state.info.cutoff = min(box_size/2 - 0.01, 1.2)
    state.info.setTemperature(300.0)
    
    atoms = []
    residues = []
    
    qO_core = 1.66260
    qD = -1.76260
    qM = -0.957100
    qH = 0.528550
    alpha = 0.0013
    
    mol_id = 0
    for i in range(n_waters):
        x = (i % n_per_dim + 0.5) * 0.31
        y = ((i // n_per_dim) % n_per_dim + 0.5) * 0.31
        z = (i // (n_per_dim * n_per_dim) + 0.5) * 0.31
        
        # Create 5 atoms per water
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
    state.activeResidueCount = n_waters
    
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
    
    drude_force = pygcmc.DrudeForce()
    
    for i in range(n_waters):
        drude_force.addParticle(i*5+1, i*5, -1, -1, -1, -1, qD, alpha, 1.0, 1.0)
    
    return state, drude_force, n_waters

# Test small system first
print("Testing 8 waters (should converge)...")
state, drude_force, _ = create_drude_water_box(2)
energy = drude_force.calculateEnergySCF(state)
print(f"Energy = {energy:.2f} kJ/mol\n")

# Test larger system
print("Testing 125 waters (may not converge)...")
state, drude_force, _ = create_drude_water_box(5)
energy = drude_force.calculateEnergySCF(state)
print(f"Energy = {energy:.2f} kJ/mol")
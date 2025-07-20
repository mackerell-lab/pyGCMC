#!/usr/bin/env python
"""Debug SCF convergence issues"""

import sys
sys.path.append('/home/zhaomt/gcmc/test107/pygcmc_dev/build/modules/bindings')

import pygcmc
import numpy as np

def create_water_box(n_waters):
    """Create a small water box"""
    n_dim = int(n_waters**(1/3) + 0.5)
    box_size = n_dim * 0.31
    
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
    for i in range(n_waters):
        x = (i % n_dim + 0.5) * 0.31
        y = ((i // n_dim) % n_dim + 0.5) * 0.31
        z = (i // (n_dim * n_dim) + 0.5) * 0.31
        
        # Oxygen
        o = pygcmc.MCAtom()
        o.x, o.y, o.z = x, y, z
        o.charge = qO_core
        o.type = 0
        atoms.append(o)
        
        # Drude
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

print("Testing SCF convergence with different tolerances...\n")

# Test with a small system
state, drude_force = create_water_box(8)

# Test different tolerances
tolerances = [100.0, 10.0, 1.0, 0.1]

for tol in tolerances:
    print(f"\n{'='*60}")
    print(f"Testing with tolerance = {tol} kJ/mol/nm")
    print(f"{'='*60}")
    
    params = pygcmc.DrudeSCFParams()
    params.tolerance = tol
    params.maxIterations = 10  # Fewer iterations to see output
    params.maxDrudeDistance = 0.02
    drude_force.setSCFParameters(params)
    
    print("\nCalculating energy (check stderr for debug output)...")
    energy = drude_force.calculateEnergySCF(state)
    print(f"Energy = {energy:.6f} kJ/mol")
    
    # Move atoms slightly to create different initial conditions
    state.atoms[0].x += 0.001
    state.atoms[2].x += 0.001
    state.atoms[3].x += 0.001
    state.atoms[4].x += 0.001
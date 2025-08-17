#!/usr/bin/env python3
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../build/modules/bindings'))

import pygcmc
from .active_pool import ActivePool

# Simple water molecule 
def create_water():
    atoms = []
    # O
    o = pygcmc.MCAtom()
    o.x, o.y, o.z = 0, 0, 0
    o.charge = -0.834
    o.type = 0
    atoms.append(o)
    # H1
    h1 = pygcmc.MCAtom()
    h1.x, h1.y, h1.z = 0.1, 0, 0
    h1.charge = 0.417
    h1.type = 1
    atoms.append(h1)
    # H2
    h2 = pygcmc.MCAtom()
    h2.x, h2.y, h2.z = 0, 0.1, 0
    h2.charge = 0.417
    h2.type = 1
    atoms.append(h2)
    return atoms

ff = pygcmc.MCForceField()
ff.numTotalTypes = 2
ff.ljSigma = [0.315, 0.0]
ff.ljEps = [0.636, 0.0]

pool = ActivePool(box=[5.0, 5.0, 5.0], cutoff=1.2, forcefield=ff)

print("Inserting 3 waters...")
res0 = pool.insert_molecule(create_water())
print(f"After res0: atoms={len(pool.state.atoms)}, metadata={[m.__dict__ for m in pool.residue_metadata]}")

res1 = pool.insert_molecule(create_water())
print(f"After res1: atoms={len(pool.state.atoms)}, metadata={[m.__dict__ for m in pool.residue_metadata]}")

res2 = pool.insert_molecule(create_water())
print(f"After res2: atoms={len(pool.state.atoms)}, metadata={[m.__dict__ for m in pool.residue_metadata]}")

print("\nDeleting res1...")
pool.delete_residue(res1)
print(f"After delete: metadata={[m.__dict__ for m in pool.residue_metadata]}")

print("\nCompacting...")
compacted = pool.compact(force=True)
print(f"Compacted {compacted} residues")
print(f"After compact: metadata={[m.__dict__ for m in pool.residue_metadata]}")
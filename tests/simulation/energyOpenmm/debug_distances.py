"""
Debug atom distances in test systems
"""

import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pygcmc
from pygcmc import MCState, MCAtom, MCResidue, MCForceField


def calculate_distances(atoms, box):
    """Calculate all pairwise distances with PBC"""
    n = len(atoms)
    distances = []
    
    for i in range(n):
        for j in range(i+1, n):
            dx = atoms[i].x - atoms[j].x
            dy = atoms[i].y - atoms[j].y
            dz = atoms[i].z - atoms[j].z
            
            # Apply minimum image convention
            dx -= box[0] * round(dx / box[0])
            dy -= box[1] * round(dy / box[1])
            dz -= box[2] * round(dz / box[2])
            
            r = math.sqrt(dx*dx + dy*dy + dz*dz)
            distances.append((i, j, r))
            
    return distances


def test_system_from_test_full_nonbonded():
    """Check the system from test_full_nonbonded_comparison.py"""
    
    box_size = 4.0
    cutoff = 1.2
    
    state = MCState()
    state.info.box = [box_size, box_size, box_size]
    state.info.cutoff = cutoff
    
    # Force field with real LJ parameters
    ff = MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.ljEps = [0.996]
    ff.ljSigma = [0.340]
    state.forcefield = ff
    
    atoms = []
    
    # From the modified test
    test_atoms = [
        ([1.0, 2.0, 2.0], 0.5),
        ([3.0, 2.0, 2.0], -0.5),
        ([2.0, 0.5, 2.0], 0.2),
        ([2.0, 3.5, 2.0], -0.2),
        ([1.5, 1.5, 3.0], 0.0),
        ([2.5, 2.5, 1.0], 0.0),
    ]
    
    for pos, charge in test_atoms:
        atom = MCAtom()
        atom.x, atom.y, atom.z = pos
        atom.charge = charge
        atom.type = 0
        atoms.append(atom)
    
    state.atoms = atoms
    
    print("System from test_full_nonbonded_comparison:")
    print(f"Box: {box_size} nm, Cutoff: {cutoff} nm")
    print(f"LJ params: eps={ff.ljEps[0]}, sigma={ff.ljSigma[0]}")
    print("\nAtom positions:")
    for i, (pos, charge) in enumerate(test_atoms):
        print(f"  Atom {i}: {pos}, charge={charge}")
    
    print("\nPairwise distances:")
    distances = calculate_distances(atoms, state.info.box)
    for i, j, r in distances:
        within_cutoff = "✓" if r < cutoff else "✗"
        print(f"  {i}-{j}: {r:.3f} nm {within_cutoff}")
    
    # Count interactions within cutoff
    n_within = sum(1 for _, _, r in distances if r < cutoff)
    print(f"\nTotal pairs within cutoff: {n_within}")


def test_system_from_test_pme_lj_only():
    """Check the system from test_pme_lj_only.py"""
    
    box_size = 3.0
    cutoff = 1.2
    
    print("\n\nSystem from test_pme_lj_only:")
    print(f"Box: {box_size} nm, Cutoff: {cutoff} nm")
    
    # Modified positions
    positions = [
        ([1.0, 1.5, 1.5], 0),
        ([1.6, 1.5, 1.5], 0),
        ([2.2, 1.5, 1.5], 0),
        ([1.5, 0.9, 1.5], 1),
        ([1.5, 2.1, 1.5], 1),
        ([1.5, 1.5, 2.2], 1),
    ]
    
    atoms = []
    for pos, atom_type in positions:
        atom = MCAtom()
        atom.x, atom.y, atom.z = pos
        atom.type = atom_type
        atoms.append(atom)
    
    print("\nAtom positions:")
    for i, (pos, t) in enumerate(positions):
        print(f"  Atom {i}: {pos}, type={t}")
    
    print("\nPairwise distances:")
    distances = calculate_distances(atoms, [box_size, box_size, box_size])
    for i, j, r in distances:
        within_cutoff = "✓" if r < cutoff else "✗"
        print(f"  {i}-{j}: {r:.3f} nm {within_cutoff}")
    
    # Count interactions within cutoff
    n_within = sum(1 for _, _, r in distances if r < cutoff)
    print(f"\nTotal pairs within cutoff: {n_within}")


if __name__ == "__main__":
    test_system_from_test_full_nonbonded()
    test_system_from_test_pme_lj_only()
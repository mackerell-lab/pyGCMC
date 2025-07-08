"""
Debug the test_lj_only_cutoff system
"""

import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pygcmc
from pygcmc import MCState, MCAtom, MCResidue, MCForceField
from pygcmc import computeSystemEnergyCutoffFixed


def create_lj_only_system(n_atoms=6):
    """Recreate the system from test_pme_lj_only.py"""
    
    box_size = 3.0  # nm
    cutoff = 1.2    # nm
    
    state = MCState()
    state.info.box = [box_size, box_size, box_size]
    state.info.cutoff = cutoff
    
    # Force field with real LJ parameters
    ff = MCForceField()
    ff.numTotalTypes = 2  # Two types for testing mixing rules
    ff.numMovementTypes = 2
    # Need to provide full interaction matrix (2x2 = 4 values)
    eps0 = 0.996   # Argon-like
    eps1 = 1.230   # Methane-like
    sigma0 = 0.340
    sigma1 = 0.373
    
    # Full matrix: [00, 01, 10, 11]
    ff.ljEps = [
        eps0,                          # 0-0
        math.sqrt(eps0 * eps1),        # 0-1
        math.sqrt(eps0 * eps1),        # 1-0
        eps1                           # 1-1
    ]
    ff.ljSigma = [
        sigma0,                        # 0-0
        (sigma0 + sigma1) / 2,         # 0-1
        (sigma0 + sigma1) / 2,         # 1-0
        sigma1                         # 1-1
    ]
    state.forcefield = ff
    
    atoms = []
    residues = []
    
    # Create atoms in a regular pattern with safe distances (> 0.4 nm apart)
    positions = [
        ([1.0, 1.5, 1.5], 0),  # Type 0
        ([1.6, 1.5, 1.5], 0),  # Type 0 - increased spacing
        ([2.2, 1.5, 1.5], 0),  # Type 0 - increased spacing
        ([1.5, 0.9, 1.5], 1),  # Type 1 - increased spacing
        ([1.5, 2.1, 1.5], 1),  # Type 1 - increased spacing
        ([1.5, 1.5, 2.2], 1),  # Type 1 - increased spacing
    ]
    
    for i, (pos, atom_type) in enumerate(positions[:n_atoms]):
        atom = MCAtom()
        atom.x, atom.y, atom.z = pos
        atom.charge = 0.0  # NO CHARGE
        atom.type = atom_type
        atoms.append(atom)
        
        res = MCResidue()
        res.active = True
        res.fixed = True
        res.atomStart = i
        res.atomCount = 1
        res.type = atom_type
        residues.append(res)
    
    state.atoms = atoms
    state.activeAtomCount = len(atoms)
    state.residues = residues
    state.activeResidueCount = len(residues)
    
    return state


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
            distances.append((i, j, r, atoms[i].type, atoms[j].type))
            
    return distances


def main():
    state = create_lj_only_system(n_atoms=6)
    
    print("Debug test_lj_only_cutoff system:")
    print(f"Box: {state.info.box[0]} nm")
    print(f"Cutoff: {state.info.cutoff} nm")
    print(f"\nLJ parameters:")
    print(f"  Type 0-0: eps={state.forcefield.ljEps[0]:.3f}, sigma={state.forcefield.ljSigma[0]:.3f}")
    print(f"  Type 0-1: eps={state.forcefield.ljEps[1]:.3f}, sigma={state.forcefield.ljSigma[1]:.3f}")
    print(f"  Type 1-0: eps={state.forcefield.ljEps[2]:.3f}, sigma={state.forcefield.ljSigma[2]:.3f}")
    print(f"  Type 1-1: eps={state.forcefield.ljEps[3]:.3f}, sigma={state.forcefield.ljSigma[3]:.3f}")
    
    print("\nAtom positions:")
    for i, atom in enumerate(state.atoms):
        print(f"  Atom {i}: ({atom.x:.1f}, {atom.y:.1f}, {atom.z:.1f}), type={atom.type}")
    
    print("\nPairwise distances:")
    distances = calculate_distances(state.atoms, state.info.box)
    n_within_cutoff = 0
    for i, j, r, ti, tj in distances:
        within = "✓" if r < state.info.cutoff else "✗"
        print(f"  {i}-{j}: {r:.3f} nm, types=({ti},{tj}) {within}")
        if r < state.info.cutoff:
            n_within_cutoff += 1
    
    print(f"\nTotal pairs within cutoff: {n_within_cutoff}")
    
    # Calculate energy using fixed function
    elec, vdw, total = computeSystemEnergyCutoffFixed(state)
    print(f"\nPyGCMC Fixed energy: {vdw:.6f} kJ/mol")
    
    # Check for any pairs that might be too close
    min_dist = min(r for _, _, r, _, _ in distances)
    print(f"\nMinimum distance: {min_dist:.3f} nm")
    if min_dist < 0.1:
        print("⚠️  WARNING: Some atoms are very close!")


if __name__ == "__main__":
    main()
"""
Debug real space pair counting in PME
"""

import pygcmc
from pygcmc import MCState, MCAtom, MCResidue, MCForceField
from pme_medium_complexity import create_medium_complexity_system
import math

# Create the system
state = create_medium_complexity_system()

print("Debug Real Space Pair Counting")
print("=" * 60)
print(f"Total atoms: {state.activeAtomCount}")
print(f"Total residues: {state.activeResidueCount}")

# Count pairs manually - the way PME real space does it
pair_count = 0
close_pairs = 0
cutoff = 1.2

# Count by residue pairs (as PME does)
for r1 in range(state.activeResidueCount):
    if not state.residues[r1].active:
        continue
    
    for r2 in range(r1 + 1, state.activeResidueCount):
        if not state.residues[r2].active:
            continue
        
        # Count atom pairs within these residues
        for i in range(state.residues[r1].atomStart, 
                      state.residues[r1].atomStart + state.residues[r1].atomCount):
            if i >= state.activeAtomCount:
                continue
                
            for j in range(state.residues[r2].atomStart,
                          state.residues[r2].atomStart + state.residues[r2].atomCount):
                if j >= state.activeAtomCount:
                    continue
                
                # Calculate distance
                dx = state.atoms[i].x - state.atoms[j].x
                dy = state.atoms[i].y - state.atoms[j].y
                dz = state.atoms[i].z - state.atoms[j].z
                
                # Apply PBC
                box = state.info.box
                dx -= box[0] * round(dx / box[0])
                dy -= box[1] * round(dy / box[1])
                dz -= box[2] * round(dz / box[2])
                
                r = math.sqrt(dx*dx + dy*dy + dz*dz)
                
                if r <= cutoff:
                    pair_count += 1
                    if r < 0.5:
                        close_pairs += 1

print(f"\nPairs counted by residue loop (PME method):")
print(f"  Total pairs within cutoff: {pair_count}")
print(f"  Close pairs (< 0.5 nm): {close_pairs}")

# Count all pairs (direct method)
all_pairs = 0
all_close = 0
for i in range(state.activeAtomCount):
    for j in range(i+1, state.activeAtomCount):
        dx = state.atoms[i].x - state.atoms[j].x
        dy = state.atoms[i].y - state.atoms[j].y
        dz = state.atoms[i].z - state.atoms[j].z
        
        # Apply PBC
        box = state.info.box
        dx -= box[0] * round(dx / box[0])
        dy -= box[1] * round(dy / box[1])
        dz -= box[2] * round(dz / box[2])
        
        r = math.sqrt(dx*dx + dy*dy + dz*dz)
        
        if r <= cutoff:
            all_pairs += 1
            if r < 0.5:
                all_close += 1

print(f"\nPairs counted by all-pairs method:")
print(f"  Total pairs within cutoff: {all_pairs}")
print(f"  Close pairs (< 0.5 nm): {all_close}")

# Check residue structure
print(f"\nResidue structure:")
for i in range(min(5, state.activeResidueCount)):
    res = state.residues[i]
    print(f"  Residue {i}: atoms {res.atomStart} to {res.atomStart + res.atomCount - 1} (count={res.atomCount})")

# Check if we're missing intra-residue pairs
print(f"\nChecking for intra-residue pairs:")
intra_pairs = 0
for r in range(state.activeResidueCount):
    if not state.residues[r].active:
        continue
    
    res = state.residues[r]
    n_atoms = res.atomCount
    n_pairs = n_atoms * (n_atoms - 1) // 2
    if n_pairs > 0:
        print(f"  Residue {r}: {n_atoms} atoms, {n_pairs} possible intra-residue pairs")
        intra_pairs += n_pairs

print(f"Total possible intra-residue pairs: {intra_pairs}")
print(f"\nMissing pairs: {all_pairs - pair_count} (these are likely intra-residue pairs)")
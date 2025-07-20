#!/usr/bin/env python
"""Debug LJ calculation"""

import sys
sys.path.append('/home/zhaomt/gcmc/test107/pygcmc_dev/build/modules/bindings')

import pygcmc

def test_lj_debug():
    print("=== Debug LJ Calculation ===\n")
    
    state = pygcmc.MCState()
    state.info.cutoff = 2.0
    
    # Force field
    state.forcefield.numTotalTypes = 4
    state.forcefield.numMovementTypes = 4
    
    # LJ parameters
    sigma_O = 0.318395  # nm
    eps_O = 0.21094 * 4.184  # kJ/mol
    
    print("Setting up LJ parameters:")
    print(f"  sigma_O = {sigma_O} nm")
    print(f"  eps_O = {eps_O} kJ/mol")
    print()
    
    # Create lists
    ljSigma = [0.0] * 16
    ljEps = [0.0] * 16
    
    # Set O-O interaction (type 0 with type 0)
    idx = 0 * 4 + 0
    print(f"Setting index {idx}: sigma={sigma_O}, eps={eps_O}")
    ljSigma[idx] = sigma_O
    ljEps[idx] = eps_O
    
    # Assign to forcefield
    state.forcefield.ljSigma = ljSigma
    state.forcefield.ljEps = ljEps
    
    # Verify it was set
    print(f"After setting: ljEps[{idx}] = {state.forcefield.ljEps[idx]}")
    
    print("\nLJ matrix (epsilon values):")
    for i in range(4):
        row = []
        for j in range(4):
            idx = i * 4 + j
            row.append(f"{state.forcefield.ljEps[idx]:.3f}")
        print(f"  Type {i}: {' '.join(row)}")
    print()
    
    # Create simple system - two O atoms
    atoms = []
    
    # O1
    o1 = pygcmc.MCAtom()
    o1.x = 0.0
    o1.y = 0.0
    o1.z = 0.0
    o1.charge = 0.0  # No charge to isolate LJ
    o1.type = 0
    atoms.append(o1)
    
    # O2
    o2 = pygcmc.MCAtom()
    o2.x = 0.3  # 3 Å
    o2.y = 0.0
    o2.z = 0.0
    o2.charge = 0.0
    o2.type = 0
    atoms.append(o2)
    
    state.atoms = atoms
    state.activeAtomCount = 2
    
    # Two residues
    res1 = pygcmc.MCResidue()
    res1.atomStart = 0
    res1.atomCount = 1
    res1.active = True
    res1.type = 0
    residues = [res1]
    
    res2 = pygcmc.MCResidue()
    res2.atomStart = 1
    res2.atomCount = 1
    res2.active = True
    res2.type = 0
    residues.append(res2)
    
    state.residues = residues
    state.activeResidueCount = 2
    
    print("System setup:")
    print(f"  Atom 1: type={o1.type}, pos=({o1.x}, {o1.y}, {o1.z})")
    print(f"  Atom 2: type={o2.type}, pos=({o2.x}, {o2.y}, {o2.z})")
    print(f"  Distance: 0.3 nm = 3.0 Å")
    print()
    
    # Calculate energy
    pygcmc.computeSystemEnergyCutoff(state)
    
    print("Energy results:")
    print(f"  Residue 1 VDW: {state.residues[0].energy_vdw:.3f} kJ/mol")
    print(f"  Residue 2 VDW: {state.residues[1].energy_vdw:.3f} kJ/mol")
    print(f"  Total VDW: {state.residues[0].energy_vdw + state.residues[1].energy_vdw:.3f} kJ/mol")
    print(f"  Total/2: {(state.residues[0].energy_vdw + state.residues[1].energy_vdw)/2:.3f} kJ/mol")
    
    # Expected LJ
    r = 0.3
    sigma_r = sigma_O / r
    lj_expected = 4 * eps_O * (sigma_r**12 - sigma_r**6)
    print(f"\n  Expected LJ: {lj_expected:.3f} kJ/mol")

if __name__ == "__main__":
    test_lj_debug()
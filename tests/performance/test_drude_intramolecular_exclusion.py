#!/usr/bin/env python
"""Test Drude with intramolecular exclusion"""

import sys
import numpy as np
sys.path.append('/home/zhaomt/gcmc/test107/pygcmc_dev/build/modules/bindings')

import pygcmc

def test_intramolecular_exclusion():
    print("=== Testing Drude with Intramolecular Exclusion ===\n")
    
    # Initialize
    pygcmc.initializeDrudeForce()
    
    # Create state
    state = pygcmc.MCState()
    state.info.cutoff = 1.0
    state.info.box = [5.0, 5.0, 5.0]
    
    # SWM4-NDP parameters
    qO = 1.71636
    qD = -1.71636
    qH = 0.55733
    qM = -1.11466
    
    # Create single water molecule
    atoms = []
    
    # Oxygen
    o = pygcmc.MCAtom()
    o.x = 0.0
    o.y = 0.0
    o.z = 0.0
    o.charge = qO
    o.type = 0
    atoms.append(o)
    
    # Drude on oxygen
    d = pygcmc.MCAtom()
    d.x = 0.0
    d.y = 0.0
    d.z = 0.0
    d.charge = qD
    d.type = 1
    atoms.append(d)
    
    # Hydrogen 1
    h1 = pygcmc.MCAtom()
    h1.x = 0.09572
    h1.y = 0.0
    h1.z = 0.0
    h1.charge = qH
    h1.type = 2
    atoms.append(h1)
    
    # Hydrogen 2
    h2 = pygcmc.MCAtom()
    h2.x = -0.023999
    h2.y = 0.092663
    h2.z = 0.0
    h2.charge = qH
    h2.type = 2
    atoms.append(h2)
    
    # M-site
    m = pygcmc.MCAtom()
    m.x = 0.00793
    m.y = 0.00986
    m.z = 0.0
    m.charge = qM
    m.type = 3
    atoms.append(m)
    
    state.atoms = atoms
    state.activeAtomCount = 5
    
    # Create residue to mark all atoms as same molecule
    res = pygcmc.MCResidue()
    res.atomStart = 0
    res.atomCount = 5
    res.active = True
    res.fixed = False
    res.type = 0
    state.residues = [res]
    state.activeResidueCount = 1
    
    # Force field
    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 4
    state.forcefield = ff
    
    # Calculate polarizability
    k_kj_nm2 = 418400.0  # 1000 kcal/mol/Å²
    polarizability = 138.935456 * qD * qD / k_kj_nm2
    
    # Add Drude particle
    pygcmc.addDrudeParticle(
        drudeIndex=1,
        parentIndex=0,
        charge=qD,
        polarizability=polarizability
    )
    
    # Set SCF parameters
    params = pygcmc.DrudeSCFParams()
    params.tolerance = 1e-8
    params.maxIterations = 100
    params.dampingFactor = 0.5
    params.maxDrudeDistance = 0.02
    pygcmc.setDrudeSCFParameters(params)
    
    print("Test 1: Single water molecule (all intramolecular)\n")
    
    energy, _ = pygcmc.computeSystemEnergyDrude(state)
    
    # Check Drude position
    dx = state.atoms[1].x - state.atoms[0].x
    dy = state.atoms[1].y - state.atoms[0].y
    dz = state.atoms[1].z - state.atoms[0].z
    dist = np.sqrt(dx*dx + dy*dy + dz*dz)
    
    print(f"Drude displacement: {dist*1000:.6f} pm")
    print(f"Energy: {energy:.2f} kJ/mol")
    print(f"Should be ~0 energy (only harmonic term, no Coulomb)")
    
    # Test 2: Add external charge (intermolecular)
    print("\n\nTest 2: Add external charge\n")
    
    # Add external charge to create intermolecular interaction
    ext = pygcmc.MCAtom()
    ext.x = 0.3  # 3 Å away
    ext.y = 0.0
    ext.z = 0.0
    ext.charge = 1.0
    ext.type = 4
    state.atoms.append(ext)
    state.activeAtomCount = 6
    
    # Create separate residue for external charge
    res2 = pygcmc.MCResidue()
    res2.atomStart = 5
    res2.atomCount = 1
    res2.active = True
    res2.fixed = False
    res2.type = 1
    state.residues.append(res2)
    state.activeResidueCount = 2
    
    # Reset Drude
    state.atoms[1].x = 0.0
    state.atoms[1].y = 0.0
    state.atoms[1].z = 0.0
    
    energy2, _ = pygcmc.computeSystemEnergyDrude(state)
    
    dx = state.atoms[1].x - state.atoms[0].x
    dist2 = abs(dx)
    
    print(f"Drude displacement: {dist2*1000:.6f} pm")
    print(f"Energy: {energy2:.2f} kJ/mol")
    print(f"Should have significant energy from intermolecular interactions")
    
    # Test 3: Two water molecules
    print("\n\nTest 3: Two water molecules\n")
    
    # Create second water
    atoms2 = []
    for i in range(5):
        a = pygcmc.MCAtom()
        a.x = state.atoms[i].x + 0.3  # 3 Å away
        a.y = state.atoms[i].y
        a.z = state.atoms[i].z
        a.charge = state.atoms[i].charge
        a.type = state.atoms[i].type
        atoms2.append(a)
    
    # Remove external charge, add second water
    state.atoms = state.atoms[:5] + atoms2
    state.activeAtomCount = 10
    
    # Update residues
    res2.atomStart = 5
    res2.atomCount = 5
    res2.type = 0
    state.residues[1] = res2
    
    # Add second Drude
    pygcmc.addDrudeParticle(
        drudeIndex=6,
        parentIndex=5,
        charge=qD,
        polarizability=polarizability
    )
    
    # Reset both Drudes
    state.atoms[1].x = 0.0
    state.atoms[6].x = 0.3
    
    energy3, _ = pygcmc.computeSystemEnergyDrude(state)
    
    # Check both Drude displacements
    dx1 = state.atoms[1].x - state.atoms[0].x
    dx2 = state.atoms[6].x - state.atoms[5].x
    
    print(f"Water 1 Drude displacement: {dx1*1000:.6f} pm")
    print(f"Water 2 Drude displacement: {dx2*1000:.6f} pm")
    print(f"Energy: {energy3:.2f} kJ/mol")
    print(f"Both Drudes should move due to intermolecular polarization")

if __name__ == "__main__":
    test_intramolecular_exclusion()
#!/usr/bin/env python
"""Test SWM4-NDP water with fixed Drude implementation"""

import sys
import numpy as np
sys.path.append('/home/zhaomt/gcmc/test107/pygcmc_dev/build/modules/bindings')

import pygcmc

def test_swm4_fixed():
    print("=== Testing SWM4-NDP Water with Fixed Drude ===\n")
    
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
    
    # Create water molecule
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
    
    # Force field
    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 4
    state.forcefield = ff
    
    # Calculate polarizability from k = 1000 kcal/mol/Å²
    k_kcal = 1000.0  # kcal/mol/Å²
    k_kj_nm2 = k_kcal * 4.184 * 100.0  # kJ/mol/nm²
    polarizability = 138.935456 * qD * qD / k_kj_nm2
    
    print(f"Drude parameters:")
    print(f"  Charge: {qD}")
    print(f"  k: {k_kcal} kcal/mol/Å² = {k_kj_nm2} kJ/mol/nm²")
    print(f"  Polarizability: {polarizability:.6e} nm³")
    print(f"  α in Å³: {polarizability*1000:.6f}")
    
    # Add Drude particle
    pygcmc.addDrudeParticle(
        drudeIndex=1,
        parentIndex=0,
        charge=qD,
        polarizability=polarizability
    )
    
    # Set SCF parameters with different hard wall distances
    print("\n\nSingle water energy with different hard walls:\n")
    
    for max_dist in [0.001, 0.002, 0.005, 0.01, 0.02]:
        params = pygcmc.DrudeSCFParams()
        params.tolerance = 1e-8
        params.maxIterations = 100
        params.dampingFactor = 0.5
        params.maxDrudeDistance = max_dist
        pygcmc.setDrudeSCFParameters(params)
        
        # Reset Drude
        state.atoms[1].x = 0.0
        state.atoms[1].y = 0.0
        state.atoms[1].z = 0.0
        
        energy, _ = pygcmc.computeSystemEnergyDrude(state)
        
        # Check displacement
        dx = state.atoms[1].x - state.atoms[0].x
        dy = state.atoms[1].y - state.atoms[0].y
        dz = state.atoms[1].z - state.atoms[0].z
        dist = np.sqrt(dx*dx + dy*dy + dz*dz)
        
        print(f"Max distance {max_dist*1000:.1f} pm:")
        print(f"  Drude displacement: {dist*1000:.3f} pm")
        print(f"  Energy: {energy:.2f} kJ/mol")
        
    # Test with external field
    print("\n\nWith external field (test charge at 5 Å):\n")
    
    # Add test charge
    test = pygcmc.MCAtom()
    test.x = 0.5  # 5 Å away
    test.y = 0.0
    test.z = 0.0
    test.charge = 1.0
    test.type = 4
    state.atoms.append(test)
    state.activeAtomCount = 6
    
    params.maxDrudeDistance = 0.02  # OpenMM default
    pygcmc.setDrudeSCFParameters(params)
    
    # Reset Drude
    state.atoms[1].x = 0.0
    state.atoms[1].y = 0.0
    state.atoms[1].z = 0.0
    
    energy, _ = pygcmc.computeSystemEnergyDrude(state)
    
    dx = state.atoms[1].x - state.atoms[0].x
    dist = abs(dx)
    
    print(f"Drude displacement: {dist*1000:.3f} pm ({dist*10:.4f} Å)")
    print(f"Energy: {energy:.2f} kJ/mol")
    
    # Calculate expected displacement
    # F = k*q1*q2/r² = 138.935456 * (-1.71636) * 1.0 / 0.5²
    F_coulomb = 138.935456 * 1.71636 * 1.0 / (0.5 * 0.5)
    # For equilibrium: F_coulomb = k_spring * x
    expected_disp = F_coulomb / k_kj_nm2
    print(f"\nExpected displacement (no hard wall): {expected_disp*1000:.3f} pm")
    print(f"Ratio actual/expected: {dist/expected_disp:.3f}")

if __name__ == "__main__":
    test_swm4_fixed()
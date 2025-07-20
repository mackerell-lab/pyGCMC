#!/usr/bin/env python
"""Diagnose Drude force calculation issues"""

import sys
import numpy as np
sys.path.append('/home/zhaomt/gcmc/test107/pygcmc_dev/build/modules/bindings')

import pygcmc

def diagnose_forces():
    print("=== Diagnosing Drude Force Calculation ===\n")
    
    # Initialize
    pygcmc.initializeDrudeForce()
    
    # Create minimal system: just O, D, and one H
    state = pygcmc.MCState()
    state.info.cutoff = 2.0
    state.info.box = [10.0, 10.0, 10.0]
    
    atoms = []
    
    # Oxygen (parent)
    o = pygcmc.MCAtom()
    o.x = 0.0
    o.y = 0.0
    o.z = 0.0
    o.charge = 1.71636
    o.type = 0
    atoms.append(o)
    
    # Drude
    d = pygcmc.MCAtom()
    d.x = 0.0
    d.y = 0.0
    d.z = 0.0
    d.charge = -1.71636
    d.type = 1
    atoms.append(d)
    
    # Hydrogen (to create field)
    h = pygcmc.MCAtom()
    h.x = 0.1  # 1 Å away
    h.y = 0.0
    h.z = 0.0
    h.charge = 0.55733
    h.type = 2
    atoms.append(h)
    
    state.atoms = atoms
    state.activeAtomCount = 3
    
    # Force field
    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 3
    state.forcefield = ff
    
    # Set up Drude with known parameters
    k_kj_nm2 = 418400.0  # 1000 kcal/mol/Å²
    polarizability = 138.935456 * 1.71636 * 1.71636 / k_kj_nm2
    
    print(f"System setup:")
    print(f"  O at origin, charge = +1.71636")
    print(f"  D at origin, charge = -1.71636")
    print(f"  H at (1Å, 0, 0), charge = +0.55733")
    print(f"  k = {k_kj_nm2} kJ/mol/nm²")
    print(f"  α = {polarizability:.6e} nm³\n")
    
    pygcmc.addDrudeParticle(
        drudeIndex=1,
        parentIndex=0,
        charge=-1.71636,
        polarizability=polarizability
    )
    
    # Test 1: Large hard wall, let SCF find equilibrium
    print("Test 1: Find natural equilibrium (large hard wall)\n")
    
    params = pygcmc.DrudeSCFParams()
    params.tolerance = 1e-10
    params.maxIterations = 200
    params.dampingFactor = 0.3
    params.maxDrudeDistance = 0.1  # 1 Å - very large
    pygcmc.setDrudeSCFParameters(params)
    
    # Calculate
    energy1, _ = pygcmc.computeSystemEnergyDrude(state)
    
    dx = state.atoms[1].x
    print(f"Drude position: x = {dx*1000:.6f} pm")
    print(f"Total energy: {energy1:.2f} kJ/mol")
    
    # Calculate forces manually
    # Force on D from H: F = k*qD*qH/r²
    r_DH = 0.1 - dx  # H is at 0.1 nm
    F_DH = 138.935456 * (-1.71636) * 0.55733 / (r_DH * r_DH)
    print(f"\nForce on D from H: {F_DH:.2f} kJ/mol/nm")
    print(f"Direction: {'toward H (attractive)' if F_DH > 0 else 'away from H (repulsive)'}")
    
    # Spring force: F = -k*x
    F_spring = -k_kj_nm2 * dx
    print(f"Spring force on D: {F_spring:.2f} kJ/mol/nm")
    
    # Net force should be zero at equilibrium
    F_net = F_DH + F_spring
    print(f"Net force: {F_net:.2f} kJ/mol/nm (should be ~0)")
    
    # Test 2: With M-site contribution
    print("\n\nTest 2: Add M-site\n")
    
    # Add M-site
    m = pygcmc.MCAtom()
    m.x = 0.00793
    m.y = 0.00986
    m.z = 0.0
    m.charge = -1.11466
    m.type = 3
    state.atoms.append(m)
    state.activeAtomCount = 4
    
    # Reset Drude
    state.atoms[1].x = 0.0
    state.atoms[1].y = 0.0
    state.atoms[1].z = 0.0
    
    energy2, _ = pygcmc.computeSystemEnergyDrude(state)
    
    dx = state.atoms[1].x
    print(f"Drude position with M: x = {dx*1000:.6f} pm")
    print(f"Total energy: {energy2:.2f} kJ/mol")
    
    # Calculate O-M distance
    r_OM = np.sqrt(0.00793**2 + 0.00986**2)
    print(f"\nO-M distance: {r_OM*1000:.3f} pm ({r_OM*10:.4f} Å)")
    
    # Energy contributions
    E_OM = 138.935456 * 1.71636 * (-1.11466) / r_OM
    print(f"O-M Coulomb energy: {E_OM:.2f} kJ/mol")
    
    # This is the problem! O and M are very close and have opposite charges

if __name__ == "__main__":
    diagnose_forces()
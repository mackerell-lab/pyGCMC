#!/usr/bin/env python
"""Debug Drude SCF issues"""

import sys
import numpy as np
sys.path.append('/home/zhaomt/gcmc/test107/pygcmc_dev/build/modules/bindings')

import pygcmc

def test_drude_debug():
    print("=== Debugging Drude SCF ===\n")
    
    # Initialize
    pygcmc.initializeDrudeForce()
    
    # Create minimal system
    state = pygcmc.MCState()
    state.info.cutoff = 2.0
    state.info.box = [10.0, 10.0, 10.0]
    
    # Just two atoms: parent and Drude
    atoms = []
    
    # Parent atom
    parent = pygcmc.MCAtom()
    parent.x = 0.0
    parent.y = 0.0
    parent.z = 0.0
    parent.charge = 0.0  # Neutral parent
    parent.type = 0
    atoms.append(parent)
    
    # Drude particle
    drude = pygcmc.MCAtom()
    drude.x = 0.0
    drude.y = 0.0
    drude.z = 0.0
    drude.charge = -1.0  # Charged Drude
    drude.type = 1
    atoms.append(drude)
    
    # External test charge to pull Drude
    external = pygcmc.MCAtom()
    external.x = 0.5  # 5 Å away
    external.y = 0.0
    external.z = 0.0
    external.charge = 1.0  # Positive to attract negative Drude
    external.type = 2
    atoms.append(external)
    
    state.atoms = atoms
    state.activeAtomCount = 3
    
    # Create residues - parent+Drude in one, external in another
    res1 = pygcmc.MCResidue()
    res1.atomStart = 0
    res1.atomCount = 2  # Parent + Drude
    res1.active = True
    res1.fixed = False
    res1.type = 0
    
    res2 = pygcmc.MCResidue()
    res2.atomStart = 2
    res2.atomCount = 1  # External
    res2.active = True
    res2.fixed = False
    res2.type = 1
    
    state.residues = [res1, res2]
    state.activeResidueCount = 2
    
    # Force field
    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 3
    state.forcefield = ff
    
    # Add Drude with reasonable parameters
    # k = 1000 kcal/mol/Å² = 418400 kJ/mol/nm²
    k_kj_nm2 = 418400.0
    polarizability = 138.935456 * 1.0 * 1.0 / k_kj_nm2  # For charge = -1
    
    print(f"Setup:")
    print(f"  Parent at origin (charge = 0)")
    print(f"  Drude at origin (charge = -1)")
    print(f"  External at 5 Å (charge = +1)")
    print(f"  k = {k_kj_nm2} kJ/mol/nm²")
    print(f"  α = {polarizability:.6e} nm³\n")
    
    pygcmc.addDrudeParticle(
        drudeIndex=1,
        parentIndex=0,
        charge=-1.0,
        polarizability=polarizability
    )
    
    # Set SCF parameters
    params = pygcmc.DrudeSCFParams()
    params.tolerance = 1e-8
    params.maxIterations = 50
    params.dampingFactor = 0.5
    params.maxDrudeDistance = 0.1  # Large limit
    pygcmc.setDrudeSCFParameters(params)
    
    # Calculate energy
    print("Running SCF optimization...\n")
    energy, _ = pygcmc.computeSystemEnergyDrude(state)
    
    # Check Drude position
    dx = state.atoms[1].x - state.atoms[0].x
    print(f"Results:")
    print(f"  Drude displacement: {dx*1000:.6f} pm ({dx*10:.6f} Å)")
    print(f"  Energy: {energy:.2f} kJ/mol")
    
    # Calculate expected displacement
    # Force from external: F = k*q1*q2/r²
    F_ext = 138.935456 * (-1.0) * 1.0 / (0.5 * 0.5)
    # At equilibrium: F_ext = k_spring * x
    x_expected = F_ext / k_kj_nm2
    print(f"\nExpected displacement: {x_expected*1000:.6f} pm")
    print(f"Ratio actual/expected: {dx/x_expected if x_expected != 0 else 'N/A'}")
    
    # Test 2: Move Drude manually and check energy
    print("\n\nTest 2: Manual Drude displacement\n")
    
    # Reset and move Drude
    state.atoms[1].x = 0.01  # 0.1 Å
    
    # Calculate energy without SCF
    energy2, _ = pygcmc.computeSystemEnergyDrude(state)
    
    # Harmonic energy
    E_harm = 0.5 * k_kj_nm2 * (0.01 * 0.01)
    # Coulomb energy with external
    r_ext = 0.5 - 0.01  # Distance to external
    E_coulomb = 138.935456 * (-1.0) * 1.0 / r_ext
    
    print(f"With Drude at 0.1 Å:")
    print(f"  Total energy: {energy2:.2f} kJ/mol")
    print(f"  Expected harmonic: {E_harm:.2f} kJ/mol")
    print(f"  Expected Coulomb: {E_coulomb:.2f} kJ/mol")
    print(f"  Expected total: {E_harm + E_coulomb:.2f} kJ/mol")

if __name__ == "__main__":
    test_drude_debug()
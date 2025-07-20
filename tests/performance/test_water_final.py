#!/usr/bin/env python
"""Final complete energy calculation for SWM4-NDP water dimer"""

import sys
import numpy as np
sys.path.append('/home/zhaomt/gcmc/test107/pygcmc_dev/build/modules/bindings')

import pygcmc

def test_water_final():
    print("=== Final Water Dimer Energy Test ===\n")
    
    # Initialize Drude
    pygcmc.initializeDrudeForce()
    
    # Create state
    state = pygcmc.MCState()
    state.info.cutoff = 2.0
    state.info.box = [10.0, 10.0, 10.0]
    
    # Force field setup first
    state.forcefield.numTotalTypes = 4  # O, D, H, M
    state.forcefield.numMovementTypes = 4
    
    # SWM4-NDP parameters
    qO_core = 1.66260
    qD = -1.76260
    qM = -0.95710
    qH = 0.52855
    alpha_nm3 = 0.0013
    k_kj_nm2 = 166018.9
    
    # LJ parameters - only O has LJ
    sigma_O = 0.318395  # nm
    eps_O = 0.21094 * 4.184  # kJ/mol
    
    # Initialize LJ arrays (4x4 = 16)
    ljSigma = [0.0] * 16
    ljEps = [0.0] * 16
    
    # Set O-O interaction (type 0 with type 0)
    ljSigma[0] = sigma_O
    ljEps[0] = eps_O
    
    # Assign to forcefield
    state.forcefield.ljSigma = ljSigma
    state.forcefield.ljEps = ljEps
    
    print("Parameters:")
    print(f"  O: charge={qO_core:.4f}, σ={sigma_O:.6f} nm, ε={eps_O:.3f} kJ/mol")
    print(f"  D: charge={qD:.4f} (Drude particle)")
    print(f"  H: charge={qH:.4f}")
    print(f"  M: charge={qM:.4f} (virtual site)")
    print(f"  Polarizability: {alpha_nm3*1000:.1f} Å³")
    print()
    
    # Test at 3 Å
    separation = 0.3  # nm
    
    # Create atoms
    atoms = []
    
    # Water 1
    atoms.append(pygcmc.MCAtom())  # O
    atoms[-1].x = 0.0
    atoms[-1].y = 0.0
    atoms[-1].z = 0.0
    atoms[-1].charge = qO_core
    atoms[-1].type = 0
    
    atoms.append(pygcmc.MCAtom())  # D
    atoms[-1].x = 0.0
    atoms[-1].y = 0.0
    atoms[-1].z = 0.0
    atoms[-1].charge = qD
    atoms[-1].type = 1
    
    atoms.append(pygcmc.MCAtom())  # H1
    atoms[-1].x = 0.09572
    atoms[-1].y = 0.0
    atoms[-1].z = 0.0
    atoms[-1].charge = qH
    atoms[-1].type = 2
    
    atoms.append(pygcmc.MCAtom())  # H2
    atoms[-1].x = -0.023999
    atoms[-1].y = 0.092663
    atoms[-1].z = 0.0
    atoms[-1].charge = qH
    atoms[-1].type = 2
    
    atoms.append(pygcmc.MCAtom())  # M
    atoms[-1].x = 0.00793
    atoms[-1].y = 0.00986
    atoms[-1].z = 0.0
    atoms[-1].charge = qM
    atoms[-1].type = 3
    
    # Water 2
    atoms.append(pygcmc.MCAtom())  # O
    atoms[-1].x = separation
    atoms[-1].y = 0.0
    atoms[-1].z = 0.0
    atoms[-1].charge = qO_core
    atoms[-1].type = 0
    
    atoms.append(pygcmc.MCAtom())  # D
    atoms[-1].x = separation
    atoms[-1].y = 0.0
    atoms[-1].z = 0.0
    atoms[-1].charge = qD
    atoms[-1].type = 1
    
    atoms.append(pygcmc.MCAtom())  # H1
    atoms[-1].x = separation + 0.09572
    atoms[-1].y = 0.0
    atoms[-1].z = 0.0
    atoms[-1].charge = qH
    atoms[-1].type = 2
    
    atoms.append(pygcmc.MCAtom())  # H2
    atoms[-1].x = separation - 0.023999
    atoms[-1].y = 0.092663
    atoms[-1].z = 0.0
    atoms[-1].charge = qH
    atoms[-1].type = 2
    
    atoms.append(pygcmc.MCAtom())  # M
    atoms[-1].x = separation + 0.00793
    atoms[-1].y = 0.00986
    atoms[-1].z = 0.0
    atoms[-1].charge = qM
    atoms[-1].type = 3
    
    state.atoms = atoms
    state.activeAtomCount = 10
    
    # Residues
    res1 = pygcmc.MCResidue()
    res1.atomStart = 0
    res1.atomCount = 5
    res1.active = True
    res1.fixed = False
    res1.type = 0
    
    res2 = pygcmc.MCResidue()
    res2.atomStart = 5
    res2.atomCount = 5
    res2.active = True
    res2.fixed = False
    res2.type = 0
    
    state.residues = [res1, res2]
    state.activeResidueCount = 2
    
    # Add Drude particles
    pygcmc.addDrudeParticle(
        drudeIndex=1,
        parentIndex=0,
        charge=qD,
        polarizability=alpha_nm3
    )
    
    pygcmc.addDrudeParticle(
        drudeIndex=6,
        parentIndex=5,
        charge=qD,
        polarizability=alpha_nm3
    )
    
    # Set SCF parameters
    params = pygcmc.DrudeSCFParams()
    params.tolerance = 1e-8
    params.maxIterations = 100
    params.dampingFactor = 0.5
    params.maxDrudeDistance = 0.02
    pygcmc.setDrudeSCFParameters(params)
    
    print("Step 1: Energy without polarization")
    print("-" * 40)
    
    # Calculate without SCF
    pygcmc.computeSystemEnergyCutoff(state)
    
    # Each interaction is counted twice, so divide by 2
    vdw_no_pol = sum(res.energy_vdw for res in state.residues) / 2
    elec_no_pol = sum(res.energy_elec for res in state.residues) / 2
    
    print(f"VDW energy: {vdw_no_pol:.2f} kJ/mol")
    print(f"Coulomb energy: {elec_no_pol:.2f} kJ/mol")
    print(f"Total (no polarization): {vdw_no_pol + elec_no_pol:.2f} kJ/mol")
    
    print("\n\nStep 2: Energy with polarization (Drude SCF)")
    print("-" * 40)
    
    # Run SCF
    drude_energy, energy_dict = pygcmc.computeSystemEnergyDrude(state)
    
    # Recalculate with moved Drude particles
    pygcmc.computeSystemEnergyCutoff(state)
    
    vdw_with_pol = sum(res.energy_vdw for res in state.residues) / 2
    elec_with_pol = sum(res.energy_elec for res in state.residues) / 2
    
    print(f"VDW energy: {vdw_with_pol:.2f} kJ/mol")
    print(f"Coulomb energy: {elec_with_pol:.2f} kJ/mol")
    print(f"Drude harmonic: {drude_energy:.2f} kJ/mol")
    print(f"Total with polarization: {vdw_with_pol + elec_with_pol + drude_energy:.2f} kJ/mol")
    
    # Show Drude displacements
    print("\nDrude displacements:")
    for i, (drude_idx, parent_idx) in enumerate([(1, 0), (6, 5)]):
        dx = state.atoms[drude_idx].x - state.atoms[parent_idx].x
        dy = state.atoms[drude_idx].y - state.atoms[parent_idx].y
        dz = state.atoms[drude_idx].z - state.atoms[parent_idx].z
        dist = np.sqrt(dx*dx + dy*dy + dz*dz)
        print(f"  Water {i+1}: {dist*1000:.3f} pm")
    
    print("\n\nFinal Summary:")
    print("=" * 50)
    print(f"Water dimer at {separation*10:.1f} Å:")
    print(f"  VDW contribution: {vdw_with_pol:.2f} kJ/mol")
    print(f"  Coulomb (initial): {elec_no_pol:.2f} kJ/mol")
    print(f"  Coulomb (polarized): {elec_with_pol:.2f} kJ/mol")
    print(f"  Polarization change: {elec_with_pol - elec_no_pol:.2f} kJ/mol")
    print(f"  Drude restraint: {drude_energy:.2f} kJ/mol")
    print(f"\n  TOTAL INTERACTION ENERGY: {vdw_with_pol + elec_with_pol + drude_energy:.2f} kJ/mol")
    print("\nExpected range for H-bonded dimer: -20 to -30 kJ/mol")
    print("Note: Our configuration is not optimally H-bonded")

if __name__ == "__main__":
    test_water_final()
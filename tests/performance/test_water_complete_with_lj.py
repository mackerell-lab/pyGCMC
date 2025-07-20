#!/usr/bin/env python
"""Complete energy calculation including LJ for SWM4-NDP water"""

import sys
import numpy as np
sys.path.append('/home/zhaomt/gcmc/test107/pygcmc_dev/build/modules/bindings')

import pygcmc

def test_water_complete_with_lj():
    print("=== Complete Energy with LJ for SWM4-NDP Water ===\n")
    
    # Initialize Drude
    pygcmc.initializeDrudeForce()
    
    # Create state
    state = pygcmc.MCState()
    state.info.cutoff = 2.0
    state.info.box = [10.0, 10.0, 10.0]
    
    # SWM4-NDP parameters
    qO_core = 1.66260
    qD = -1.76260
    qM = -0.95710
    qH = 0.52855
    alpha_nm3 = 0.0013
    k_kj_nm2 = 166018.9
    
    # LJ parameters for SWM4-NDP (from OpenMM test)
    # Only O has LJ, others are zero
    sigma_O = 0.318395  # nm
    eps_O = 0.21094 * 4.184  # kJ/mol (converted from kcal/mol)
    
    print("Parameters:")
    print(f"  O: σ = {sigma_O} nm, ε = {eps_O:.3f} kJ/mol")
    print(f"  H, D, M: no LJ interactions")
    print()
    
    # Test at 3 Å separation
    separation = 0.3  # nm
    
    # Create two waters
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
    
    # Force field - 4 types: O, D, H, M
    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 4
    
    # Initialize LJ parameters (4x4 matrix)
    # Matrix is stored as [i*numTypes + j] for interaction between type i and j
    ff.ljSigma = [0.0] * 16
    ff.ljEps = [0.0] * 16
    
    # Set LJ parameters for O-O interaction (type 0 with type 0)
    # Need to set both (0,0) position
    ff.ljSigma[0 * 4 + 0] = sigma_O  # O-O interaction
    ff.ljEps[0 * 4 + 0] = eps_O      # O-O interaction
    
    # For debugging, let's also check if we need switching parameters
    ff.vdwSwitch = 0.9  # Default switching distance
    ff.vdwCutoff = 2.0  # Match state cutoff
    
    state.forcefield = ff
    
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
    
    print("Step 1: Calculate without Drude SCF")
    print("-" * 40)
    
    # First calculate energy with standard method (includes LJ)
    pygcmc.computeSystemEnergyCutoff(state)
    
    # Sum energies
    total_vdw = sum(res.energy_vdw for res in state.residues)
    total_elec = sum(res.energy_elec for res in state.residues)
    
    print(f"VDW energy: {total_vdw:.2f} kJ/mol")
    print(f"Coulomb energy: {total_elec:.2f} kJ/mol")
    print(f"Total (no polarization): {total_vdw + total_elec:.2f} kJ/mol")
    
    # Calculate expected LJ
    r_OO = separation
    sigma_r = sigma_O / r_OO
    lj_expected = 4 * eps_O * (sigma_r**12 - sigma_r**6)
    print(f"\nExpected O-O LJ: {lj_expected:.2f} kJ/mol")
    
    print("\n\nStep 2: Run Drude SCF")
    print("-" * 40)
    
    # Run SCF
    drude_energy, energy_dict = pygcmc.computeSystemEnergyDrude(state)
    
    # Recalculate with moved Drude particles
    pygcmc.computeSystemEnergyCutoff(state)
    
    total_vdw_scf = sum(res.energy_vdw for res in state.residues)
    total_elec_scf = sum(res.energy_elec for res in state.residues)
    
    print(f"VDW energy: {total_vdw_scf:.2f} kJ/mol (unchanged)")
    print(f"Coulomb energy: {total_elec_scf:.2f} kJ/mol")
    print(f"Drude harmonic: {drude_energy:.2f} kJ/mol")
    print(f"Total with polarization: {total_vdw_scf + total_elec_scf + drude_energy:.2f} kJ/mol")
    
    # Show Drude displacements
    print("\nDrude displacements:")
    for i, (drude_idx, parent_idx) in enumerate([(1, 0), (6, 5)]):
        dx = state.atoms[drude_idx].x - state.atoms[parent_idx].x
        dy = state.atoms[drude_idx].y - state.atoms[parent_idx].y
        dz = state.atoms[drude_idx].z - state.atoms[parent_idx].z
        dist = np.sqrt(dx*dx + dy*dy + dz*dz)
        print(f"  Water {i+1}: {dist*1000:.3f} pm")
    
    print("\n\nSummary:")
    print("-" * 40)
    print(f"Water dimer at {separation*10:.1f} Å:")
    print(f"  LJ contribution: {total_vdw_scf:.2f} kJ/mol")
    print(f"  Coulomb (no pol): {total_elec:.2f} kJ/mol")
    print(f"  Coulomb (with pol): {total_elec_scf:.2f} kJ/mol")
    print(f"  Polarization energy: {total_elec_scf - total_elec:.2f} kJ/mol")
    print(f"  Drude restraint: {drude_energy:.2f} kJ/mol")
    print(f"  Total interaction: {total_vdw_scf + total_elec_scf + drude_energy:.2f} kJ/mol")
    print("\nExpected range: -20 to -30 kJ/mol for hydrogen-bonded dimer")

if __name__ == "__main__":
    test_water_complete_with_lj()
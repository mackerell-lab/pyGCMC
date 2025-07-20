#!/usr/bin/env python
"""Comprehensive test of Drude SCF water model"""

import sys
import numpy as np
sys.path.append('/home/zhaomt/gcmc/test107/pygcmc_dev/build/modules/bindings')

import pygcmc

def test_drude_water_comprehensive():
    print("=== Comprehensive Drude SCF Water Model Test ===\n")
    
    # Initialize
    pygcmc.initializeDrudeForce()
    
    # Create state
    state = pygcmc.MCState()
    state.info.cutoff = 2.0  # 20 Å cutoff
    state.info.box = [10.0, 10.0, 10.0]  # 100 Å box
    
    # SWM4-NDP parameters from sold.psf
    qO_core = 1.66260    # 氧核心电荷
    qD = -1.76260        # Drude电荷  
    qM = -0.95710        # M-site电荷
    qH = 0.52855         # 氢电荷
    alpha_A3 = 1.30000   # Å³ from PSF
    thole = 1.3          # Thole参数
    
    # Convert units
    alpha_nm3 = alpha_A3 * 0.001
    CCELEC = 332.0716  # kcal*Å/mol/e²
    k_kcal = qD**2 * CCELEC / (2 * alpha_A3)
    k_kj_nm2 = k_kcal * 4.184 * 100.0
    
    print("SWM4-NDP Parameters:")
    print(f"  Charges: O={qO_core}, D={qD}, M={qM}, H={qH}")
    print(f"  Net O charge: {qO_core + qD:.3f}")
    print(f"  Polarizability: {alpha_A3} Å³ = {alpha_nm3} nm³")
    print(f"  Force constant: {k_kcal:.1f} kcal/mol/Å² = {k_kj_nm2:.1f} kJ/mol/nm²")
    print(f"  Thole parameter: {thole}\n")
    
    # Test 1: Single water molecule
    print("Test 1: Single Water Molecule")
    print("-" * 40)
    
    atoms = []
    
    # Water 1 at origin
    # Oxygen
    o1 = pygcmc.MCAtom()
    o1.x = 0.0
    o1.y = 0.0
    o1.z = 0.0
    o1.charge = qO_core
    o1.type = 0
    atoms.append(o1)
    
    # Drude on oxygen
    d1 = pygcmc.MCAtom()
    d1.x = 0.0
    d1.y = 0.0
    d1.z = 0.0
    d1.charge = qD
    d1.type = 1
    atoms.append(d1)
    
    # Hydrogen 1
    h11 = pygcmc.MCAtom()
    h11.x = 0.09572
    h11.y = 0.0
    h11.z = 0.0
    h11.charge = qH
    h11.type = 2
    atoms.append(h11)
    
    # Hydrogen 2
    h12 = pygcmc.MCAtom()
    h12.x = -0.023999
    h12.y = 0.092663
    h12.z = 0.0
    h12.charge = qH
    h12.type = 2
    atoms.append(h12)
    
    # M-site
    m1 = pygcmc.MCAtom()
    m1.x = 0.00793
    m1.y = 0.00986
    m1.z = 0.0
    m1.charge = qM
    m1.type = 3
    atoms.append(m1)
    
    # Create residue for water 1
    res1 = pygcmc.MCResidue()
    res1.atomStart = 0
    res1.atomCount = 5
    res1.active = True
    res1.fixed = False
    res1.type = 0
    
    state.atoms = atoms
    state.activeAtomCount = 5
    state.residues = [res1]
    state.activeResidueCount = 1
    
    # Force field
    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 4
    state.forcefield = ff
    
    # Add Drude particle
    pygcmc.addDrudeParticle(
        drudeIndex=1,
        parentIndex=0,
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
    
    # Calculate energy
    energy1, _ = pygcmc.computeSystemEnergyDrude(state)
    
    # Check Drude displacement
    dx = state.atoms[1].x - state.atoms[0].x
    dy = state.atoms[1].y - state.atoms[0].y
    dz = state.atoms[1].z - state.atoms[0].z
    dist1 = np.sqrt(dx*dx + dy*dy + dz*dz)
    
    # Calculate dipole moment
    # μ = Σ q_i * r_i
    dipole_x = 0.0
    dipole_y = 0.0
    dipole_z = 0.0
    for i in range(5):
        dipole_x += state.atoms[i].charge * state.atoms[i].x
        dipole_y += state.atoms[i].charge * state.atoms[i].y
        dipole_z += state.atoms[i].charge * state.atoms[i].z
    dipole_mag = np.sqrt(dipole_x**2 + dipole_y**2 + dipole_z**2)
    
    # Convert to Debye (1 e·nm = 4.80321 D)
    dipole_debye = dipole_mag * 4.80321 * 10  # Factor of 10 for nm to Å
    
    print(f"Energy: {energy1:.2f} kJ/mol")
    print(f"Drude displacement: {dist1*1000:.3f} pm")
    print(f"Dipole moment: {dipole_debye:.3f} D")
    print(f"Expected dipole (SWM4-NDP): ~2.4 D")
    
    # Test 2: Water dimer
    print("\n\nTest 2: Water Dimer")
    print("-" * 40)
    
    # Add second water at 3 Å
    offset = 0.3  # nm
    
    # Water 2
    o2 = pygcmc.MCAtom()
    o2.x = offset
    o2.y = 0.0
    o2.z = 0.0
    o2.charge = qO_core
    o2.type = 0
    atoms.append(o2)
    
    d2 = pygcmc.MCAtom()
    d2.x = offset
    d2.y = 0.0
    d2.z = 0.0
    d2.charge = qD
    d2.type = 1
    atoms.append(d2)
    
    h21 = pygcmc.MCAtom()
    h21.x = offset + 0.09572
    h21.y = 0.0
    h21.z = 0.0
    h21.charge = qH
    h21.type = 2
    atoms.append(h21)
    
    h22 = pygcmc.MCAtom()
    h22.x = offset - 0.023999
    h22.y = 0.092663
    h22.z = 0.0
    h22.charge = qH
    h22.type = 2
    atoms.append(h22)
    
    m2 = pygcmc.MCAtom()
    m2.x = offset + 0.00793
    m2.y = 0.00986
    m2.z = 0.0
    m2.charge = qM
    m2.type = 3
    atoms.append(m2)
    
    # Update state
    state.atoms = atoms
    state.activeAtomCount = 10
    
    # Add residue for water 2
    res2 = pygcmc.MCResidue()
    res2.atomStart = 5
    res2.atomCount = 5
    res2.active = True
    res2.fixed = False
    res2.type = 0
    
    state.residues = [res1, res2]
    state.activeResidueCount = 2
    
    # Add second Drude
    pygcmc.addDrudeParticle(
        drudeIndex=6,
        parentIndex=5,
        charge=qD,
        polarizability=alpha_nm3
    )
    
    # Add Thole screening between waters
    pygcmc.addDrudeScreenedPair(0, 1, thole)
    
    # Reset Drude positions
    state.atoms[1].x = 0.0
    state.atoms[6].x = offset
    
    # Calculate energy
    energy2, _ = pygcmc.computeSystemEnergyDrude(state)
    
    # Check Drude displacements
    dx1 = state.atoms[1].x - state.atoms[0].x
    dy1 = state.atoms[1].y - state.atoms[0].y
    dz1 = state.atoms[1].z - state.atoms[0].z
    dist_w1 = np.sqrt(dx1*dx1 + dy1*dy1 + dz1*dz1)
    
    dx2 = state.atoms[6].x - state.atoms[5].x
    dy2 = state.atoms[6].y - state.atoms[5].y
    dz2 = state.atoms[6].z - state.atoms[5].z
    dist_w2 = np.sqrt(dx2*dx2 + dy2*dy2 + dz2*dz2)
    
    # Interaction energy
    interaction_energy = energy2 - 2*energy1
    
    print(f"Total energy: {energy2:.2f} kJ/mol")
    print(f"Interaction energy: {interaction_energy:.2f} kJ/mol")
    print(f"Expected range: -20 to -30 kJ/mol for 3 Å separation")
    print(f"Water 1 Drude displacement: {dist_w1*1000:.3f} pm")
    print(f"Water 2 Drude displacement: {dist_w2*1000:.3f} pm")
    
    # Calculate individual dipole moments
    dipole1_x = sum(state.atoms[i].charge * state.atoms[i].x for i in range(5))
    dipole1_y = sum(state.atoms[i].charge * state.atoms[i].y for i in range(5))
    dipole1_z = sum(state.atoms[i].charge * state.atoms[i].z for i in range(5))
    dipole1_mag = np.sqrt(dipole1_x**2 + dipole1_y**2 + dipole1_z**2) * 48.0321
    
    dipole2_x = sum(state.atoms[i].charge * state.atoms[i].x for i in range(5, 10))
    dipole2_y = sum(state.atoms[i].charge * state.atoms[i].y for i in range(5, 10))
    dipole2_z = sum(state.atoms[i].charge * state.atoms[i].z for i in range(5, 10))
    dipole2_mag = np.sqrt(dipole2_x**2 + dipole2_y**2 + dipole2_z**2) * 48.0321
    
    print(f"Water 1 dipole: {dipole1_mag:.3f} D")
    print(f"Water 2 dipole: {dipole2_mag:.3f} D")
    
    # Test 3: External field response
    print("\n\nTest 3: Response to External Field")
    print("-" * 40)
    
    # Remove second water, add external charge
    state.atoms = atoms[:5]
    state.activeAtomCount = 5
    state.residues = [res1]
    state.activeResidueCount = 1
    
    # Clear and re-add single Drude
    pygcmc.clearDrudeForce()
    pygcmc.initializeDrudeForce()
    pygcmc.addDrudeParticle(
        drudeIndex=1,
        parentIndex=0,
        charge=qD,
        polarizability=alpha_nm3
    )
    pygcmc.setDrudeSCFParameters(params)
    
    # Add external charge
    ext = pygcmc.MCAtom()
    ext.x = 0.5  # 5 Å away
    ext.y = 0.0
    ext.z = 0.0
    ext.charge = 10.0  # Strong charge
    ext.type = 4
    state.atoms.append(ext)
    state.activeAtomCount = 6
    
    # External charge in separate residue
    res_ext = pygcmc.MCResidue()
    res_ext.atomStart = 5
    res_ext.atomCount = 1
    res_ext.active = True
    res_ext.fixed = False
    res_ext.type = 1
    state.residues = [res1, res_ext]
    state.activeResidueCount = 2
    
    # Reset Drude
    state.atoms[1].x = 0.0
    state.atoms[1].y = 0.0
    state.atoms[1].z = 0.0
    
    # Calculate
    energy_field, _ = pygcmc.computeSystemEnergyDrude(state)
    
    # Check Drude displacement
    dx_field = state.atoms[1].x - state.atoms[0].x
    dist_field = abs(dx_field)
    
    # Calculate induced dipole
    dipole_field_x = sum(state.atoms[i].charge * state.atoms[i].x for i in range(5))
    dipole_field_mag = abs(dipole_field_x) * 48.0321
    
    print(f"External charge: +10e at 5 Å")
    print(f"Drude displacement: {dist_field*1000:.3f} pm")
    print(f"Direction: {'toward' if dx_field < 0 else 'away from'} external charge")
    print(f"Induced dipole: {dipole_field_mag:.3f} D")
    print(f"Energy: {energy_field:.2f} kJ/mol")
    
    # Summary
    print("\n\nSummary of Results")
    print("-" * 40)
    print(f"Single water:")
    print(f"  Energy: {energy1:.2f} kJ/mol (should be ~0 if intramolecular excluded)")
    print(f"  Dipole: {dipole_debye:.3f} D (expected ~2.4 D)")
    print(f"Water dimer:")
    print(f"  Interaction: {interaction_energy:.2f} kJ/mol (expected -20 to -30)")
    print(f"  Mutual polarization observed: {'Yes' if dist_w1 > 1e-6 else 'No'}")
    print(f"External field:")
    print(f"  Response: {'Correct' if dx_field < 0 else 'Incorrect'} (D should move toward +charge)")
    print(f"  Polarizability working: {'Yes' if dist_field > 1e-6 else 'No'}")

if __name__ == "__main__":
    test_drude_water_comprehensive()
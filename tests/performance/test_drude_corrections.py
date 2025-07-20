#!/usr/bin/env python
"""Test Drude corrections based on CHARMM/OpenMM analysis"""

import sys
import numpy as np
sys.path.append('/home/zhaomt/gcmc/test107/pygcmc_dev/build/modules/bindings')

import pygcmc

def test_drude_corrections():
    print("=== Testing Drude Corrections Based on CHARMM/OpenMM ===\n")
    
    # Initialize
    pygcmc.initializeDrudeForce()
    
    # Create state
    state = pygcmc.MCState()
    state.info.cutoff = 2.0
    state.info.box = [10.0, 10.0, 10.0]
    
    # Test 1: Check polarizability and charge relationship
    print("Test 1: Polarizability-Charge-Force Constant Relationship\n")
    
    # CHARMM formula: q_D = sqrt(2*alpha*k_drude/CCELEC)
    # Where CCELEC = 332.0716 (CHARMM electrostatic constant in kcal*Å/mol/e²)
    # For SWM4-NDP: alpha = 0.97825258 Å³, k = 1000 kcal/mol/Å²
    
    alpha_A3 = 0.97825258  # Å³
    k_kcal = 1000.0  # kcal/mol/Å²
    CCELEC = 332.0716  # kcal*Å/mol/e²
    
    # Calculate expected Drude charge
    q_D_expected = -np.sqrt(2 * alpha_A3 * k_kcal / CCELEC)
    print(f"From CHARMM formula:")
    print(f"  α = {alpha_A3} Å³")
    print(f"  k = {k_kcal} kcal/mol/Å²")
    print(f"  Expected q_D = {q_D_expected:.6f}")
    print(f"  SWM4-NDP q_D = -1.71636")
    print(f"  Ratio: {-1.71636/q_D_expected:.3f}")
    
    # Test 2: Verify Thole screening
    print("\n\nTest 2: Thole Screening Function\n")
    
    # Create two atoms with Drude particles
    atoms = []
    
    # First atom pair
    parent1 = pygcmc.MCAtom()
    parent1.x = 0.0
    parent1.y = 0.0
    parent1.z = 0.0
    parent1.charge = 1.0
    parent1.type = 0
    atoms.append(parent1)
    
    drude1 = pygcmc.MCAtom()
    drude1.x = 0.0
    drude1.y = 0.0
    drude1.z = 0.0
    drude1.charge = -1.0
    drude1.type = 1
    atoms.append(drude1)
    
    # Second atom pair at 3 Å
    parent2 = pygcmc.MCAtom()
    parent2.x = 0.3  # 3 Å
    parent2.y = 0.0
    parent2.z = 0.0
    parent2.charge = 1.0
    parent2.type = 0
    atoms.append(parent2)
    
    drude2 = pygcmc.MCAtom()
    drude2.x = 0.3
    drude2.y = 0.0
    drude2.z = 0.0
    drude2.charge = -1.0
    drude2.type = 1
    atoms.append(drude2)
    
    state.atoms = atoms
    state.activeAtomCount = 4
    
    # Create residues (each dipole is separate molecule)
    res1 = pygcmc.MCResidue()
    res1.atomStart = 0
    res1.atomCount = 2
    res1.active = True
    res1.fixed = False
    res1.type = 0
    
    res2 = pygcmc.MCResidue()
    res2.atomStart = 2
    res2.atomCount = 2
    res2.active = True
    res2.fixed = False
    res2.type = 0
    
    state.residues = [res1, res2]
    state.activeResidueCount = 2
    
    # Force field
    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 2
    state.forcefield = ff
    
    # Test with isotropic polarizability
    polarizability = 0.001  # nm³
    
    idx1 = pygcmc.addDrudeParticle(
        drudeIndex=1,
        parentIndex=0,
        charge=-1.0,
        polarizability=polarizability
    )
    
    idx2 = pygcmc.addDrudeParticle(
        drudeIndex=3,
        parentIndex=2,
        charge=-1.0,
        polarizability=polarizability
    )
    
    # Add screened pair with standard Thole parameter
    thole = 1.3  # Standard value for SWM4-NDP
    pygcmc.addDrudeScreenedPair(idx1, idx2, thole)
    
    # Set SCF parameters
    params = pygcmc.DrudeSCFParams()
    params.tolerance = 1e-8
    params.maxIterations = 100
    params.dampingFactor = 0.5
    params.maxDrudeDistance = 0.02
    pygcmc.setDrudeSCFParameters(params)
    
    # Calculate energy
    energy, _ = pygcmc.computeSystemEnergyDrude(state)
    
    # Check Drude displacements
    dx1 = state.atoms[1].x - state.atoms[0].x
    dx2 = state.atoms[3].x - state.atoms[2].x
    
    print(f"At 3 Å separation with Thole = {thole}:")
    print(f"  Drude 1 displacement: {dx1*1000:.3f} pm")
    print(f"  Drude 2 displacement: {dx2*1000:.3f} pm")
    print(f"  Energy: {energy:.2f} kJ/mol")
    
    # Calculate Thole screening parameter
    r = 0.3  # nm
    u = r * thole / (polarizability ** (1.0/3.0))
    screening = 1.0 - (1.0 + 0.5*u) * np.exp(-u)
    print(f"\nThole screening calculation:")
    print(f"  u = r × thole / α^(1/3) = {u:.3f}")
    print(f"  Screening factor = {screening:.3f}")
    print(f"  Effective interaction = {screening*100:.1f}% of bare Coulomb")
    
    # Test 3: Compare with no screening
    print("\n\nTest 3: Effect of Thole Screening\n")
    
    # Clear and recreate without screening
    pygcmc.clearDrudeForce()
    pygcmc.initializeDrudeForce()
    
    # Add same Drude particles but no screening
    pygcmc.addDrudeParticle(
        drudeIndex=1,
        parentIndex=0,
        charge=-1.0,
        polarizability=polarizability
    )
    
    pygcmc.addDrudeParticle(
        drudeIndex=3,
        parentIndex=2,
        charge=-1.0,
        polarizability=polarizability
    )
    
    # Reset positions
    state.atoms[1].x = 0.0
    state.atoms[3].x = 0.3
    
    pygcmc.setDrudeSCFParameters(params)
    energy_no_screen, _ = pygcmc.computeSystemEnergyDrude(state)
    
    dx1_no = state.atoms[1].x - state.atoms[0].x
    dx2_no = state.atoms[3].x - state.atoms[2].x
    
    print(f"Without Thole screening:")
    print(f"  Drude 1 displacement: {dx1_no*1000:.3f} pm")
    print(f"  Drude 2 displacement: {dx2_no*1000:.3f} pm")
    print(f"  Energy: {energy_no_screen:.2f} kJ/mol")
    print(f"\nScreening reduces displacement by: {(1 - abs(dx1/dx1_no))*100:.1f}%")

if __name__ == "__main__":
    test_drude_corrections()
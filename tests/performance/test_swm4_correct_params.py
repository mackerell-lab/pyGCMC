#!/usr/bin/env python
"""Test SWM4-NDP with correct parameters from PSF file"""

import sys
import numpy as np
sys.path.append('/home/zhaomt/gcmc/test107/pygcmc_dev/build/modules/bindings')

import pygcmc

def test_swm4_correct_params():
    print("=== Testing SWM4-NDP with Correct Parameters ===\n")
    
    # Initialize
    pygcmc.initializeDrudeForce()
    
    # Create state
    state = pygcmc.MCState()
    state.info.cutoff = 2.0
    state.info.box = [10.0, 10.0, 10.0]
    
    # SWM4-NDP parameters from sold.psf
    qO_core = 1.66260    # 氧核心电荷
    qD = -1.76260        # Drude电荷  
    qM = -0.95710        # M-site电荷
    qH = 0.52855         # 氢电荷
    alpha_A3 = 1.30000   # Å³ from PSF
    thole = 1.3          # Thole参数
    
    # Convert alpha to nm³
    alpha_nm3 = alpha_A3 * 0.001
    
    # Calculate force constant from charge and polarizability
    # q_D = -sqrt(2*alpha*k/CCELEC)
    # k = q_D² * CCELEC / (2*alpha)
    CCELEC = 332.0716  # kcal*Å/mol/e²
    k_kcal = qD**2 * CCELEC / (2 * alpha_A3)
    k_kj_nm2 = k_kcal * 4.184 * 100.0  # Convert to kJ/mol/nm²
    
    print(f"Calculated parameters:")
    print(f"  α = {alpha_A3} Å³ = {alpha_nm3} nm³")
    print(f"  q_D = {qD}")
    print(f"  k = {k_kcal:.1f} kcal/mol/Å² = {k_kj_nm2:.1f} kJ/mol/nm²")
    print(f"  (Note: NOT 1000 kcal/mol/Å²!)\n")
    
    # Verify with CHARMM formula
    q_D_check = -np.sqrt(2 * alpha_A3 * k_kcal / CCELEC)
    print(f"Verification: q_D from CHARMM formula = {q_D_check:.5f}")
    print(f"Matches PSF value: {abs(q_D_check - qD) < 0.0001}\n")
    
    # Create single water molecule
    atoms = []
    
    # Oxygen
    o = pygcmc.MCAtom()
    o.x = 0.0
    o.y = 0.0
    o.z = 0.0
    o.charge = qO_core
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
    
    # Create residue
    res = pygcmc.MCResidue()
    res.atomStart = 0
    res.atomCount = 5
    res.active = True
    res.fixed = False
    res.type = 0
    
    state.atoms = atoms
    state.activeAtomCount = 5
    state.residues = [res]
    state.activeResidueCount = 1
    
    # Force field
    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 4
    state.forcefield = ff
    
    # Add Drude particle with CORRECT polarizability
    # Note: We pass the polarizability in nm³
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
    
    print("Test 1: Single water molecule\n")
    energy, _ = pygcmc.computeSystemEnergyDrude(state)
    
    # Check Drude displacement
    dx = state.atoms[1].x - state.atoms[0].x
    dy = state.atoms[1].y - state.atoms[0].y
    dz = state.atoms[1].z - state.atoms[0].z
    dist = np.sqrt(dx*dx + dy*dy + dz*dz)
    
    print(f"Energy: {energy:.2f} kJ/mol")
    print(f"Drude displacement: {dist*1000:.3f} pm")
    
    # Test 2: Add external charge
    print("\n\nTest 2: Water in external field\n")
    
    # Add test charge
    ext = pygcmc.MCAtom()
    ext.x = 0.5  # 5 Å away
    ext.y = 0.0
    ext.z = 0.0
    ext.charge = 1.0
    ext.type = 4
    state.atoms.append(ext)
    state.activeAtomCount = 6
    
    # Create separate residue
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
    
    print(f"Energy with external charge: {energy2:.2f} kJ/mol")
    print(f"Drude displacement: {dist2*1000:.3f} pm ({dist2*10:.4f} Å)")
    
    # Calculate expected displacement
    # Force from external on Drude: F = k*qD*qext/r²
    r_ext = 0.5  # nm
    F_ext = 138.935456 * qD * 1.0 / (r_ext * r_ext)
    x_expected = F_ext / k_kj_nm2
    
    print(f"\nExpected displacement: {x_expected*1000:.3f} pm")
    print(f"Ratio actual/expected: {dist2/x_expected:.3f}")

if __name__ == "__main__":
    test_swm4_correct_params()
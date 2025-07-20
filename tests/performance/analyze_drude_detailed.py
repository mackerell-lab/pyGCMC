#!/usr/bin/env python
"""Detailed analysis of Drude energy components"""

import sys
import numpy as np
import math
sys.path.append('/home/zhaomt/gcmc/test107/pygcmc_dev/build/modules/bindings')

import pygcmc

def analyze_water_dimer_step_by_step():
    """Analyze water dimer energy step by step"""
    print("=== Detailed Water Dimer Analysis ===\n")
    
    # Create simple dimer
    state = pygcmc.MCState()
    state.info.box = [10.0, 10.0, 10.0]
    state.info.cutoff = 5.0
    
    # Initialize Drude
    pygcmc.initializeDrudeForce()
    
    # Use very relaxed parameters to ensure convergence
    scf_params = pygcmc.DrudeSCFParams()
    scf_params.tolerance = 10.0  
    scf_params.maxIterations = 200
    scf_params.dampingFactor = 0.3
    scf_params.maxDrudeDistance = 0.02
    pygcmc.setDrudeSCFParameters(scf_params)
    
    atoms = []
    residues = []
    
    # SWM4-NDP parameters
    qO = 1.71636
    qD = -1.71636
    qH = 0.55733
    qM = -1.11466
    rOH = 0.09572  # nm
    aHOH = 104.52 * math.pi / 180
    
    # Create two waters at 3.0 Å separation
    water_positions = [(5.0, 5.0, 5.0), (5.3, 5.0, 5.0)]
    
    for water_id, (x, y, z) in enumerate(water_positions):
        # Create water atoms
        # Oxygen
        o = pygcmc.MCAtom()
        o.x, o.y, o.z = x, y, z
        o.charge = qO
        o.type = 0
        atoms.append(o)
        
        # Drude - initially at parent
        d = pygcmc.MCAtom()
        d.x, d.y, d.z = x, y, z
        d.charge = qD
        d.type = 1
        atoms.append(d)
        
        # H1 - along x
        h1 = pygcmc.MCAtom()
        h1.x = x + rOH
        h1.y = y
        h1.z = z
        h1.charge = qH
        h1.type = 2
        atoms.append(h1)
        
        # H2 - in xy plane
        h2 = pygcmc.MCAtom()
        h2.x = x + rOH * math.cos(aHOH)
        h2.y = y + rOH * math.sin(aHOH)
        h2.z = z
        h2.charge = qH
        h2.type = 2
        atoms.append(h2)
        
        # M-site
        m = pygcmc.MCAtom()
        w_O = 0.786646558
        w_H = 0.106676721
        m.x = w_O * o.x + w_H * h1.x + w_H * h2.x
        m.y = w_O * o.y + w_H * h1.y + w_H * h2.y
        m.z = w_O * o.z + w_H * h1.z + w_H * h2.z
        m.charge = qM
        m.type = 3
        atoms.append(m)
        
        # Create residue
        res = pygcmc.MCResidue()
        res.atomStart = 5 * water_id
        res.atomCount = 5
        res.active = True
        res.type = 0
        residues.append(res)
        
        # Add Drude
        pygcmc.addDrudeParticle(
            drudeIndex=5*water_id + 1,
            parentIndex=5*water_id,
            charge=qD,
            polarizability=0.000978253
        )
    
    state.atoms = atoms
    state.activeAtomCount = len(atoms)
    state.residues = residues
    state.activeResidueCount = len(residues)
    
    # Add Thole screening
    pygcmc.addDrudeScreenedPair(0, 1, 1.3)
    
    # Step 1: Check initial configuration
    print("Step 1: Initial Configuration")
    print("-" * 40)
    print(f"O-O distance: {water_positions[1][0] - water_positions[0][0]:.3f} nm")
    
    # Step 2: Calculate energy without SCF (Drude at parent)
    print("\nStep 2: Energy with Drude at parent positions")
    print("-" * 40)
    
    result = pygcmc.computeSystemEnergyDrude(state)
    if isinstance(result, tuple):
        energy, components = result
    else:
        energy = result
        components = {}
    
    print(f"Total Drude energy: {energy:.4f} kJ/mol")
    print(f"Components: {components}")
    
    # Check Drude positions after SCF
    print("\nDrude positions after SCF:")
    for i in range(2):
        o_idx = 5 * i
        d_idx = 5 * i + 1
        dx = state.atoms[d_idx].x - state.atoms[o_idx].x
        dy = state.atoms[d_idx].y - state.atoms[o_idx].y
        dz = state.atoms[d_idx].z - state.atoms[o_idx].z
        dist = math.sqrt(dx*dx + dy*dy + dz*dz)
        print(f"  Water {i+1}: Drude displacement = {dist*1000:.3f} pm")
        print(f"    Direction: ({dx/dist:.3f}, {dy/dist:.3f}, {dz/dist:.3f})")
    
    # Step 3: Calculate interaction energy
    print("\nStep 3: Interaction Energy")
    print("-" * 40)
    
    # Energy of dimer
    E_dimer = energy
    
    # Energy of isolated waters
    # Water 1 alone
    state1 = pygcmc.MCState()
    state1.info.box = state.info.box
    state1.info.cutoff = state.info.cutoff
    state1.atoms = atoms[:5]
    state1.activeAtomCount = 5
    state1.residues = [residues[0]]
    state1.activeResidueCount = 1
    
    pygcmc.initializeDrudeForce()
    pygcmc.setDrudeSCFParameters(scf_params)
    pygcmc.addDrudeParticle(
        drudeIndex=1,
        parentIndex=0,
        charge=qD,
        polarizability=0.000978253
    )
    
    result1 = pygcmc.computeSystemEnergyDrude(state1)
    E_water1 = result1[0] if isinstance(result1, tuple) else result1
    
    # Water 2 alone
    state2 = pygcmc.MCState()
    state2.info.box = state.info.box
    state2.info.cutoff = state.info.cutoff
    state2.atoms = atoms[5:]
    state2.activeAtomCount = 5
    res2 = pygcmc.MCResidue()
    res2.atomStart = 0
    res2.atomCount = 5
    res2.active = True
    res2.type = 0
    state2.residues = [res2]
    state2.activeResidueCount = 1
    
    pygcmc.initializeDrudeForce()
    pygcmc.setDrudeSCFParameters(scf_params)
    pygcmc.addDrudeParticle(
        drudeIndex=1,
        parentIndex=0,
        charge=qD,
        polarizability=0.000978253
    )
    
    result2 = pygcmc.computeSystemEnergyDrude(state2)
    E_water2 = result2[0] if isinstance(result2, tuple) else result2
    
    print(f"E(dimer): {E_dimer:.4f} kJ/mol")
    print(f"E(water1): {E_water1:.4f} kJ/mol")
    print(f"E(water2): {E_water2:.4f} kJ/mol")
    print(f"Interaction energy: {E_dimer - E_water1 - E_water2:.4f} kJ/mol")
    
    # Step 4: Analyze force constants
    print("\nStep 4: Force Constant Analysis")
    print("-" * 40)
    
    # Calculate expected force constant
    ONE_4PI_EPS0 = 138.935456  # kJ/mol·nm·e^-2
    alpha = 0.000978253  # nm^3
    k_calc = ONE_4PI_EPS0 * qD * qD / alpha
    
    print(f"Drude charge: {qD} e")
    print(f"Polarizability: {alpha} nm³")
    print(f"Calculated k: {k_calc:.1f} kJ/mol/nm²")
    print(f"In kcal/mol/Å²: {k_calc / 418.4:.1f}")
    
    # Expected from CHARMM
    print(f"\nExpected from CHARMM:")
    print(f"  k = 1000 kcal/mol/Å² = 418400 kJ/mol/nm²")
    print(f"  Our value matches: {abs(k_calc - 418400) < 1000}")

def test_simple_systems():
    """Test very simple systems to verify implementation"""
    print("\n\n=== Simple System Tests ===\n")
    
    # Test 1: Two opposite charges
    print("Test 1: Two opposite charges (no Drude)")
    state = pygcmc.MCState()
    state.info.box = [10.0, 10.0, 10.0]
    state.info.cutoff = 5.0
    
    atoms = []
    # +1 charge at origin
    a1 = pygcmc.MCAtom()
    a1.x, a1.y, a1.z = 5.0, 5.0, 5.0
    a1.charge = 1.0
    a1.type = 0
    atoms.append(a1)
    
    # -1 charge at 0.5 nm
    a2 = pygcmc.MCAtom()
    a2.x, a2.y, a2.z = 5.5, 5.0, 5.0
    a2.charge = -1.0
    a2.type = 0
    atoms.append(a2)
    
    state.atoms = atoms
    state.activeAtomCount = 2
    
    # Create residues (each in separate residue)
    res1 = pygcmc.MCResidue()
    res1.atomStart = 0
    res1.atomCount = 1
    res1.active = True
    res1.type = 0
    
    res2 = pygcmc.MCResidue()
    res2.atomStart = 1
    res2.atomCount = 1
    res2.active = True
    res2.type = 0
    
    state.residues = [res1, res2]
    state.activeResidueCount = 2
    
    # Set empty force field
    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.ljSigma = [0.0]
    ff.ljEps = [0.0]
    state.forcefield = ff
    
    # Calculate Coulomb energy
    pygcmc.computeSystemEnergyCutoff(state)
    coulomb = sum(res.energy_elec for res in state.residues)
    
    # Expected: E = k*q1*q2/r = 138.935456 * 1 * (-1) / 0.5 = -277.87 kJ/mol
    expected = -138.935456 / 0.5
    
    print(f"  Calculated: {coulomb:.2f} kJ/mol")
    print(f"  Expected: {expected:.2f} kJ/mol")
    print(f"  Match: {abs(coulomb - expected) < 0.1}")

if __name__ == "__main__":
    analyze_water_dimer_step_by_step()
    test_simple_systems()
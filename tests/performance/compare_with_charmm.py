#!/usr/bin/env python
"""Compare PyGCMC Drude results with CHARMM reference values"""

import sys
import numpy as np
import math
sys.path.append('/home/zhaomt/gcmc/test107/pygcmc_dev/build/modules/bindings')

import pygcmc

def create_water_dimer():
    """Create water dimer for comparison with CHARMM"""
    state = pygcmc.MCState()
    state.info.box = [10.0, 10.0, 10.0]
    state.info.cutoff = 5.0
    
    # Initialize Drude force with relaxed parameters for convergence
    pygcmc.initializeDrudeForce()
    
    # Use slightly relaxed parameters for better convergence
    scf_params = pygcmc.DrudeSCFParams()
    scf_params.tolerance = 5.0  # Slightly relaxed from OpenMM default
    scf_params.maxIterations = 100  # More iterations
    scf_params.dampingFactor = 0.5
    scf_params.maxDrudeDistance = 0.02
    pygcmc.setDrudeSCFParameters(scf_params)
    
    atoms = []
    residues = []
    
    # SWM4-NDP parameters (from CHARMM)
    qO = 1.71636
    qD = -1.71636
    qH = 0.55733
    qM = -1.11466
    rOH = 0.09572  # nm
    aHOH = 104.52 * math.pi / 180
    
    # Water 1 at origin
    water1_pos = [5.0, 5.0, 5.0]
    
    # Water 2 - positioned for optimal H-bond
    # O-O distance ~0.28 nm (2.8 Å)
    water2_pos = [5.28, 5.0, 5.0]
    
    # Create two waters
    for water_id, (x, y, z) in enumerate([water1_pos, water2_pos]):
        # Oxygen
        o = pygcmc.MCAtom()
        o.x, o.y, o.z = x, y, z
        o.charge = qO
        o.type = 0
        atoms.append(o)
        
        # Drude
        d = pygcmc.MCAtom()
        d.x, d.y, d.z = x, y, z
        d.charge = qD
        d.type = 1
        atoms.append(d)
        
        # For water 2, rotate to form H-bond
        if water_id == 1:
            # Point H toward water 1's oxygen
            h1_x = x - rOH  # H pointing back
            h1_y = y
            h1_z = z
            
            h2_x = x - rOH * math.cos(aHOH)
            h2_y = y + rOH * math.sin(aHOH)
            h2_z = z
        else:
            # Normal orientation
            h1_x = x + rOH
            h1_y = y
            h1_z = z
            
            h2_x = x + rOH * math.cos(aHOH)
            h2_y = y + rOH * math.sin(aHOH)
            h2_z = z
        
        # Hydrogen 1
        h1 = pygcmc.MCAtom()
        h1.x, h1.y, h1.z = h1_x, h1_y, h1_z
        h1.charge = qH
        h1.type = 2
        atoms.append(h1)
        
        # Hydrogen 2
        h2 = pygcmc.MCAtom()
        h2.x, h2.y, h2.z = h2_x, h2_y, h2_z
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
        
        # Add Drude particle
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
    
    # Set force field parameters for LJ
    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 4  # O, D, H, M
    ff.numMovementTypes = 4
    
    # LJ parameters (simplified - only O-O interaction)
    # From OpenMM: sigma=0.318395 nm, epsilon=0.88257 kJ/mol
    ljSigma = [0.0] * 16
    ljEps = [0.0] * 16
    
    # O-O interaction
    ljSigma[0] = 0.318395
    ljEps[0] = 0.88257
    
    ff.ljSigma = ljSigma
    ff.ljEps = ljEps
    state.forcefield = ff
    
    # Add Thole screening
    pygcmc.addDrudeScreenedPair(0, 1, 1.3)
    
    return state

def main():
    print("=== PyGCMC vs CHARMM Comparison ===\n")
    
    # Test 1: Water dimer
    print("1. Water Dimer Test")
    print("-" * 40)
    
    state = create_water_dimer()
    
    # Calculate energy
    result = pygcmc.computeSystemEnergyDrude(state)
    
    if isinstance(result, tuple):
        energy = result[0]
        components = result[1]
    else:
        energy = result
        components = {}
    
    print(f"\nPyGCMC Results:")
    print(f"  Total Drude energy: {energy:.4f} kJ/mol")
    print(f"  Components: {components}")
    
    # Also calculate Coulomb and LJ separately
    pygcmc.computeSystemEnergyCutoff(state)
    coulomb_energy = sum(res.energy_elec for res in state.residues)
    lj_energy = sum(res.energy_vdw for res in state.residues)
    
    print(f"\nAdditional energy components:")
    print(f"  Coulomb energy: {coulomb_energy:.4f} kJ/mol")
    print(f"  LJ energy: {lj_energy:.4f} kJ/mol")
    print(f"  Total non-Drude: {coulomb_energy + lj_energy:.4f} kJ/mol")
    
    print(f"\nCHARMM Reference Values (SWM4-NDP):")
    print(f"  Water dimer at 2.8 Å:")
    print(f"    Interaction energy: ~-21 kJ/mol")
    print(f"    Includes polarization effects")
    
    # Test 2: Single water properties
    print("\n\n2. Single Water Molecule Test")
    print("-" * 40)
    
    # Create single water
    state_single = pygcmc.MCState()
    state_single.info.box = [10.0, 10.0, 10.0]
    state_single.info.cutoff = 5.0
    
    pygcmc.initializeDrudeForce()
    scf_params = pygcmc.DrudeSCFParams()
    scf_params.tolerance = 5.0
    scf_params.maxIterations = 100
    pygcmc.setDrudeSCFParameters(scf_params)
    
    # Just take first water from dimer
    state_single.atoms = state.atoms[:5]
    state_single.activeAtomCount = 5
    state_single.residues = [state.residues[0]]
    state_single.activeResidueCount = 1
    
    # Re-add Drude particle
    pygcmc.initializeDrudeForce()
    pygcmc.setDrudeSCFParameters(scf_params)
    pygcmc.addDrudeParticle(
        drudeIndex=1,
        parentIndex=0,
        charge=-1.71636,
        polarizability=0.000978253
    )
    
    result_single = pygcmc.computeSystemEnergyDrude(state_single)
    
    if isinstance(result_single, tuple):
        energy_single = result_single[0]
    else:
        energy_single = result_single
    
    print(f"\nPyGCMC Results:")
    print(f"  Single water energy: {energy_single:.4f} kJ/mol")
    print(f"  (Should be 0.0 for isolated molecule)")
    
    # Test 3: Performance metrics
    print("\n\n3. Performance Metrics")
    print("-" * 40)
    
    import time
    
    # Time 100 energy evaluations
    n_evals = 100
    start = time.time()
    for _ in range(n_evals):
        pygcmc.computeSystemEnergyDrude(state)
    elapsed = time.time() - start
    
    print(f"\nWater dimer energy calculation:")
    print(f"  {n_evals} evaluations in {elapsed:.3f} s")
    print(f"  Average: {elapsed/n_evals*1000:.2f} ms per evaluation")
    print(f"  Rate: {n_evals/elapsed:.1f} evaluations/s")
    
    # Summary
    print("\n\n=== SUMMARY ===")
    print("\nKey observations:")
    print("1. Single water should have 0 energy (only intramolecular)")
    print("2. Water dimer should be around -21 kJ/mol")
    print("3. Performance depends on SCF convergence")
    
    if abs(energy_single) < 0.1:
        print("\n✓ Single water energy correct (near 0)")
    else:
        print("\n✗ Single water energy incorrect")
    
    if -30 < energy < -15:
        print("✓ Dimer energy in reasonable range")
    else:
        print("✗ Dimer energy seems off")

if __name__ == "__main__":
    main()
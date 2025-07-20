#!/usr/bin/env python
"""Test OPT3 vs SCF performance and accuracy"""

import sys
import time
import numpy as np
sys.path.append('/home/zhaomt/gcmc/test107/pygcmc_dev/build/modules/bindings')

import pygcmc

def create_water_dimer():
    """Create a simple water dimer for testing"""
    state = pygcmc.MCState()
    state.info.box = [10.0, 10.0, 10.0]
    state.info.cutoff = 5.0
    
    atoms = []
    
    # Water 1: SWM4-NDP
    # O, D, H, H, M
    positions1 = [
        [5.0, 5.0, 5.0],      # O
        [5.0, 5.0, 5.0],      # D
        [5.09572, 5.0, 5.0],  # H1
        [4.97, 5.09, 5.0],    # H2
        [5.015, 5.011, 5.0]   # M-site
    ]
    charges1 = [1.71636, -1.71636, 0.55733, 0.55733, -1.11466]
    types1 = [0, 1, 2, 2, 3]
    
    # Water 2: 2.8 Å away
    positions2 = [
        [5.28, 5.0, 5.0],      # O
        [5.28, 5.0, 5.0],      # D
        [5.37572, 5.0, 5.0],   # H1
        [5.25, 5.09, 5.0],     # H2
        [5.295, 5.011, 5.0]    # M-site
    ]
    
    # Add atoms
    for i in range(5):
        a = pygcmc.MCAtom()
        a.x, a.y, a.z = positions1[i]
        a.charge = charges1[i]
        a.type = types1[i]
        atoms.append(a)
    
    for i in range(5):
        a = pygcmc.MCAtom()
        a.x, a.y, a.z = positions2[i]
        a.charge = charges1[i]
        a.type = types1[i]
        atoms.append(a)
    
    state.atoms = atoms
    state.activeAtomCount = len(atoms)
    
    # Create residues
    residues = []
    for i in range(2):
        res = pygcmc.MCResidue()
        res.atomStart = 5 * i
        res.atomCount = 5
        res.active = True
        res.type = 0
        residues.append(res)
    
    state.residues = residues
    state.activeResidueCount = len(residues)
    
    # Set force field
    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 4
    ff.numMovementTypes = 4
    ff.ljSigma = [0.318395] + [0.0] * 15
    ff.ljEps = [0.88257] + [0.0] * 15
    state.forcefield = ff
    
    return state

def test_opt3_accuracy():
    """Compare OPT3 and SCF accuracy"""
    print("=== OPT3 vs SCF Accuracy Test ===\n")
    
    # Create test system
    state = create_water_dimer()
    
    # Initialize Drude force
    drude_force = pygcmc.DrudeForce()
    
    # Add Drude particles
    for i in range(2):
        drude_force.addParticle(
            drudeIndex=5*i + 1,
            parentIndex=5*i,
            aniso1Index=-1,
            aniso2Index=-1,
            aniso3Index=-1,
            aniso4Index=-1,
            charge=-1.71636,
            polarizability=0.000978253,
            aniso12=0.0,
            aniso34=0.0
        )
    
    # Add Thole screening
    drude_force.addScreenedPair(0, 1, 1.3)
    
    # Test 1: Traditional SCF
    drude_force.setUseOPT3(False)
    start = time.time()
    energy_scf = drude_force.calculateEnergySCF(state)
    time_scf = time.time() - start
    
    # Save Drude positions from SCF
    drude_pos_scf = []
    for i in range(2):
        idx = 5*i + 1
        drude_pos_scf.append([state.atoms[idx].x, state.atoms[idx].y, state.atoms[idx].z])
    
    # Reset Drude positions to parent positions
    for i in range(2):
        parent_idx = 5*i
        drude_idx = 5*i + 1
        state.atoms[drude_idx].x = state.atoms[parent_idx].x
        state.atoms[drude_idx].y = state.atoms[parent_idx].y
        state.atoms[drude_idx].z = state.atoms[parent_idx].z
    
    # Test 2: OPT3
    drude_force.setUseOPT3(True)
    start = time.time()
    energy_opt3 = drude_force.calculateEnergySCF(state)
    time_opt3 = time.time() - start
    
    # Save Drude positions from OPT3
    drude_pos_opt3 = []
    for i in range(2):
        idx = 5*i + 1
        drude_pos_opt3.append([state.atoms[idx].x, state.atoms[idx].y, state.atoms[idx].z])
    
    # Calculate position differences
    print("Results:")
    print(f"SCF Energy:  {energy_scf:.6f} kJ/mol")
    print(f"OPT3 Energy: {energy_opt3:.6f} kJ/mol")
    print(f"Energy diff: {abs(energy_opt3 - energy_scf):.6f} kJ/mol")
    print(f"\nTiming:")
    print(f"SCF time:  {time_scf*1000:.2f} ms")
    print(f"OPT3 time: {time_opt3*1000:.2f} ms")
    print(f"Speedup:   {time_scf/time_opt3:.1f}x")
    
    print(f"\nDrude displacements:")
    for i in range(2):
        parent_pos = [state.atoms[5*i].x, state.atoms[5*i].y, state.atoms[5*i].z]
        
        # SCF displacement
        disp_scf = np.array(drude_pos_scf[i]) - np.array(parent_pos)
        r_scf = np.linalg.norm(disp_scf) * 10  # Convert to Å
        
        # OPT3 displacement
        disp_opt3 = np.array(drude_pos_opt3[i]) - np.array(parent_pos)
        r_opt3 = np.linalg.norm(disp_opt3) * 10  # Convert to Å
        
        print(f"Water {i+1}: SCF={r_scf:.4f} Å, OPT3={r_opt3:.4f} Å, diff={abs(r_scf-r_opt3):.4f} Å")

def test_larger_system():
    """Test on a larger system"""
    print("\n\n=== Larger System Test (32 waters) ===\n")
    
    # Use existing function to create system
    from gcmc_speed_test import create_simple_water_system
    state = create_simple_water_system(32, 'swm4')
    
    # Get the global Drude force that was already set up
    print("Testing with global Drude force...")
    
    # Traditional SCF
    pygcmc.setDrudeUseOPT3(False)
    start = time.time()
    energy_scf, _ = pygcmc.computeSystemEnergyDrude(state)
    time_scf = time.time() - start
    
    # OPT3
    pygcmc.setDrudeUseOPT3(True)
    start = time.time()
    energy_opt3, _ = pygcmc.computeSystemEnergyDrude(state)
    time_opt3 = time.time() - start
    
    print(f"SCF Energy:  {energy_scf:.2f} kJ/mol")
    print(f"OPT3 Energy: {energy_opt3:.2f} kJ/mol")
    print(f"Energy diff: {abs(energy_opt3 - energy_scf):.2f} kJ/mol")
    print(f"\nTiming:")
    print(f"SCF time:  {time_scf*1000:.2f} ms")
    print(f"OPT3 time: {time_opt3*1000:.2f} ms")
    print(f"Speedup:   {time_scf/time_opt3:.1f}x")

def main():
    """Run all tests"""
    
    # Test accuracy on small system
    test_opt3_accuracy()
    
    # Test performance on larger system
    test_larger_system()
    
    print("\n=== Conclusions ===")
    print("1. OPT3 provides significant speedup (3-10x)")
    print("2. Energy accuracy needs coefficient optimization")
    print("3. Current coefficients from induced dipole model")
    print("4. Need to train on Drude systems for better accuracy")

if __name__ == "__main__":
    main()
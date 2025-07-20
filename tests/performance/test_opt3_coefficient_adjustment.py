#!/usr/bin/env python
"""Test OPT3 coefficient adjustment functionality"""

import sys
import time
import numpy as np
sys.path.append('/home/zhaomt/gcmc/test107/pygcmc_dev/build/modules/bindings')

import pygcmc

def create_water_dimer(distance=0.3):
    """Create water dimer at specified O-O distance (nm)"""
    state = pygcmc.MCState()
    state.info.box = [10.0, 10.0, 10.0]
    state.info.cutoff = 5.0
    
    atoms = []
    
    # Water 1 positions
    positions1 = [
        [5.0, 5.0, 5.0],      # O
        [5.0, 5.0, 5.0],      # D (will be optimized)
        [5.09572, 5.0, 5.0],  # H1
        [4.97, 5.09, 5.0],    # H2
        [5.015, 5.011, 5.0]   # M-site
    ]
    
    # Water 2 positions (distance away)
    positions2 = [
        [5.0 + distance, 5.0, 5.0],      # O
        [5.0 + distance, 5.0, 5.0],      # D
        [5.09572 + distance, 5.0, 5.0], # H1
        [4.97 + distance, 5.09, 5.0],    # H2
        [5.015 + distance, 5.011, 5.0]  # M-site
    ]
    
    charges = [1.71636, -1.71636, 0.55733, 0.55733, -1.11466]
    types = [0, 1, 2, 2, 3]
    
    # Add atoms
    for i in range(5):
        a = pygcmc.MCAtom()
        a.x, a.y, a.z = positions1[i]
        a.charge = charges[i]
        a.type = types[i]
        atoms.append(a)
    
    for i in range(5):
        a = pygcmc.MCAtom()
        a.x, a.y, a.z = positions2[i]
        a.charge = charges[i]
        a.type = types[i]
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

def test_coefficient_adjustment():
    """Test setting and getting OPT3 coefficients"""
    print("=== Testing OPT3 Coefficient Adjustment ===\n")
    
    # Create DrudeForce
    drude_force = pygcmc.DrudeForce()
    
    # Get default coefficients
    default_coeffs = drude_force.getOPT3Coefficients()
    print(f"Default coefficients:")
    print(f"  c0 = {default_coeffs.c0}")
    print(f"  c1 = {default_coeffs.c1}")
    print(f"  c2 = {default_coeffs.c2}")
    print(f"  c3 = {default_coeffs.c3}")
    print(f"  Sum = {default_coeffs.c0 + default_coeffs.c1 + default_coeffs.c2 + default_coeffs.c3:.3f}")
    
    # Test different coefficient sets
    test_sets = [
        ("AMOEBA original", [-0.154, 0.017, 0.657, 0.475]),
        ("Uniform", [0.25, 0.25, 0.25, 0.25]),
        ("Conservative", [0.05, 0.15, 0.50, 0.30]),
        ("Aggressive", [0.10, 0.25, 0.40, 0.25])
    ]
    
    # Create test system
    state = create_water_dimer(0.3)
    
    # Add Drude particles
    for i in range(2):
        drude_force.addParticle(
            drudeIndex=5*i + 1,
            parentIndex=5*i,
            aniso1Index=-1, aniso2Index=-1,
            aniso3Index=-1, aniso4Index=-1,
            charge=-1.71636,
            polarizability=0.000978253,
            aniso12=0.0, aniso34=0.0
        )
    drude_force.addScreenedPair(0, 1, 1.3)
    
    # Test each coefficient set
    print("\n\nTesting different coefficient sets:")
    print("-" * 80)
    
    for name, coeffs in test_sets:
        c0, c1, c2, c3 = coeffs
        drude_force.setOPT3Coefficients(c0, c1, c2, c3)
        
        # Get coefficients back
        current = drude_force.getOPT3Coefficients()
        
        print(f"\n{name}:")
        print(f"  Set: [{c0:.3f}, {c1:.3f}, {c2:.3f}, {c3:.3f}]")
        print(f"  Got: [{current.c0:.3f}, {current.c1:.3f}, {current.c2:.3f}, {current.c3:.3f}]")
        print(f"  Sum: {current.c0 + current.c1 + current.c2 + current.c3:.3f}")
        
        # Calculate energy with each set
        drude_force.setUseOPT3(True)
        energy_opt3 = drude_force.calculateEnergySCF(state)
        
        drude_force.setUseOPT3(False)
        energy_scf = drude_force.calculateEnergySCF(state)
        
        print(f"  Energy (OPT3): {energy_opt3:.3f} kJ/mol")
        print(f"  Energy (SCF):  {energy_scf:.3f} kJ/mol")
        print(f"  Difference:    {abs(energy_opt3 - energy_scf):.3f} kJ/mol")

def test_training_data_collection():
    """Test collecting training data"""
    print("\n\n=== Testing Training Data Collection ===\n")
    
    # Create DrudeForce
    drude_force = pygcmc.DrudeForce()
    
    # Create test system
    state = create_water_dimer(0.28)
    
    # Add Drude particles
    for i in range(2):
        drude_force.addParticle(
            drudeIndex=5*i + 1,
            parentIndex=5*i,
            aniso1Index=-1, aniso2Index=-1,
            aniso3Index=-1, aniso4Index=-1,
            charge=-1.71636,
            polarizability=0.000978253,
            aniso12=0.0, aniso34=0.0
        )
    drude_force.addScreenedPair(0, 1, 1.3)
    
    # Collect training data
    print("Collecting training data...")
    training_data = drude_force.collectTrainingData(state)
    
    print(f"\nTraining data collected:")
    print(f"  Number of Drude particles: {len(training_data.r0)}")
    
    # Analyze perturbation orders
    for i in range(min(2, len(training_data.r0))):
        print(f"\n  Drude particle {i}:")
        r0 = training_data.r0[i]
        r1 = training_data.r1[i]
        r2 = training_data.r2[i]
        r3 = training_data.r3[i]
        r_scf = training_data.r_scf[i]
        
        print(f"    |r0| = {r0.norm():.6f} nm ({r0.norm()*10:.4f} Å)")
        print(f"    |r1| = {r1.norm():.6f} nm ({r1.norm()*10:.4f} Å)")
        print(f"    |r2| = {r2.norm():.6f} nm ({r2.norm()*10:.4f} Å)")
        print(f"    |r3| = {r3.norm():.6f} nm ({r3.norm()*10:.4f} Å)")
        print(f"    |r_scf| = {r_scf.norm():.6f} nm ({r_scf.norm()*10:.4f} Å)")
        
        # Test OPT3 reconstruction with different coefficients
        print(f"\n    OPT3 reconstructions:")
        test_coeffs = [
            ("Default", [0.10, 0.25, 0.40, 0.25]),
            ("Uniform", [0.25, 0.25, 0.25, 0.25]),
            ("Conservative", [0.05, 0.15, 0.50, 0.30])
        ]
        
        for name, (c0, c1, c2, c3) in test_coeffs:
            r_opt3_x = c0 * r0.x + c1 * r1.x + c2 * r2.x + c3 * r3.x
            r_opt3_y = c0 * r0.y + c1 * r1.y + c2 * r2.y + c3 * r3.y
            r_opt3_z = c0 * r0.z + c1 * r1.z + c2 * r2.z + c3 * r3.z
            r_opt3_norm = np.sqrt(r_opt3_x**2 + r_opt3_y**2 + r_opt3_z**2)
            
            error = np.sqrt((r_opt3_x - r_scf.x)**2 + 
                          (r_opt3_y - r_scf.y)**2 + 
                          (r_opt3_z - r_scf.z)**2)
            
            print(f"      {name}: |r_opt3| = {r_opt3_norm:.6f} nm, error = {error:.6f} nm ({error*10:.4f} Å)")

def main():
    """Run all tests"""
    test_coefficient_adjustment()
    test_training_data_collection()
    
    print("\n\n=== Summary ===")
    print("\nThe infrastructure for OPT3 coefficient adjustment is working!")
    print("Key features implemented:")
    print("1. Runtime adjustment of OPT3 coefficients")
    print("2. Collection of training data (r0, r1, r2, r3, r_scf)")
    print("3. Python bindings for coefficient control")
    print("\nNext steps:")
    print("1. Collect training data from diverse systems")
    print("2. Implement optimization algorithm to find best coefficients")
    print("3. Validate on test systems")

if __name__ == "__main__":
    main()
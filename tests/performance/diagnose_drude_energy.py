#!/usr/bin/env python
"""Deep diagnosis of Drude energy anomaly"""

import sys
import numpy as np
import math
sys.path.append('/home/zhaomt/gcmc/test107/pygcmc_dev/build/modules/bindings')

import pygcmc

def analyze_single_water():
    """Analyze a single water molecule to understand the energy problem"""
    print("=== Single Water Analysis ===\n")
    
    # Create state with single water
    state = pygcmc.MCState()
    state.info.box = [10.0, 10.0, 10.0]  # Large box
    state.info.cutoff = 5.0
    
    # Initialize Drude force
    pygcmc.initializeDrudeForce()
    
    # Create single water at origin
    atoms = []
    
    # SWM4-NDP parameters
    qO = 1.71636
    qD = -1.71636
    qH = 0.55733
    qM = -1.11466
    rOH = 0.09572  # nm
    aHOH = 104.52 * math.pi / 180  # radians
    
    # Oxygen
    o = pygcmc.MCAtom()
    o.x, o.y, o.z = 5.0, 5.0, 5.0  # Center of box
    o.charge = qO
    o.type = 0
    atoms.append(o)
    
    # Drude - test different initial positions
    d = pygcmc.MCAtom()
    d.x, d.y, d.z = 5.0, 5.0, 5.0  # Start at parent
    d.charge = qD
    d.type = 1
    atoms.append(d)
    
    # Hydrogen 1
    h1 = pygcmc.MCAtom()
    h1.x = 5.0 + rOH
    h1.y = 5.0
    h1.z = 5.0
    h1.charge = qH
    h1.type = 2
    atoms.append(h1)
    
    # Hydrogen 2
    h2 = pygcmc.MCAtom()
    h2.x = 5.0 + rOH * math.cos(aHOH)
    h2.y = 5.0 + rOH * math.sin(aHOH)
    h2.z = 5.0
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
    
    state.atoms = atoms
    state.activeAtomCount = len(atoms)
    
    # Create residue for the water molecule
    res = pygcmc.MCResidue()
    res.atomStart = 0
    res.atomCount = 5
    res.active = True
    res.type = 0
    state.residues = [res]
    state.activeResidueCount = 1
    
    print("Initial atom positions:")
    for i, atom in enumerate(atoms):
        atom_type = ["O", "D", "H1", "H2", "M"][i]
        print(f"  {atom_type}: ({atom.x:.4f}, {atom.y:.4f}, {atom.z:.4f}) q={atom.charge:.4f}")
    
    # Add Drude particle
    pygcmc.addDrudeParticle(
        drudeIndex=1,
        parentIndex=0,
        charge=qD,
        polarizability=0.000978253  # nm^3
    )
    
    print(f"\nNumber of Drude particles: {pygcmc.getNumDrudeParticles()}")
    
    # Test energy calculation with different Drude positions
    print("\n--- Testing different Drude displacements ---")
    
    displacements = [0.0, 0.0001, 0.001, 0.01, 0.1]  # nm
    
    for disp in displacements:
        # Move Drude
        state.atoms[1].x = 5.0 + disp
        
        # Calculate energy
        result = pygcmc.computeSystemEnergyDrude(state)
        
        if isinstance(result, tuple):
            energy, components = result
            print(f"\nDrude displacement: {disp} nm")
            print(f"  Total energy: {energy:.2f} kJ/mol")
            print(f"  Components: {components}")
        else:
            print(f"\nDrude displacement: {disp} nm")
            print(f"  Total energy: {result:.2f} kJ/mol")
    
    # Now test SCF convergence
    print("\n--- Testing SCF Convergence ---")
    
    # Reset Drude to a displaced position
    state.atoms[1].x = 5.01  # 0.01 nm displacement
    
    # Set different SCF parameters
    scf_params = [
        {'tolerance': 1.0, 'max_iter': 50},
        {'tolerance': 10.0, 'max_iter': 100},
        {'tolerance': 100.0, 'max_iter': 200},
    ]
    
    for params in scf_params:
        scf_obj = pygcmc.DrudeSCFParams()
        scf_obj.tolerance = params['tolerance']
        scf_obj.maxIterations = params['max_iter']
        pygcmc.setDrudeSCFParameters(scf_obj)
        
        print(f"\nSCF with tolerance={params['tolerance']}, max_iter={params['max_iter']}")
        
        # Calculate energy (which triggers SCF)
        result = pygcmc.computeSystemEnergyDrude(state)
        
        if isinstance(result, tuple):
            energy, components = result
            print(f"  Energy: {energy:.2f} kJ/mol")
            print(f"  Final Drude position: ({state.atoms[1].x:.6f}, {state.atoms[1].y:.6f}, {state.atoms[1].z:.6f})")
            print(f"  Drude displacement: {state.atoms[1].x - state.atoms[0].x:.6f} nm")
        else:
            energy = result
            print(f"  Energy: {energy:.2f} kJ/mol")

def analyze_force_calculation():
    """Analyze force calculation in detail"""
    print("\n\n=== Force Calculation Analysis ===\n")
    
    # Constants
    ONE_4PI_EPS0 = 138.935456  # kJ*nm/mol/e^2
    k_drude = 418400.0  # kJ/mol/nm^2 (1000 kcal/mol/A^2)
    
    # Charges
    qD = -1.71636
    qH = 0.55733
    qM = -1.11466
    
    # Typical distances in SWM4-NDP water
    r_DH = 0.11  # nm (Drude to Hydrogen)
    r_DM = 0.02  # nm (Drude to M-site, very close!)
    
    print("Force contributions on Drude particle:")
    
    # 1. Harmonic force (assuming 0.01 nm displacement)
    disp = 0.01  # nm
    F_harm = -k_drude * disp
    print(f"\n1. Harmonic force (disp={disp} nm):")
    print(f"   F = -k * r = -{k_drude} * {disp} = {F_harm:.2f} kJ/mol/nm")
    
    # 2. Coulomb from H1
    F_DH = ONE_4PI_EPS0 * qD * qH / (r_DH**2)
    print(f"\n2. Coulomb from H (r={r_DH} nm):")
    print(f"   F = k*q1*q2/r^2 = {ONE_4PI_EPS0} * {qD} * {qH} / {r_DH}^2 = {F_DH:.2f} kJ/mol/nm")
    
    # 3. Coulomb from M
    F_DM = ONE_4PI_EPS0 * qD * qM / (r_DM**2)
    print(f"\n3. Coulomb from M (r={r_DM} nm):")
    print(f"   F = k*q1*q2/r^2 = {ONE_4PI_EPS0} * {qD} * {qM} / {r_DM}^2 = {F_DM:.2f} kJ/mol/nm")
    
    print(f"\nTotal force magnitude: ~{abs(F_harm + 2*F_DH + F_DM):.0f} kJ/mol/nm")
    
    print("\n*** Key insight: D-M interaction creates HUGE repulsive force! ***")
    print(f"    Both D and M have negative charges ({qD}, {qM})")
    print(f"    At close distance ({r_DM} nm), force is {F_DM:.0f} kJ/mol/nm")
    print("    This pushes Drude away from equilibrium position")

def propose_solutions():
    """Propose solutions to the energy problem"""
    print("\n\n=== Proposed Solutions ===\n")
    
    print("1. **Intramolecular Exclusions**")
    print("   - Must exclude ALL intramolecular non-bonded interactions")
    print("   - Only keep: Drude-parent harmonic, and Drude with other molecules")
    print("   - Current issue: D-H and D-M intramolecular forces are included")
    
    print("\n2. **SCF Algorithm Improvements**")
    print("   - Use smaller initial Drude displacement (1e-6 nm instead of 0)")
    print("   - Implement adaptive damping based on force magnitude")
    print("   - Add convergence checks based on position change, not just force")
    
    print("\n3. **Force Constant Adjustment**")
    print("   - Current k=418400 kJ/mol/nm^2 may be too stiff")
    print("   - Consider using polarizability-based force constant")
    print("   - k = q^2 / (4πε₀ * α) where α is polarizability")
    
    print("\n4. **Initialization Strategy**")
    print("   - Start with Drude exactly at parent position")
    print("   - Do preliminary SCF with only intermolecular forces")
    print("   - Then include full system")
    
    print("\n5. **Debugging Steps**")
    print("   - Print force components during SCF")
    print("   - Monitor Drude displacement at each iteration")
    print("   - Check which interactions dominate the energy")

def main():
    """Run all diagnostic tests"""
    analyze_single_water()
    analyze_force_calculation()
    propose_solutions()
    
    print("\n\n=== Summary ===")
    print("The main problem is intramolecular D-H and D-M interactions")
    print("These should be excluded in SWM4-NDP model")
    print("The huge negative energy comes from incorrect force calculations")
    print("This pushes Drude particles far from equilibrium")

if __name__ == "__main__":
    main()
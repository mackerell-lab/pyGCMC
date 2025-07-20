#!/usr/bin/env python
"""Complete energy calculation for SWM4-NDP water"""

import sys
import numpy as np
sys.path.append('/home/zhaomt/gcmc/test107/pygcmc_dev/build/modules/bindings')

import pygcmc

def calculate_coulomb_energy(state):
    """Calculate all intermolecular Coulomb interactions"""
    ONE_4PI_EPS0 = 138.935456  # kJ·nm/mol/e²
    energy = 0.0
    
    for res1_idx in range(state.activeResidueCount):
        res1 = state.residues[res1_idx]
        for res2_idx in range(res1_idx + 1, state.activeResidueCount):
            res2 = state.residues[res2_idx]
            
            # Calculate all pairwise interactions between residues
            for i in range(res1.atomStart, res1.atomStart + res1.atomCount):
                for j in range(res2.atomStart, res2.atomStart + res2.atomCount):
                    atom1 = state.atoms[i]
                    atom2 = state.atoms[j]
                    
                    r_vec = np.array([atom2.x - atom1.x, 
                                     atom2.y - atom1.y, 
                                     atom2.z - atom1.z])
                    r = np.linalg.norm(r_vec)
                    
                    if r > 1e-10:
                        energy += ONE_4PI_EPS0 * atom1.charge * atom2.charge / r
    
    return energy

def test_water_complete_energy():
    print("=== Complete Energy Calculation for SWM4-NDP Water ===\n")
    
    # Initialize
    pygcmc.initializeDrudeForce()
    
    # Create state
    state = pygcmc.MCState()
    state.info.cutoff = 2.0
    state.info.box = [10.0, 10.0, 10.0]
    
    # Parameters from PSF
    qO_core = 1.66260
    qD = -1.76260
    qM = -0.95710
    qH = 0.52855
    alpha_nm3 = 0.0013
    k_kj_nm2 = 166018.9
    
    # Test at 3 Å separation
    separation = 0.3  # nm
    
    # Create two waters
    atoms = []
    
    # Water 1
    o1 = pygcmc.MCAtom()
    o1.x = 0.0
    o1.y = 0.0
    o1.z = 0.0
    o1.charge = qO_core
    o1.type = 0
    atoms.append(o1)
    
    d1 = pygcmc.MCAtom()
    d1.x = 0.0
    d1.y = 0.0
    d1.z = 0.0
    d1.charge = qD
    d1.type = 1
    atoms.append(d1)
    
    h1_1 = pygcmc.MCAtom()
    h1_1.x = 0.09572
    h1_1.y = 0.0
    h1_1.z = 0.0
    h1_1.charge = qH
    h1_1.type = 2
    atoms.append(h1_1)
    
    h2_1 = pygcmc.MCAtom()
    h2_1.x = -0.023999
    h2_1.y = 0.092663
    h2_1.z = 0.0
    h2_1.charge = qH
    h2_1.type = 2
    atoms.append(h2_1)
    
    m1 = pygcmc.MCAtom()
    m1.x = 0.00793
    m1.y = 0.00986
    m1.z = 0.0
    m1.charge = qM
    m1.type = 3
    atoms.append(m1)
    
    # Water 2
    o2 = pygcmc.MCAtom()
    o2.x = separation
    o2.y = 0.0
    o2.z = 0.0
    o2.charge = qO_core
    o2.type = 0
    atoms.append(o2)
    
    d2 = pygcmc.MCAtom()
    d2.x = separation
    d2.y = 0.0
    d2.z = 0.0
    d2.charge = qD
    d2.type = 1
    atoms.append(d2)
    
    h1_2 = pygcmc.MCAtom()
    h1_2.x = separation + 0.09572
    h1_2.y = 0.0
    h1_2.z = 0.0
    h1_2.charge = qH
    h1_2.type = 2
    atoms.append(h1_2)
    
    h2_2 = pygcmc.MCAtom()
    h2_2.x = separation - 0.023999
    h2_2.y = 0.092663
    h2_2.z = 0.0
    h2_2.charge = qH
    h2_2.type = 2
    atoms.append(h2_2)
    
    m2 = pygcmc.MCAtom()
    m2.x = separation + 0.00793
    m2.y = 0.00986
    m2.z = 0.0
    m2.charge = qM
    m2.type = 3
    atoms.append(m2)
    
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
    
    # Force field
    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 4
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
    
    print("Before SCF:")
    print("-" * 40)
    
    # Calculate initial Coulomb energy
    coulomb_initial = calculate_coulomb_energy(state)
    print(f"Coulomb energy (D at parent): {coulomb_initial:.2f} kJ/mol")
    
    # Run SCF
    drude_energy, energy_dict = pygcmc.computeSystemEnergyDrude(state)
    
    print("\nAfter SCF:")
    print("-" * 40)
    
    # Calculate final Coulomb energy
    coulomb_final = calculate_coulomb_energy(state)
    print(f"Coulomb energy (D moved): {coulomb_final:.2f} kJ/mol")
    print(f"Drude harmonic energy: {drude_energy:.2f} kJ/mol")
    print(f"Total energy: {coulomb_final + drude_energy:.2f} kJ/mol")
    
    # Show Drude displacements
    print("\nDrude displacements:")
    for i, drude_idx in enumerate([1, 6]):
        parent_idx = [0, 5][i]
        dx = state.atoms[drude_idx].x - state.atoms[parent_idx].x
        dy = state.atoms[drude_idx].y - state.atoms[parent_idx].y
        dz = state.atoms[drude_idx].z - state.atoms[parent_idx].z
        dist = np.sqrt(dx*dx + dy*dy + dz*dz)
        print(f"  Water {i+1}: {dist*1000:.3f} pm")
    
    print("\nEnergy breakdown:")
    print("-" * 40)
    print(f"Initial Coulomb (no polarization): {coulomb_initial:.2f} kJ/mol")
    print(f"Polarization energy: {coulomb_final - coulomb_initial:.2f} kJ/mol")
    print(f"Drude restraint energy: {drude_energy:.2f} kJ/mol")
    print(f"Net interaction energy: {coulomb_final + drude_energy:.2f} kJ/mol")

if __name__ == "__main__":
    test_water_complete_energy()
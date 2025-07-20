#!/usr/bin/env python
"""Detailed analysis of water dimer interaction"""

import sys
import numpy as np
sys.path.append('/home/zhaomt/gcmc/test107/pygcmc_dev/build/modules/bindings')

import pygcmc

def test_water_dimer_detailed():
    print("=== Detailed Water Dimer Analysis ===\n")
    
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
    
    # Create water molecules at different distances
    distances = [0.25, 0.3, 0.35, 0.4, 0.5, 0.6]  # nm
    
    for dist in distances:
        # Reset state
        state.atoms = []
        state.residues = []
        
        # Water 1
        atoms = []
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
        
        # Water 2 - displaced along x
        o2 = pygcmc.MCAtom()
        o2.x = dist
        o2.y = 0.0
        o2.z = 0.0
        o2.charge = qO_core
        o2.type = 0
        atoms.append(o2)
        
        d2 = pygcmc.MCAtom()
        d2.x = dist
        d2.y = 0.0
        d2.z = 0.0
        d2.charge = qD
        d2.type = 1
        atoms.append(d2)
        
        h1_2 = pygcmc.MCAtom()
        h1_2.x = dist + 0.09572
        h1_2.y = 0.0
        h1_2.z = 0.0
        h1_2.charge = qH
        h1_2.type = 2
        atoms.append(h1_2)
        
        h2_2 = pygcmc.MCAtom()
        h2_2.x = dist - 0.023999
        h2_2.y = 0.092663
        h2_2.z = 0.0
        h2_2.charge = qH
        h2_2.type = 2
        atoms.append(h2_2)
        
        m2 = pygcmc.MCAtom()
        m2.x = dist + 0.00793
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
        
        # Clear and re-add Drude particles
        pygcmc.clearDrudeForce()
        pygcmc.initializeDrudeForce()
        
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
        
        # Calculate energy
        energy, converged = pygcmc.computeSystemEnergyDrude(state)
        
        # Calculate some key distances
        r_OO = dist
        r_OM1 = np.sqrt((dist - 0.00793)**2 + 0.00986**2)
        r_MO2 = np.sqrt((dist + 0.00793)**2 + 0.00986**2)
        
        # Calculate expected dominant contributions
        ONE_4PI_EPS0 = 138.935456
        qO_net = qO_core + qD  # Net charge on oxygen
        
        # O-O interaction (net charges)
        E_OO = ONE_4PI_EPS0 * qO_net * qO_net / r_OO
        
        # O-M interactions
        E_OM1 = ONE_4PI_EPS0 * qO_net * qM / r_OM1
        E_MO2 = ONE_4PI_EPS0 * qO_net * qM / r_MO2
        
        # M-M interaction
        r_MM = dist
        E_MM = ONE_4PI_EPS0 * qM * qM / r_MM
        
        print(f"Distance: {dist*10:.1f} Å")
        print(f"  Total energy: {energy:.2f} kJ/mol")
        print(f"  Converged: {converged}")
        print(f"  Key contributions (estimate):")
        print(f"    O-O (net): {E_OO:.2f} kJ/mol")
        print(f"    O-M cross: {E_OM1 + E_MO2:.2f} kJ/mol")
        print(f"    M-M: {E_MM:.2f} kJ/mol")
        print(f"    Sum: {E_OO + E_OM1 + E_MO2 + E_MM:.2f} kJ/mol")
        print()
    
    print("\nAnalysis:")
    print("-" * 50)
    print("The interaction energy includes:")
    print("1. Direct Coulomb between all intermolecular atom pairs")
    print("2. Polarization effects through Drude displacement")
    print("3. Virtual site M participates in intermolecular interactions")
    print("\nNote: Real water dimer also has:")
    print("- H-bonding directionality")
    print("- van der Waals interactions")
    print("- Many-body polarization effects")

if __name__ == "__main__":
    test_water_dimer_detailed()
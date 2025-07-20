#!/usr/bin/env python
"""Test Drude hard wall constraint - minimal version"""

import sys
import numpy as np
sys.path.append('/home/zhaomt/gcmc/test107/pygcmc_dev/build/modules/bindings')

import pygcmc

def test_drude_hardwall():
    print("=== Testing Drude Hard Wall Constraint ===\n")
    
    # Initialize
    pygcmc.initializeDrudeForce()
    
    # Create a simple system
    state = pygcmc.MCState()
    state.info.cutoff = 1.0  # 10 Å cutoff
    state.info.box = [5.0, 5.0, 5.0]  # 50 Å box
    
    # Parent atom
    parent = pygcmc.MCAtom()
    parent.x = 0.0
    parent.y = 0.0
    parent.z = 0.0
    parent.charge = 1.0
    parent.type = 0
    
    # Drude particle (start at parent)
    drude = pygcmc.MCAtom()
    drude.x = 0.0
    drude.y = 0.0
    drude.z = 0.0
    drude.charge = -1.0
    drude.type = 1
    
    # External test charge to pull Drude
    external = pygcmc.MCAtom()
    external.x = 0.2  # 2 Å away
    external.y = 0.0
    external.z = 0.0
    external.charge = 5.0  # Strong positive charge
    external.type = 2
    
    # Set atoms
    state.atoms = [parent, drude, external]
    state.activeAtomCount = 3
    
    # Force field (minimal)
    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 3
    state.forcefield = ff
    
    # Set up Drude particle
    # k = 1000 kcal/mol/Å² = 418400 kJ/mol/nm²
    polarizability = 138.935456 * 1.0 * 1.0 / 418400.0
    
    pygcmc.addDrudeParticle(
        drudeIndex=1,
        parentIndex=0,
        charge=-1.0,
        polarizability=polarizability
    )
    
    print("Setup: Parent at origin, external charge at 2 Å")
    print("External charge will pull Drude away from parent\n")
    
    # Test different hard wall distances
    for max_dist in [0.005, 0.01, 0.02, 0.05]:  # nm
        print(f"\nHard wall at {max_dist} nm ({max_dist*10} Å):")
        
        # Reset Drude to parent
        state.atoms[1].x = 0.0
        
        # Set SCF parameters
        params = pygcmc.DrudeSCFParams()
        params.tolerance = 1e-8
        params.maxIterations = 100
        params.dampingFactor = 0.3
        params.maxDrudeDistance = max_dist
        pygcmc.setDrudeSCFParameters(params)
        
        # Calculate energy
        energy, _ = pygcmc.computeSystemEnergyDrude(state)
        
        # Check Drude position
        dx = state.atoms[1].x - state.atoms[0].x
        dy = state.atoms[1].y - state.atoms[0].y
        dz = state.atoms[1].z - state.atoms[0].z
        dist = np.sqrt(dx*dx + dy*dy + dz*dz)
        
        print(f"  Drude displacement: {dist:.6f} nm ({dist*10:.4f} Å)")
        print(f"  Energy: {energy:.2f} kJ/mol")
        
        # Check if hard wall is active
        if abs(dist - max_dist) < 1e-6:
            print(f"  ✓ Hard wall constraint ACTIVE")
        else:
            print(f"  - Hard wall not reached")
    
    # Test without hard wall
    print("\n\nWithout hard wall (max_dist = 1.0 nm):")
    params.maxDrudeDistance = 1.0
    pygcmc.setDrudeSCFParameters(params)
    
    state.atoms[1].x = 0.0  # Reset
    energy, _ = pygcmc.computeSystemEnergyDrude(state)
    
    dx = state.atoms[1].x - state.atoms[0].x
    dist = abs(dx)
    
    print(f"  Drude displacement: {dist:.6f} nm ({dist*10:.4f} Å)")
    print(f"  Energy: {energy:.2f} kJ/mol")
    print(f"  Much larger than 0.02 nm OpenMM default!")

if __name__ == "__main__":
    test_drude_hardwall()
#!/usr/bin/env python
"""Test Drude implementation with hard wall constraint - simplified version"""

import sys
import numpy as np
sys.path.append('/home/zhaomt/gcmc/test107/pygcmc_dev/build/modules/bindings')

import pygcmc

def test_drude_hardwall():
    print("=== Testing Drude with Hard Wall Constraint ===\n")
    
    # Initialize
    pygcmc.initializeDrudeForce()
    
    # Create a simple system
    state = pygcmc.MCState()
    state.info.cutoff = 1.0  # 10 Å cutoff
    state.info.box = [5.0, 5.0, 5.0]  # 50 Å box
    
    # Create atoms list
    atoms = []
    
    # Parent atom
    parent = pygcmc.MCAtom()
    parent.x = 0.0
    parent.y = 0.0
    parent.z = 0.0
    parent.charge = 1.0  # Positive charge
    parent.sigma = 0.3
    parent.epsilon = 0.1
    parent.type = 0
    atoms.append(parent)
    
    # Drude particle (start displaced)
    drude = pygcmc.MCAtom()
    drude.x = 0.05  # 0.5 Å displacement
    drude.y = 0.0
    drude.z = 0.0
    drude.charge = -1.0  # Negative charge
    drude.sigma = 0.0
    drude.epsilon = 0.0
    drude.type = 1
    atoms.append(drude)
    
    # Set atoms
    state.atoms = atoms
    state.activeAtomCount = 2
    
    # Set up force field
    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 2
    ff.numMovementTypes = 2
    ff.ljSigma = [0.3, 0.0, 0.0, 0.0]
    ff.ljEps = [0.1, 0.0, 0.0, 0.0]
    state.forcefield = ff
    
    # Set up Drude particle
    # k = 1000 kcal/mol/Å² = 418400 kJ/mol/nm²
    # α = q²/(4πε₀k)
    polarizability = 138.935456 * 1.0 * 1.0 / (418400.0)
    
    drude_idx = pygcmc.addDrudeParticle(
        drudeIndex=1,
        parentIndex=0,
        charge=-1.0,
        polarizability=polarizability
    )
    
    print("Initial Drude displacement: 0.05 nm (0.5 Å)\n")
    
    # Test with different hard wall distances
    print("Testing different hard wall constraints:\n")
    
    for max_dist in [0.002, 0.01, 0.02, 0.05]:  # nm
        print(f"Max Drude distance: {max_dist} nm ({max_dist*10} Å)")
        
        # Reset Drude position far from parent
        state.atoms[1].x = 0.05  # 0.5 Å away
        
        # Set SCF parameters with hard wall
        params = pygcmc.DrudeSCFParams()
        params.tolerance = 1e-8
        params.maxIterations = 100
        params.dampingFactor = 0.5
        params.maxDrudeDistance = max_dist
        pygcmc.setDrudeSCFParameters(params)
        
        # Calculate energy
        energy, energy_dict = pygcmc.computeSystemEnergyDrude(state)
        
        # Check final Drude position
        dr = state.atoms[1].x - state.atoms[0].x
        dist = abs(dr)
        
        print(f"  Final displacement: {dist:.6f} nm ({dist*10:.4f} Å)")
        print(f"  Energy: {energy:.2f} kJ/mol")
        print(f"  Constraint active: {'YES' if dist >= max_dist*0.999 else 'NO'}")
        print()
    
    # Test extreme case with external field
    print("\nTesting extreme external field:\n")
    
    # Add very strong external charge
    test_charge = pygcmc.MCAtom()
    test_charge.x = 0.1  # 1 Å away
    test_charge.y = 0.0
    test_charge.z = 0.0
    test_charge.charge = 10.0  # Very strong positive charge
    test_charge.sigma = 0.0
    test_charge.epsilon = 0.0
    test_charge.type = 2
    
    state.atoms.append(test_charge)
    state.activeAtomCount = 3
    
    # Reset Drude to parent
    state.atoms[1].x = 0.0
    
    # Set hard wall at 0.02 nm (OpenMM default)
    params.maxDrudeDistance = 0.02
    pygcmc.setDrudeSCFParameters(params)
    
    energy, _ = pygcmc.computeSystemEnergyDrude(state)
    
    dr = state.atoms[1].x - state.atoms[0].x
    dist = abs(dr)
    
    print(f"With 10e charge at 1 Å:")
    print(f"  Drude displacement: {dist:.6f} nm ({dist*10:.4f} Å)")
    print(f"  Expected: exactly 0.02 nm (hard wall limit)")
    print(f"  Energy: {energy:.2f} kJ/mol")
    print(f"  Hard wall working: {'YES' if abs(dist - 0.02) < 1e-6 else 'NO!!!'}")
    
    # Show that without hard wall, Drude would move further
    print("\n\nFor comparison, without hard wall constraint:")
    params.maxDrudeDistance = 1.0  # Very large limit
    pygcmc.setDrudeSCFParameters(params)
    
    state.atoms[1].x = 0.0  # Reset
    energy, _ = pygcmc.computeSystemEnergyDrude(state)
    
    dr = state.atoms[1].x - state.atoms[0].x
    dist = abs(dr)
    
    print(f"  Drude displacement: {dist:.6f} nm ({dist*10:.4f} Å)")
    print(f"  Much larger than 0.02 nm limit!")

if __name__ == "__main__":
    test_drude_hardwall()
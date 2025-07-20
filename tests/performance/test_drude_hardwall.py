#!/usr/bin/env python
"""Test Drude implementation with hard wall constraint"""

import sys
import numpy as np
sys.path.append('/home/zhaomt/gcmc/test107/pygcmc_dev/build/modules/bindings')

import pygcmc

def test_drude_hardwall():
    print("=== Testing Drude with Hard Wall Constraint ===\n")
    
    # Initialize
    pygcmc.initializeDrudeForce()
    
    # Create a simple water molecule
    state = pygcmc.MCState()
    state.info.cutoff = 1.0  # 10 Å cutoff
    state.info.box = [5.0, 5.0, 5.0]  # 50 Å box
    
    # SWM4-NDP water parameters
    qO = 1.71636
    qD = -1.71636
    qH = 0.55733
    qM = -1.11466
    
    # Add atoms
    # Oxygen at origin
    state.add_atom(pygcmc.MCAtom(
        x=0.0, y=0.0, z=0.0,
        charge=qO,
        sigma=0.318395, epsilon=0.21094
    ))
    
    # Drude on oxygen (start slightly displaced)
    state.add_atom(pygcmc.MCAtom(
        x=0.001, y=0.001, z=0.001,  # 0.01 Å displacement
        charge=qD,
        sigma=0.0, epsilon=0.0
    ))
    
    # Hydrogens
    state.add_atom(pygcmc.MCAtom(
        x=0.09572, y=0.0, z=0.0,
        charge=qH,
        sigma=0.0, epsilon=0.0
    ))
    
    state.add_atom(pygcmc.MCAtom(
        x=-0.02395, y=0.09273, z=0.0,
        charge=qH,
        sigma=0.0, epsilon=0.0
    ))
    
    # M-site
    state.add_atom(pygcmc.MCAtom(
        x=0.00793, y=0.00986, z=0.0,
        charge=qM,
        sigma=0.0, epsilon=0.0
    ))
    
    # Set up Drude particle
    # k = 1000 kcal/mol/Å² = 418400 kJ/mol/nm²
    polarizability = 138.935456 * qD * qD / (418400.0)
    
    drude_idx = pygcmc.addDrudeParticle(
        drudeIndex=1,
        parentIndex=0,
        charge=qD,
        polarizability=polarizability
    )
    
    # Test with different hard wall distances
    print("Testing different hard wall constraints:\n")
    
    for max_dist in [0.002, 0.01, 0.02, 0.05]:  # nm
        print(f"\nMax Drude distance: {max_dist} nm ({max_dist*10} Å)")
        
        # Reset Drude position far from parent
        state.atoms[1].x = 0.05  # 0.5 Å away
        state.atoms[1].y = 0.05
        state.atoms[1].z = 0.05
        
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
        dr = np.array([
            state.atoms[1].x - state.atoms[0].x,
            state.atoms[1].y - state.atoms[0].y,
            state.atoms[1].z - state.atoms[0].z
        ])
        dist = np.linalg.norm(dr)
        
        print(f"  Initial displacement: 0.0866 nm")
        print(f"  Final displacement: {dist:.6f} nm ({dist*10:.4f} Å)")
        print(f"  Energy: {energy:.2f} kJ/mol")
        print(f"  Constraint active: {'YES' if dist >= max_dist*0.999 else 'NO'}")
    
    # Test SCF convergence with hard wall
    print("\n\nTesting SCF convergence with external field:\n")
    
    # Apply external electric field by adding a test charge
    state.add_atom(model.MCAtom(
        x=0.5, y=0.0, z=0.0,  # 5 Å away
        charge=1.0,
        sigma=0.0, epsilon=0.0
    ))
    
    # Reset Drude to parent position
    state.atoms[1].x = state.atoms[0].x
    state.atoms[1].y = state.atoms[0].y
    state.atoms[1].z = state.atoms[0].z
    
    # Set hard wall at 0.02 nm (OpenMM default)
    params.maxDrudeDistance = 0.02
    pygcmc.setDrudeSCFParameters(params)
    
    # Calculate with external field
    energy, _ = pygcmc.computeSystemEnergyDrude(state)
    
    dr = np.array([
        state.atoms[1].x - state.atoms[0].x,
        state.atoms[1].y - state.atoms[0].y,
        state.atoms[1].z - state.atoms[0].z
    ])
    dist = np.linalg.norm(dr)
    
    print(f"With external field:")
    print(f"  Drude displacement: {dist:.6f} nm ({dist*10:.4f} Å)")
    print(f"  Energy: {energy:.2f} kJ/mol")
    print(f"  Hard wall active: {'YES' if dist >= 0.02*0.999 else 'NO'}")
    
    # Test extreme case
    print("\n\nTesting extreme external field:\n")
    
    # Move test charge very close
    state.atoms[5].x = 0.1  # 1 Å away
    state.atoms[5].charge = 10.0  # Very strong
    
    energy, _ = pygcmc.computeSystemEnergyDrude(state)
    
    dr = np.array([
        state.atoms[1].x - state.atoms[0].x,
        state.atoms[1].y - state.atoms[0].y,
        state.atoms[1].z - state.atoms[0].z
    ])
    dist = np.linalg.norm(dr)
    
    print(f"With extreme field (10e at 1 Å):")
    print(f"  Drude displacement: {dist:.6f} nm ({dist*10:.4f} Å)")
    print(f"  Should be exactly: 0.02 nm (hard wall limit)")
    print(f"  Energy: {energy:.2f} kJ/mol")
    print(f"  Constraint working: {'YES' if abs(dist - 0.02) < 1e-6 else 'NO!!!'}")

if __name__ == "__main__":
    test_drude_hardwall()
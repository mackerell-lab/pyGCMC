#!/usr/bin/env python
"""Debug cavity bias calculation"""
import numpy as np
import pygcmc

def test_cavity():
    """Test cavity bias calculation for empty box"""
    state = pygcmc.MCState()
    state.info.box = np.array([4.0, 4.0, 4.0])  # nm
    
    # Setup ideal gas forcefield
    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.maxTypes = 1
    ff.ljEps = [0.0]  # No LJ interactions
    ff.ljSigma = [0.35]  # But valid sigma
    state.forcefield = ff
    
    # Setup movement module
    mover = pygcmc.movement.MovementModule()
    params = pygcmc.movement.MovementParams()
    params.temperature = 298.15
    params.chemicalPotential = -5.0
    params.useCavityBias = True
    params.cavityGridSpacing = 0.25
    params.probeRadius = 0.10
    params.seed = 12345
    mover.setParams(params)
    
    # Try insertion
    print("Empty box test:")
    result = mover.attemptCavityBiasInsertion(state)
    print(f"  accepted: {result.accepted}")
    print(f"  cavityBiasFactor: {result.cavityBiasFactor}")
    
    # Add a particle manually to test non-empty
    if result.accepted:
        print("\nAfter one particle:")
        result2 = mover.attemptCavityBiasInsertion(state)
        print(f"  accepted: {result2.accepted}")
        print(f"  cavityBiasFactor: {result2.cavityBiasFactor}")
        
        # Add more particles
        for i in range(5):
            mover.attemptCavityBiasInsertion(state)
        
        print(f"\nAfter ~6 particles (n={state.activeResidueCount}):")
        result3 = mover.attemptCavityBiasInsertion(state)
        print(f"  cavityBiasFactor: {result3.cavityBiasFactor}")

if __name__ == "__main__":
    test_cavity()
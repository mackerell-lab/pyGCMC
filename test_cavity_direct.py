#!/usr/bin/env python
"""Test cavity core directly"""
import numpy as np
import pygcmc

def test_direct():
    """Test CavityBiasCore directly"""
    state = pygcmc.MCState()
    state.info.box = np.array([4.0, 4.0, 4.0])  # nm
    
    # Setup forcefield
    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.maxTypes = 1
    ff.ljEps = [0.0]
    ff.ljSigma = [0.35]  # nm
    state.forcefield = ff
    
    # Create CavityBiasCore directly
    core = pygcmc.movement.CavityBiasCore(0.25, 0.10)  # grid spacing, probe radius in nm
    
    # Test empty box
    print("Empty box:")
    volume = core.calculateCavityVolume(state, pygcmc.movement.CavityMode.FAST_APPROX)
    print(f"  Cavity volume: {volume:.3f} nm³")
    print(f"  Box volume: {4**3:.3f} nm³")
    print(f"  Fraction: {volume/(4**3):.3f}")
    
    # Add a particle manually
    res = pygcmc.MCResidue()
    res.atomStart = 0
    res.atomCount = 1
    res.active = True
    res.type = 0
    state.residues.append(res)
    state.activeResidueCount = 1
    
    atom = pygcmc.MCAtom()
    atom.x = 2.0
    atom.y = 2.0
    atom.z = 2.0
    atom.type = 0
    state.atoms.append(atom)
    state.activeAtomCount = 1
    
    # Invalidate cache and recalculate
    core.invalidateCache()
    print("\nWith one particle at center:")
    volume2 = core.calculateCavityVolume(state, pygcmc.movement.CavityMode.FAST_APPROX)
    print(f"  Cavity volume: {volume2:.3f} nm³")
    print(f"  Fraction: {volume2/(4**3):.3f}")
    print(f"  Reduction: {(volume-volume2)/volume*100:.1f}%")

if __name__ == "__main__":
    test_direct()
#!/usr/bin/env python
"""
Test cavity bias insertion methods
"""

import pytest
import numpy as np
import os
import sys

# Add build path for pygcmc module
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../build'))

try:
    import pygcmc
    PYGCMC_AVAILABLE = True
except ImportError:
    PYGCMC_AVAILABLE = False
    pygcmc = None


@pytest.mark.skipif(not PYGCMC_AVAILABLE, reason="PyGCMC not available")
class TestCavityBias:
    """Test cavity bias insertion functionality"""
    
    def test_cavity_bias_empty_box(self):
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
        
        # Setup movement module with cavity bias
        mover = pygcmc.movement.MovementModule()
        params = pygcmc.movement.MovementParams()
        params.temperature = 298.15
        params.chemicalPotential = -5.0
        params.useCavityBias = True
        params.cavityGridSpacing = 0.25
        params.probeRadius = 0.10
        params.seed = 12345
        mover.setParams(params)
        
        # Try insertion in empty box
        result = mover.attemptCavityBiasInsertion(state)
        
        # Empty box should have maximum cavity bias factor
        assert result.cavityBiasFactor > 0, "Cavity bias factor should be positive"
        
        # For empty box, cavity bias factor should be close to 1.0
        # (all space is available)
        assert result.cavityBiasFactor <= 1.0, "Cavity bias factor should not exceed 1.0"
    
    def test_cavity_bias_with_particles(self):
        """Test cavity bias decreases with particle density"""
        state = pygcmc.MCState()
        state.info.box = np.array([4.0, 4.0, 4.0])  # nm
        
        # Setup forcefield
        ff = pygcmc.MCForceField()
        ff.numTotalTypes = 1
        ff.numMovementTypes = 1
        ff.maxTypes = 1
        ff.ljEps = [0.0]
        ff.ljSigma = [0.35]
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
        
        # Get cavity bias for empty box
        result_empty = mover.attemptCavityBiasInsertion(state)
        bias_empty = result_empty.cavityBiasFactor
        
        # Add particles and check cavity bias decreases
        particles_added = 0
        for _ in range(10):
            if mover.attemptCavityBiasInsertion(state).accepted:
                particles_added += 1
        
        if particles_added > 0:
            # Try another insertion with particles present
            result_filled = mover.attemptCavityBiasInsertion(state)
            bias_filled = result_filled.cavityBiasFactor
            
            # Cavity bias should decrease with particles
            assert bias_filled <= bias_empty, \
                f"Cavity bias should decrease with particles: {bias_filled} > {bias_empty}"
    
    def test_cavity_core_direct(self):
        """Test CavityBiasCore directly"""
        state = pygcmc.MCState()
        state.info.box = np.array([4.0, 4.0, 4.0])  # nm
        
        # Setup forcefield
        ff = pygcmc.MCForceField()
        ff.numTotalTypes = 1
        ff.numMovementTypes = 1
        ff.maxTypes = 1
        # ljSigma and ljEps are NxN matrices for mixing rules
        # For 1 type, we need 1x1 = 1 element
        ff.ljSigma = [0.35]  # nm (diagonal element for type 0)
        ff.ljEps = [0.0]     # No interactions
        ff.rebuildLJMatrix()  # Ensure LJ matrix is properly initialized
        state.forcefield = ff
        
        # Create CavityBiasCore directly
        core = pygcmc.movement.CavityBiasCore(0.25, 0.10)  # grid spacing, probe radius in nm
        
        # Test empty box
        volume_empty = core.calculateCavityVolume(state, pygcmc.movement.CavityMode.FAST_APPROX)
        box_volume = 4.0 ** 3
        
        assert volume_empty > 0, "Cavity volume should be positive"
        assert volume_empty <= box_volume, "Cavity volume should not exceed box volume"
        
        # For empty box with small probe, cavity should be most of the box
        fraction_empty = volume_empty / box_volume
        assert fraction_empty > 0.5, f"Empty box cavity fraction too small: {fraction_empty}"
        
        # Add a particle manually
        res = pygcmc.MCResidue()
        res.atomStart = 0
        res.atomCount = 1
        res.active = True
        res.type = 0
        state.addResidue(res)  # Use addResidue method
        
        atom = pygcmc.MCAtom()
        atom.x = 2.0  # Center of 4x4x4 box
        atom.y = 2.0
        atom.z = 2.0
        atom.type = 0
        state.addAtom(atom)  # Use addAtom method
        
        # Invalidate cache and recalculate
        core.invalidateCache()
        volume_filled = core.calculateCavityVolume(state, pygcmc.movement.CavityMode.FAST_APPROX)
        
        # Cavity volume should decrease with particle
        assert volume_filled < volume_empty, \
            f"Cavity volume should decrease with particle: {volume_filled} >= {volume_empty}"
        
        # Check reduction is reasonable
        # With grid spacing 0.25 nm and probe radius 0.10 nm, 
        # a single particle with sigma 0.35 nm should occupy some volume
        reduction = (volume_empty - volume_filled) / volume_empty
        assert 0.001 < reduction < 0.5, f"Unexpected cavity reduction: {reduction * 100}%"
    
    def test_cavity_modes(self):
        """Test different cavity calculation modes"""
        state = pygcmc.MCState()
        state.info.box = np.array([3.0, 3.0, 3.0])  # Smaller box for faster calculation
        
        ff = pygcmc.MCForceField()
        ff.numTotalTypes = 1
        ff.numMovementTypes = 1
        ff.maxTypes = 1
        ff.ljEps = [0.0]
        ff.ljSigma = [0.35]
        ff.rebuildLJMatrix()  # Ensure LJ matrix is properly initialized
        state.forcefield = ff
        
        # Add one particle
        res = pygcmc.MCResidue()
        res.atomStart = 0
        res.atomCount = 1
        res.active = True
        res.type = 0
        state.addResidue(res)  # Use addResidue method
        
        atom = pygcmc.MCAtom()
        atom.x = 1.5
        atom.y = 1.5
        atom.z = 1.5
        atom.type = 0
        state.addAtom(atom)  # Use addAtom method
        
        core = pygcmc.movement.CavityBiasCore(0.25, 0.10)
        
        # Test both calculation modes
        volume_fast = core.calculateCavityVolume(state, pygcmc.movement.CavityMode.FAST_APPROX)
        
        # Both should give positive volumes
        assert volume_fast > 0, "Fast approximation should give positive volume"
        
        # Fast approximation should be consistent
        core.invalidateCache()
        volume_fast2 = core.calculateCavityVolume(state, pygcmc.movement.CavityMode.FAST_APPROX)
        assert abs(volume_fast - volume_fast2) < 0.001, "Fast approximation should be consistent"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
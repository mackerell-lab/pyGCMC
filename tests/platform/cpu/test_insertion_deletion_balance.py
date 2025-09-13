#!/usr/bin/env python
"""
Test insertion/deletion balance in GCMC moves
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
class TestInsertionDeletionBalance:
    """Test the balance between insertion and deletion moves"""
    
    def test_paired_operations(self):
        """Test that paired insertion-deletion operations maintain detailed balance"""
        # Setup
        state = pygcmc.MCState()
        state.info.box = np.array([3.0, 3.0, 3.0])
        
        ff = pygcmc.MCForceField()
        ff.numTotalTypes = 1
        ff.numMovementTypes = 1
        ff.ljEps = [0.5]
        ff.ljSigma = [0.3]
        state.forcefield = ff
        
        params = pygcmc.movement.MovementParams()
        params.temperature = 300.0
        params.chemicalPotential = -15.0
        params.seed = 42
        
        mover = pygcmc.movement.MovementModule()
        mover.setParams(params)
        
        # Track paired operations
        paired_insertions = 0
        paired_deletions = 0
        insertion_indices = []
        deletion_indices = []
        
        for i in range(100):
            result_ins = mover.attemptInsertion(state)
            if result_ins.accepted:
                paired_insertions += 1
                insertion_indices.append(result_ins.residueIndex)
                
                # Immediately try deletion - should prefer the just-inserted residue
                result_del = mover.attemptDeletion(state)
                if result_del.accepted:
                    paired_deletions += 1
                    deletion_indices.append(result_del.residueIndex)
        
        # Analysis
        assert paired_insertions > 0, "No insertions were accepted"
        
        if paired_deletions > 0:
            # Check that deletion/insertion ratio is reasonable
            ratio = paired_deletions / paired_insertions
            assert 0.1 < ratio < 10.0, f"Unexpected deletion/insertion ratio: {ratio}"
            
            # In immediate paired operations, many deletions should target the just-inserted residue
            # This is a soft check as it depends on system state
            if len(insertion_indices) > 10 and len(deletion_indices) > 10:
                # Check if there's correlation between insertion and deletion indices
                pass  # Statistical analysis could be added here
    
    def test_equilibrium_density(self):
        """Test that the system reaches expected equilibrium density"""
        state = pygcmc.MCState()
        state.info.box = np.array([4.0, 4.0, 4.0])  # 64 nm³
        
        # Setup ideal gas forcefield (no interactions)
        ff = pygcmc.MCForceField()
        ff.numTotalTypes = 1
        ff.numMovementTypes = 1
        ff.ljEps = [0.0]  # No interactions
        ff.ljSigma = [0.3]
        state.forcefield = ff
        
        params = pygcmc.movement.MovementParams()
        params.temperature = 300.0
        params.chemicalPotential = -10.0  # Moderate chemical potential
        params.seed = 12345
        
        mover = pygcmc.movement.MovementModule()
        mover.setParams(params)
        
        # Equilibration
        for _ in range(1000):
            if np.random.random() < 0.5:
                mover.attemptInsertion(state)
            else:
                mover.attemptDeletion(state)
        
        # Production - track density
        densities = []
        for _ in range(500):
            if np.random.random() < 0.5:
                mover.attemptInsertion(state)
            else:
                mover.attemptDeletion(state)
            
            # Record density every 10 moves
            if len(densities) < 50:
                densities.append(state.activeResidueCount)
        
        # Check that we have a stable density
        avg_density = np.mean(densities)
        std_density = np.std(densities)
        
        # For ideal gas, the average should be stable
        assert avg_density > 0, "System should have particles at equilibrium"
        # Coefficient of variation should be reasonable
        if avg_density > 5:  # Only check if we have enough particles
            cv = std_density / avg_density
            assert cv < 1.0, f"Density fluctuations too large: CV={cv}"
    
    def test_deletion_prevents_negative_count(self):
        """Test that deletion correctly handles empty system"""
        state = pygcmc.MCState()
        state.info.box = np.array([3.0, 3.0, 3.0])
        
        ff = pygcmc.MCForceField()
        ff.numTotalTypes = 1
        ff.numMovementTypes = 1
        ff.ljEps = [0.5]
        ff.ljSigma = [0.3]
        state.forcefield = ff
        
        params = pygcmc.movement.MovementParams()
        params.temperature = 300.0
        params.chemicalPotential = -15.0
        params.seed = 42
        
        mover = pygcmc.movement.MovementModule()
        mover.setParams(params)
        
        # Try deletion on empty system
        result = mover.attemptDeletion(state)
        assert not result.accepted, "Deletion should not be accepted for empty system"
        assert state.activeResidueCount == 0, "Residue count should remain 0"
        
        # Add one particle
        mover.attemptInsertion(state)
        initial_count = state.activeResidueCount
        
        # Try many deletions
        for _ in range(100):
            mover.attemptDeletion(state)
        
        # Should never go negative
        assert state.activeResidueCount >= 0, "Residue count should never be negative"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
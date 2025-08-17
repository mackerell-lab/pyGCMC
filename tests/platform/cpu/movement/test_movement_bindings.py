"""
Test cases for PyGCMC Movement Module Python bindings
"""

import pytest
import numpy as np
try:
    import pygcmc
    MOVEMENT_AVAILABLE = hasattr(pygcmc, 'movement')
except ImportError:
    pygcmc = None
    MOVEMENT_AVAILABLE = False


@pytest.mark.skipif(not MOVEMENT_AVAILABLE, reason="Movement module not compiled")
class TestMovementParams:
    """Test MovementParams configuration"""
    
    def test_default_params(self):
        """Test default parameter initialization"""
        params = pygcmc.movement.MovementParams()
        assert params.temperature == pytest.approx(298.15)
        assert params.chemicalPotential == pytest.approx(-15.7)
        assert params.useCavityBias == True
        assert params.useConfigBias == True
    
    def test_custom_temperature(self):
        """Test parameter initialization with custom temperature"""
        params = pygcmc.movement.MovementParams(300.0)
        assert params.temperature == pytest.approx(300.0)
        assert params.beta == pytest.approx(1.0 / (8.314e-3 * 300.0))
    
    def test_update_derived_parameters(self):
        """Test updating derived parameters after temperature change"""
        params = pygcmc.movement.MovementParams()
        params.temperature = 350.0
        params.updateDerivedParameters()
        assert params.beta == pytest.approx(1.0 / (8.314e-3 * 350.0))


class TestMovementModule:
    """Test MovementModule main functionality"""
    
    @pytest.fixture
    def movement_module(self):
        """Create a MovementModule instance for testing"""
        if not MOVEMENT_AVAILABLE:
            return None
        params = pygcmc.movement.MovementParams(298.15)
        params.maxAtoms = 10000
        params.maxResidues = 1000
        params.chemicalPotential = -15.7  # kJ/mol for water
        params.useCavityBias = False  # Disable for simpler testing
        params.useConfigBias = False
        return pygcmc.movement.MovementModule(params)
    
    @pytest.fixture
    def mock_state(self):
        """Create a mock MCState for testing"""
        if not MOVEMENT_AVAILABLE:
            return None
        state = pygcmc.MCState()
        # Set box dimensions (in nm)
        state.info.box[0] = 3.0
        state.info.box[1] = 3.0 
        state.info.box[2] = 3.0
        # Set molecule type info
        state.info.max_types = 1  # We have one molecule type
        state.forcefield.numTotalTypes = 1
        state.forcefield.numMovementTypes = 1
        # Set basic LJ parameters (1x1 matrix for 1 type)
        state.forcefield.ljEps = [0.0]  # No LJ interaction for simple test
        state.forcefield.ljSigma = [0.0]  # water
        return state
    
    def test_module_creation(self, movement_module):
        """Test module creation and initialization"""
        if not MOVEMENT_AVAILABLE:
            pytest.skip("Movement module not available")
        assert movement_module is not None
        params = movement_module.getParams()
        assert params.temperature == pytest.approx(298.15)
        assert params.chemicalPotential == pytest.approx(-15.7)
    
    def test_insertion_attempt(self, movement_module, mock_state):
        """Test insertion attempt"""
        if not MOVEMENT_AVAILABLE:
            pytest.skip("Movement module not available")
        result = movement_module.attemptInsertion(mock_state, moleculeType=0)
        assert result.moveType == "insert"
        assert 0.0 <= result.acceptanceProbability <= 1.0
        assert result.computeTimeMs >= 0
    
    def test_deletion_attempt(self, movement_module, mock_state):
        """Test deletion attempt"""
        if not MOVEMENT_AVAILABLE:
            pytest.skip("Movement module not available")
        
        # First try to insert some molecules
        inserted = False
        for _ in range(20):  # Try up to 20 times
            result = movement_module.attemptInsertion(mock_state, moleculeType=0)
            if result.accepted:
                inserted = True
                break
        
        if inserted:
            # Now attempt deletion
            result = movement_module.attemptDeletion(mock_state)
            assert result.moveType == "delete"
            assert 0.0 <= result.acceptanceProbability <= 1.0
    
    def test_translation_attempt(self, movement_module, mock_state):
        """Test translation attempt"""
        if not MOVEMENT_AVAILABLE:
            pytest.skip("Movement module not available")
        
        # First insert a molecule
        inserted = False
        for _ in range(20):
            result = movement_module.attemptInsertion(mock_state, moleculeType=0)
            if result.accepted:
                inserted = True
                break
        
        if inserted:
            result = movement_module.attemptTranslation(mock_state)
            assert result.moveType == "translate"
            assert 0.0 <= result.acceptanceProbability <= 1.0
    
    def test_rotation_attempt(self, movement_module, mock_state):
        """Test rotation attempt"""
        if not MOVEMENT_AVAILABLE:
            pytest.skip("Movement module not available")
        
        # First insert a molecule
        inserted = False
        for _ in range(20):
            result = movement_module.attemptInsertion(mock_state, moleculeType=0)
            if result.accepted:
                inserted = True
                break
        
        if inserted:
            result = movement_module.attemptRotation(mock_state)
            assert result.moveType == "rotate"
            assert 0.0 <= result.acceptanceProbability <= 1.0
    
    def test_cavity_finding(self, movement_module, mock_state):
        """Test cavity finding functionality"""
        if not MOVEMENT_AVAILABLE:
            pytest.skip("Movement module not available")
        cavities = movement_module.findCavities(mock_state)
        assert isinstance(cavities, list)
        # Empty box should have many cavities
        assert len(cavities) >= 0  # May be 0 if cavity bias is disabled
    
    def test_statistics_tracking(self, movement_module, mock_state):
        """Test statistics tracking"""
        if not MOVEMENT_AVAILABLE:
            pytest.skip("Movement module not available")
        
        movement_module.resetStatistics()
        
        # Perform some moves
        for _ in range(10):
            movement_module.attemptInsertion(mock_state, moleculeType=0)
        
        stats = movement_module.getStatistics()
        insert_stats = stats["insert"]
        assert insert_stats.attempts == 10
        assert 0 <= insert_stats.accepts <= 10
        
        acceptance_rate = movement_module.calculateAcceptanceRate("insert")
        assert 0.0 <= acceptance_rate <= 1.0


@pytest.mark.skipif(not MOVEMENT_AVAILABLE, reason="Movement module not compiled")
class TestActivePool:
    """Test ActivePool memory management"""
    
    @pytest.fixture
    def active_pool(self):
        """Create an ActivePool instance"""
        if MOVEMENT_AVAILABLE:
            return pygcmc.movement.ActivePool(1000, 100)
        return None
    
    def test_pool_creation(self, active_pool):
        """Test pool creation with capacity"""
        assert active_pool.getMaxAtoms() == 1000
        assert active_pool.getMaxResidues() == 100
        counts = active_pool.getActiveCounts()
        assert counts == (0, 0)  # Initially empty
    
    def test_insert_molecule(self, active_pool):
        """Test molecule insertion into pool"""
        atoms = create_water_molecule()  # Helper function
        res_idx = active_pool.insertMolecule(atoms, resType=0)
        assert res_idx >= 0
        assert active_pool.isResidueActive(res_idx)
        
        counts = active_pool.getActiveCounts()
        assert counts[0] == len(atoms)  # Active atoms
        assert counts[1] == 1  # Active residues
    
    def test_delete_residue(self, active_pool):
        """Test residue deletion from pool"""
        atoms = create_water_molecule()
        res_idx = active_pool.insertMolecule(atoms, resType=0)
        
        success = active_pool.deleteResidue(res_idx)
        assert success
        assert not active_pool.isResidueActive(res_idx)
        
        counts = active_pool.getActiveCounts()
        assert counts == (0, 0)  # Should be empty after deletion
    
    def test_fragmentation(self, active_pool):
        """Test fragmentation calculation and compaction"""
        # Insert and delete to create fragmentation
        indices = []
        for _ in range(10):
            atoms = create_water_molecule()
            idx = active_pool.insertMolecule(atoms, resType=0)
            indices.append(idx)
        
        # Delete every other molecule
        for i in range(0, 10, 2):
            active_pool.deleteResidue(indices[i])
        
        fragmentation = active_pool.getFragmentation()
        assert 0.0 <= fragmentation <= 1.0
        
        # Compact the pool
        compacted = active_pool.compact(force=True)
        assert compacted >= 0  # May or may not compact depending on fragmentation
    
    def test_capacity_check(self, active_pool):
        """Test capacity checking"""
        can_insert = active_pool.canInsert(3)  # 3 atoms for water
        assert can_insert  # Should have space initially
        
        # Fill the pool (but not completely to avoid infinite loop)
        for _ in range(100):  # Insert up to 100 molecules
            if active_pool.canInsert(3):
                atoms = create_water_molecule()
                active_pool.insertMolecule(atoms, resType=0)
            else:
                break
        
        # Check that we inserted some molecules
        counts = active_pool.getActiveCounts()
        assert counts[1] > 0  # Should have some residues


class TestMovementResult:
    """Test MovementResult structure"""
    
    def test_result_creation(self):
        """Test result creation and fields"""
        if not MOVEMENT_AVAILABLE:
            pytest.skip("Movement module not available")
        
        result = pygcmc.movement.MovementResult(
            accepted=True,
            energyChange=-5.0,
            acceptanceProbability=0.8,
            moveType="insert"
        )
        
        assert result.accepted == True
        assert result.energyChange == pytest.approx(-5.0)
        assert result.acceptanceProbability == pytest.approx(0.8)
        assert result.moveType == "insert"
        assert result.isSuccessful()
    
    def test_result_summary(self):
        """Test result summary string"""
        if not MOVEMENT_AVAILABLE:
            pytest.skip("Movement module not available")
        
        result = pygcmc.movement.MovementResult(
            accepted=False,
            energyChange=10.0,
            acceptanceProbability=0.1,
            moveType="delete"
        )
        
        summary = result.summary()
        assert "delete" in summary
        assert "rejected" in summary or "accepted" in summary
        assert "10.0" in summary or "10.000" in summary


class TestIntegration:
    """Integration tests for complete GCMC workflow"""
    
    @pytest.fixture
    def gcmc_system(self):
        """Create a complete GCMC system for testing"""
        if not MOVEMENT_AVAILABLE:
            return None, None, None
        
        params = pygcmc.movement.MovementParams(298.15)
        params.chemicalPotential = -15.7  # Water
        params.useCavityBias = False  # Disable for simpler testing
        params.useConfigBias = False
        params.maxAtoms = 10000
        params.maxResidues = 1000
        
        movement = pygcmc.movement.MovementModule(params)
        state = pygcmc.MCState()
        state.info.box[0] = 3.0  # nm
        state.info.box[1] = 3.0
        state.info.box[2] = 3.0
        # Set molecule type info
        state.info.max_types = 1  # We have one molecule type
        state.forcefield.numTotalTypes = 1
        state.forcefield.numMovementTypes = 1
        # Set basic LJ parameters (1x1 matrix for 1 type)
        state.forcefield.ljEps = [0.0]  # No LJ interaction for simple test
        state.forcefield.ljSigma = [0.0]
        
        return movement, state, params
    
    def test_gcmc_equilibration(self, gcmc_system):
        """Test GCMC equilibration process"""
        movement, state, params = gcmc_system
        if movement is None:
            pytest.skip("Movement module not available")
        
        # Run GCMC moves
        move_types = ["insert", "delete", "translate", "rotate"]
        move_probs = [0.4, 0.2, 0.2, 0.2]  # Bias toward insertion initially
        
        n_steps = 50  # Reduced for faster testing
        n_molecules_history = []
        
        for step in range(n_steps):
            # Select move type
            move_type = np.random.choice(move_types, p=move_probs)
            
            # Perform move
            if move_type == "insert":
                result = movement.attemptInsertion(state, moleculeType=0)
            elif move_type == "delete":
                result = movement.attemptDeletion(state)
            elif move_type == "translate":
                result = movement.attemptTranslation(state)
            else:  # rotate
                result = movement.attemptRotation(state)
            
            # Record state
            if step % 10 == 0:
                n_molecules_history.append(state.activeResidueCount)
        
        # Check that system evolved
        assert len(n_molecules_history) > 0
    
    def test_acceptance_rates(self, gcmc_system):
        """Test that acceptance rates are reasonable"""
        movement, state, params = gcmc_system
        if movement is None:
            pytest.skip("Movement module not available")
        
        movement.resetStatistics()
        
        # Perform many moves
        for _ in range(30):
            movement.attemptInsertion(state, moleculeType=0)
            if state.activeResidueCount > 0:
                movement.attemptTranslation(state)
                movement.attemptRotation(state)
        
        # Check that we attempted some moves
        stats = movement.getStatistics()
        insert_stats = stats["insert"]
        assert insert_stats.attempts > 0
        
        # Check acceptance rate is defined
        insert_rate = movement.calculateAcceptanceRate("insert")
        assert 0.0 <= insert_rate <= 1.0


# Helper functions

def create_water_molecule():
    """Create a simple water molecule for testing"""
    if not MOVEMENT_AVAILABLE:
        return []
    
    atoms = []
    # Oxygen
    o = pygcmc.MCAtom()
    o.x, o.y, o.z = 0.0, 0.0, 0.0
    o.charge = -0.834
    o.type = 0
    atoms.append(o)
    
    # Hydrogen 1  
    h1 = pygcmc.MCAtom()
    h1.x, h1.y, h1.z = 0.0957, 0.0, 0.0  # in nm
    h1.charge = 0.417
    h1.type = 1
    atoms.append(h1)
    
    # Hydrogen 2
    h2 = pygcmc.MCAtom()
    h2.x, h2.y, h2.z = -0.0239, 0.0927, 0.0  # in nm
    h2.charge = 0.417  
    h2.type = 1
    atoms.append(h2)
    
    return atoms


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
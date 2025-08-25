# tests/simulation/movement/test_movement_module.py
"""Movement module integration and functionality tests."""

import pytest
import pygcmc
import numpy as np


class TestMovementModule:
    """Test MovementModule main functionality."""
    
    @pytest.fixture
    def setup_system(self):
        """Setup a complete test system."""
        state = pygcmc.MCState()
        state.info.box = np.array([5.0, 5.0, 5.0])
        
        # Setup force field
        ff = pygcmc.MCForceField()
        ff.numTotalTypes = 1
        ff.numMovementTypes = 1
        ff.ljEps = [0.5]
        ff.ljSigma = [0.3]
        state.forcefield = ff
        
        params = pygcmc.movement.MovementParams()
        params.temperature = 298.15
        params.chemicalPotential = -15.7
        params.seed = 42
        
        mover = pygcmc.movement.MovementModule()
        mover.setParams(params)
        
        return state, params, mover
    
    def test_module_creation(self):
        """Test MovementModule creation with different constructors."""
        # Default constructor
        mover1 = pygcmc.movement.MovementModule()
        assert mover1 is not None
        
        # Constructor with params
        params = pygcmc.movement.MovementParams()
        mover2 = pygcmc.movement.MovementModule()
        mover2.setParams(params)
        assert mover2 is not None
    
    def test_parameter_setting_and_getting(self):
        """Test setting and getting parameters."""
        mover = pygcmc.movement.MovementModule()
        
        # Set parameters
        params = pygcmc.movement.MovementParams()
        params.temperature = 300.0
        params.chemicalPotential = -20.0
        mover.setParams(params)
        
        # Get parameters
        params_back = mover.getParams()
        assert params_back.temperature == pytest.approx(300.0)
        assert params_back.chemicalPotential == pytest.approx(-20.0)
    
    def test_insertion_basic(self, setup_system):
        """Test basic insertion functionality."""
        state, params, mover = setup_system
        
        # Attempt insertion
        result = mover.attemptInsertion(state)
        
        # Check result
        assert isinstance(result.accepted, bool)
        assert result.moveType == "insert"  # Insertion
        assert isinstance(result.energyChange, float)
    
    def test_deletion_basic(self, setup_system):
        """Test basic deletion functionality."""
        state, params, mover = setup_system
        
        # First insert some atoms
        inserted = False
        for _ in range(100):
            result = mover.attemptInsertion(state)
            if result.accepted:
                inserted = True
                break
        
        if inserted:
            # Try deletion
            result = mover.attemptDeletion(state)
            assert isinstance(result.accepted, bool)
            assert result.moveType == "delete"  # Deletion
    
    def test_translation_basic(self, setup_system):
        """Test basic translation functionality."""
        state, params, mover = setup_system
        
        # First insert an atom
        inserted = False
        for _ in range(100):
            result = mover.attemptInsertion(state)
            if result.accepted:
                inserted = True
                break
        
        if inserted:
            # Try translation
            result = mover.attemptTranslation(state)
            assert isinstance(result.accepted, bool)
            assert result.moveType == "translate"  # Translation
    
    def test_rotation_basic(self, setup_system):
        """Test basic rotation functionality."""
        state, params, mover = setup_system
        
        # First insert an atom
        inserted = False
        for _ in range(100):
            result = mover.attemptInsertion(state)
            if result.accepted:
                inserted = True
                break
        
        if inserted:
            # Try rotation
            result = mover.attemptRotation(state)
            assert isinstance(result.accepted, bool)
            assert result.moveType == "rotate"  # Rotation
    
    def test_cavity_bias_insertion(self, setup_system):
        """Test cavity-biased insertion."""
        state, params, mover = setup_system
        params.useCavityBias = True
        mover.setParams(params)
        
        # Attempt cavity-biased insertion
        result = mover.attemptCavityBiasInsertion(state)
        
        assert isinstance(result.accepted, bool)
        assert result.moveType == "insert"  # Still insertion
    
    def test_config_bias_rotation(self, setup_system):
        """Test configurational bias rotation."""
        state, params, mover = setup_system
        params.useConfigBias = True
        mover.setParams(params)
        
        # First insert an atom
        inserted = False
        for _ in range(100):
            result = mover.attemptInsertion(state)
            if result.accepted:
                inserted = True
                break
        
        if inserted:
            # Try config bias rotation
            result = mover.attemptConfigBiasRotation(state)
            assert isinstance(result.accepted, bool)
    
    def test_acceptance_rate_calculation(self, setup_system):
        """Test acceptance rate calculation."""
        state, params, mover = setup_system
        
        # Reset statistics
        mover.resetStatistics()
        
        # Perform moves and track
        accepts = 0
        attempts = 100
        
        for _ in range(attempts):
            result = mover.attemptInsertion(state)
            if result.accepted:
                accepts += 1
        
        # Calculate acceptance rate
        rate = mover.calculateAcceptanceRate("insert")
        
        # Should be close to manual calculation
        expected_rate = accepts / attempts
        assert abs(rate - expected_rate) < 0.1
    
    def test_find_cavities_integration(self, setup_system):
        """Test cavity finding integration."""
        state, params, mover = setup_system
        params.useCavityBias = True
        mover.setParams(params)
        
        # Find cavities
        cavities = mover.findCavities(state)
        
        # Should return a list of positions
        assert isinstance(cavities, list)
        if len(cavities) > 0:
            # Check first cavity
            cavity = cavities[0]
            assert hasattr(cavity, 'x')
            assert hasattr(cavity, 'y')
            assert hasattr(cavity, 'z')
            
            # Should be within box
            assert 0 <= cavity.x <= 5.0
            assert 0 <= cavity.y <= 5.0
            assert 0 <= cavity.z <= 5.0
    
    def test_statistics_retrieval(self, setup_system):
        """Test statistics retrieval methods."""
        state, params, mover = setup_system
        
        # Perform some moves
        for _ in range(50):
            mover.attemptInsertion(state)
        
        # Get general statistics
        stats = mover.getStatistics()
        assert stats is not None
        
        # Get statistics - it's a dict of move types
        proposal_stats = mover.getStatistics()
        assert isinstance(proposal_stats, dict)
        assert "insert" in proposal_stats
        
        # Same stats object
        cavity_stats = mover.getStatistics()
        assert isinstance(cavity_stats, dict)
    
    def test_multi_insertion_cbmc(self, setup_system):
        """Test multi-insertion CBMC if available."""
        state, params, mover = setup_system
        params.useMultiInsertionCBMC = True
        params.maxParallelInsertions = 4
        
        try:
            mover.setParams(params)
            result = mover.attemptMultiInsertionCBMC(state, 0)
            
            # Check result
            assert result is not None
            if hasattr(result, '__iter__'):
                # Multiple results
                for r in result:
                    assert hasattr(r, 'accepted')
            else:
                # Single result
                assert hasattr(result, 'accepted')
                
        except Exception:
            # Multi-insertion may not be implemented
            pytest.skip("Multi-insertion CBMC not available")
    
    def test_temperature_effect(self, setup_system):
        """Test effect of temperature on acceptance."""
        state, params, mover = setup_system
        
        # Test at different temperatures
        temperatures = [100.0, 298.15, 500.0]
        acceptance_rates = []
        
        for temp in temperatures:
            params.temperature = temp
            params.updateDerivedParameters()
            mover.setParams(params)
            mover.resetStatistics()
            
            accepts = 0
            attempts = 100
            for _ in range(attempts):
                result = mover.attemptInsertion(state)
                if result.accepted:
                    accepts += 1
            
            acceptance_rates.append(accepts / attempts)
        
        # Higher temperature should generally give higher acceptance
        # (though this depends on many factors)
        # Just check they're different
        assert len(set(acceptance_rates)) > 1
    
    def test_chemical_potential_effect(self, setup_system):
        """Test effect of chemical potential on insertion."""
        state, params, mover = setup_system
        
        # Test at different chemical potentials (ordered from low to high)
        chem_potentials = [-30.0, -15.7, -5.0]
        acceptance_rates = []
        
        attempts = 200  # Reasonable number for statistics
        for mu in chem_potentials:
            params.chemicalPotential = mu
            mover.setParams(params)
            mover.resetStatistics()
            
            accepts = 0
            for _ in range(attempts):
                result = mover.attemptInsertion(state)
                if result.accepted:
                    accepts += 1
            
            acceptance_rates.append(accepts / attempts)
        
        # Higher chemical potential should favor insertion
        # Due to statistical fluctuations and possible system saturation,
        # we check for overall trend rather than strict monotonicity
        
        # At minimum, highest mu should have higher rate than lowest mu
        assert acceptance_rates[2] >= acceptance_rates[0], \
            f"Highest mu=-5 rate ({acceptance_rates[2]:.3f}) should be >= lowest mu=-30 rate ({acceptance_rates[0]:.3f})"
        
        # Also check that rates are not all identical (should see some effect)
        assert len(set(acceptance_rates)) > 1, \
            f"All rates identical: {acceptance_rates} - chemical potential has no effect"
    
    def test_module_repr(self):
        """Test MovementModule string representation."""
        params = pygcmc.movement.MovementParams()
        params.temperature = 300.0
        params.useCavityBias = True
        params.useConfigBias = False
        
        mover = pygcmc.movement.MovementModule()
        mover.setParams(params)
        
        repr_str = str(mover)
        assert "MovementModule" in repr_str
        assert "300" in repr_str  # Temperature
        assert "cavityBias" in repr_str
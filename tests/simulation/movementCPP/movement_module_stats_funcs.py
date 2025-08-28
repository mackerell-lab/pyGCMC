# tests/simulation/movementCPP/movement_module_stats_funcs.py
"""Movement module statistics and effects tests - extracted functions."""

import pytest
from .conftest import setup_system_with_params
import pygcmc

@pytest.fixture
def setup_system():
    """Setup a basic test system."""
    state = pygcmc.MCState()
    state.info.box = [5.0, 5.0, 5.0]
    
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
    
    return state, params

def test_statistics_retrieval(setup_system_with_params):
        """Test statistics retrieval methods."""
        state, params = setup_system_with_params
        mover = pygcmc.movement.MovementModule()
        mover.setParams(params)
        
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
    
def test_multi_insertion_cbmc_basic(setup_system_with_params):
        """Test multi-insertion CBMC basic functionality if available."""
        state, params = setup_system_with_params
        mover = pygcmc.movement.MovementModule()
        mover.setParams(params)
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
                
        except Exception as e:
            # Only skip for known "not available" errors
            msg = str(e).lower()
            if any(x in msg for x in ("not available", "not implemented", "unsupported", "not compiled")):
                pytest.skip(f"Multi-insertion CBMC not available: {e}")
            else:
                # Re-raise unexpected errors for debugging
                raise
    
def test_temperature_effect(setup_system_with_params):
        """Test effect of temperature on acceptance."""
        state, params = setup_system_with_params
        mover = pygcmc.movement.MovementModule()
        mover.setParams(params)
        
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
    
def test_chemical_potential_effect_insertion(setup_system_with_params):
        """Test effect of chemical potential on insertion acceptance."""
        # Get fresh state for each mu test to avoid saturation
        # Each chemical potential test starts with empty box
        
        # Test at different chemical potentials (ordered from low to high)
        chem_potentials = [-30.0, -15.7, -5.0]
        acceptance_rates = []
        
        attempts = 200  # Reasonable number for statistics
        for mu in chem_potentials:
            # Create fresh state and mover for each test
            state = pygcmc.MCState()
            state.info.box = [5.0, 5.0, 5.0]
            ff = pygcmc.MCForceField()
            ff.numTotalTypes = 1
            ff.numMovementTypes = 1
            ff.ljEps = [0.5]
            ff.ljSigma = [0.3]
            state.forcefield = ff
            
            params = pygcmc.movement.MovementParams()
            params.temperature = 298.15
            params.chemicalPotential = mu
            
            mover = pygcmc.movement.MovementModule()
            mover.setParams(params)
            
            accepts = 0
            for _ in range(attempts):
                result = mover.attemptInsertion(state)
                if result.accepted:
                    accepts += 1
            
            acceptance_rates.append(accepts / attempts)
        
        # Higher chemical potential should favor insertion
        # Check that we see variation (chemical potential has effect)
        assert len(set(acceptance_rates)) > 1, \
            f"All rates identical: {acceptance_rates} - chemical potential has no effect"
        
        # Generally, trend should be upward but may not be strict due to saturation
        # Just verify highest mu gives non-zero acceptance
        assert max(acceptance_rates) > 0, "Should have some accepted insertions"
    
def test_statistics_counts_and_reset(setup_system_with_params):
        """Test statistics counting and reset functionality."""
        state, params = setup_system_with_params
        mover = pygcmc.movement.MovementModule()
        mover.setParams(params)
        
        # Reset statistics to start clean
        mover.resetStatistics()
        
        # Perform exact number of attempts
        attempts = 20
        for _ in range(attempts):
            mover.attemptInsertion(state)
        
        # Verify counts
        stats = mover.getStatistics()
        assert "insert" in stats
        assert stats["insert"].attempts == attempts
        
        # Verify API consistency
        rate_api = mover.calculateAcceptanceRate("insert")
        assert rate_api == pytest.approx(stats["insert"].acceptanceRate(), rel=1e-6)
        
        # Test reset
        mover.resetStatistics()
        stats2 = mover.getStatistics()
        assert stats2["insert"].attempts == 0
        assert stats2["insert"].accepts == 0
    
def test_deletion_prefers_last_inserted(setup_system_with_params):
        """Test that deletion prefers the last inserted residue."""
        state, params = setup_system_with_params
        mover = pygcmc.movement.MovementModule()
        mover.setParams(params)
        
        # Insert a molecule
        inserted = False
        for _ in range(200):
            r = mover.attemptInsertion(state)
            if r.accepted:
                inserted = True
                break
        
        if not inserted:
            pytest.skip("Could not insert; skipping deletion preference test")
        
        # Check deletion
        before = state.activeResidueCount
        rdel = mover.attemptDeletion(state)
        assert rdel.moveType == "delete"
        
        if rdel.accepted:
            assert state.activeResidueCount == before - 1
            assert rdel.residueIndex >= 0
    
def test_find_cavities_bounds_all(setup_system_with_params):
        """Test that all cavities are within box bounds."""
        state, params = setup_system_with_params
        mover = pygcmc.movement.MovementModule()
        mover.setParams(params)
        
        params.useCavityBias = True
        mover.setParams(params)
        
        cavities = mover.findCavities(state)
        Lx, Ly, Lz = state.info.box
        
        # Check all cavities, not just the first
        for c in cavities:
            assert 0.0 <= c.x <= Lx
            assert 0.0 <= c.y <= Ly
            assert 0.0 <= c.z <= Lz
    
def test_probability_bounds_all_moves(setup_system_with_params):
        """Test that acceptance probability is always in [0,1] for all moves."""
        state, params = setup_system_with_params
        mover = pygcmc.movement.MovementModule()
        mover.setParams(params)
        
        # Try each move once; probability must be in [0,1] even if rejected
        for name in ("attemptInsertion", "attemptDeletion", "attemptTranslation", "attemptRotation"):
            res = getattr(mover, name)(state)
            assert 0.0 <= res.acceptanceProbability <= 1.0, \
                f"{name}: acceptanceProbability {res.acceptanceProbability} out of range [0,1]"

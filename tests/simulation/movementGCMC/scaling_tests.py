# tests/simulation/movementGCMC/scaling_tests.py
"""
GCMC scaling and ensemble tests
"""
import pytest
import pygcmc
import numpy as np
from .conftest import run_gcmc_steps


def test_volume_scaling(small_state):
    """Test that molecule count scales with volume"""
    # Create states with different volumes
    state_small = small_state  # 2x2x2 = 8 nm³
    
    state_large = pygcmc.MCState()
    state_large.info = pygcmc.MCInfo()
    state_large.info.box = [4.0, 4.0, 4.0]  # 4x4x4 = 64 nm³
    state_large.info.setTemperature(300.0)
    state_large.info.cutoff = 1.2
    state_large.info.volume = 64.0
    
    # Setup force field for large box
    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 10
    ff.maxTypes = 10
    ff.ljSigma = [0.0] * 100
    ff.ljEps = [0.0] * 100
    ff.ljSigma[0] = 0.3151
    ff.ljEps[0] = 0.6364
    state_large.forcefield = ff
    state_large.residues = []
    state_large.atoms = []
    
    # Same parameters for both
    params = pygcmc.movement.MovementParams()
    params.temperature = 300.0
    params.chemicalPotential = -15.7
    params.maxTranslation = 0.1
    params.maxRotation = 0.2
    
    # Run simulations
    mover_small = pygcmc.movement.MovementModule()
    mover_small.setParams(params)
    
    mover_large = pygcmc.movement.MovementModule()
    mover_large.setParams(params)
    
    # Insert molecules
    for _ in range(200):
        mover_small.attemptInsertion(state_small)
        mover_large.attemptInsertion(state_large)
    
    # Volume ratio is 64/8 = 8
    # Molecule count should scale approximately with volume
    volume_ratio = 64.0 / 8.0
    
    if state_small.activeResidueCount > 0:
        molecule_ratio = state_large.activeResidueCount / state_small.activeResidueCount
        
        # Should be roughly proportional (with wide tolerance for MC noise)
        assert 0.2 * volume_ratio <= molecule_ratio <= 5.0 * volume_ratio, \
            f"Molecule ratio {molecule_ratio} not proportional to volume ratio {volume_ratio}"


def test_chemical_potential_scaling():
    """Test that chemical potential affects molecule count correctly"""
    # Create identical states
    state1 = pygcmc.MCState()
    state2 = pygcmc.MCState()
    
    for state in [state1, state2]:
        state.info = pygcmc.MCInfo()
        state.info.box = [3.0, 3.0, 3.0]
        state.info.setTemperature(300.0)
        state.info.cutoff = 1.2
        state.info.volume = 27.0
        
        ff = pygcmc.MCForceField()
        ff.numTotalTypes = 10
        ff.maxTypes = 10
        ff.ljSigma = [0.0] * 100
        ff.ljEps = [0.0] * 100
        ff.ljSigma[0] = 0.3151
        ff.ljEps[0] = 0.6364
        state.forcefield = ff
        state.residues = []
        state.atoms = []
    
    # Different chemical potentials
    params_low = pygcmc.movement.MovementParams()
    params_low.temperature = 300.0
    params_low.chemicalPotential = -18.0  # Low
    
    params_high = pygcmc.movement.MovementParams()
    params_high.temperature = 300.0
    params_high.chemicalPotential = -13.0  # High
    
    mover_low = pygcmc.movement.MovementModule()
    mover_low.setParams(params_low)
    
    mover_high = pygcmc.movement.MovementModule()
    mover_high.setParams(params_high)
    
    # Run insertions
    for _ in range(300):
        mover_low.attemptInsertion(state1)
        mover_high.attemptInsertion(state2)
    
    # Higher chemical potential should give more molecules
    assert state2.activeResidueCount >= state1.activeResidueCount, \
        f"Higher μ ({params_high.chemicalPotential}) gave fewer molecules than lower μ ({params_low.chemicalPotential})"
    
    # Calculate expected ratio (rough approximation)
    beta = 1.0 / (8.314e-3 * 300.0)  # 1/kT in mol/kJ
    delta_mu = params_high.chemicalPotential - params_low.chemicalPotential
    expected_ratio = np.exp(beta * delta_mu)
    
    if state1.activeResidueCount > 0:
        actual_ratio = state2.activeResidueCount / state1.activeResidueCount
        
        # Should be in right order of magnitude (large tolerance for MC)
        assert 0.1 * expected_ratio <= actual_ratio <= 10.0 * expected_ratio


def test_detailed_balance():
    """Test that detailed balance is maintained"""
    state = pygcmc.MCState()
    state.info = pygcmc.MCInfo()
    state.info.box = [3.0, 3.0, 3.0]
    state.info.setTemperature(300.0)
    state.info.cutoff = 1.2
    state.info.volume = 27.0
    
    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 10
    ff.maxTypes = 10
    ff.ljSigma = [0.0] * 100
    ff.ljEps = [0.0] * 100
    ff.ljSigma[0] = 0.3151
    ff.ljEps[0] = 0.6364
    state.forcefield = ff
    state.residues = []
    state.atoms = []
    
    params = pygcmc.movement.MovementParams()
    params.temperature = 300.0
    params.chemicalPotential = -15.7
    
    mover = pygcmc.movement.MovementModule()
    mover.setParams(params)
    mover.resetStatistics()
    
    # Run many steps to reach equilibrium
    for _ in range(1000):
        r = np.random.random()
        if r < 0.5:
            mover.attemptInsertion(state)
        elif state.activeResidueCount > 0:
            mover.attemptDeletion(state)
    
    # Get statistics
    stats = mover.getStatistics()
    
    # In equilibrium, insertion and deletion rates should balance
    # This is a weak test due to finite sampling
    if hasattr(stats, 'insert') and hasattr(stats, 'delete'):
        if stats.insert.attempts > 100 and stats.delete.attempts > 100:
            insert_rate = stats.insert.accepts / stats.insert.attempts
            delete_rate = stats.delete.accepts / stats.delete.attempts
            
            # Rates should be positive
            assert insert_rate > 0
            assert delete_rate > 0
            
            # Can't expect exact balance due to finite sampling
            # Just check they're in same order of magnitude
            ratio = insert_rate / delete_rate if delete_rate > 0 else 0
            assert 0.1 <= ratio <= 10.0, \
                f"Insert/delete balance off: {insert_rate:.3f}/{delete_rate:.3f}"
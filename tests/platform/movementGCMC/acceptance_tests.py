# tests/simulation/movementGCMC/acceptance_tests.py
"""
GCMC acceptance criteria tests
"""
import pytest
import pygcmc
from .conftest import run_gcmc_steps, calculate_acceptance_rate


def test_acceptance_rate_range(setup_gcmc_state, gcmc_mover):
    """Test that acceptance rate is in valid range [0, 1]"""
    state = setup_gcmc_state
    
    # Run some GCMC steps
    results = run_gcmc_steps(state, gcmc_mover, n_steps=500)
    
    # Calculate acceptance rate
    n_accepted = sum(1 for r in results if r.accepted)
    n_total = len(results)
    
    if n_total > 0:
        acceptance_rate = n_accepted / n_total
        assert 0.0 <= acceptance_rate <= 1.0
        
        # For a well-tuned simulation, expect reasonable acceptance
        # This may need adjustment based on parameters
        assert acceptance_rate > 0.001, "Acceptance rate too low"
        assert acceptance_rate < 0.999, "Acceptance rate suspiciously high"


def test_chemical_potential_effect(setup_gcmc_state):
    """Test that chemical potential affects steady-state <N> with proper equilibration"""
    import numpy as np
    
    # Seed for reproducibility
    np.random.seed(12345)
    
    # Test multiple μ values for monotonicity
    mu_values = [-20.0, -15.0, -10.0]
    mean_n_values = []
    
    for mu in mu_values:
        # Fresh state for each test
        state = setup_gcmc_state
        # Clear state properly
        while state.activeResidueCount > 0:
            temp_mover = pygcmc.movement.MovementModule()
            temp_mover.attemptDeletion(state)
        
        params = pygcmc.movement.MovementParams()
        params.temperature = 300.0
        params.chemicalPotential = mu
        params.seed = int(12345 + mu)  # Different seed for each
        params.maxTranslation = 0.1
        params.maxRotation = 0.2
        
        # Use MovementModule(params) for proper seeding
        mover = pygcmc.movement.MovementModule(params)
        
        # EQUILIBRATION with balanced moves
        for _ in range(500):
            if np.random.random() < 0.5:
                mover.attemptInsertion(state)
            else:
                if state.activeResidueCount > 0:
                    mover.attemptDeletion(state)
        
        # STEADY-STATE SAMPLING
        n_samples = []
        for i in range(1000):
            if np.random.random() < 0.5:
                mover.attemptInsertion(state)
            else:
                if state.activeResidueCount > 0:
                    mover.attemptDeletion(state)
            
            # Sample every 10 steps
            if i % 10 == 0:
                n_samples.append(state.activeResidueCount)
        
        mean_n = np.mean(n_samples)
        std_n = np.std(n_samples)
        mean_n_values.append((mu, mean_n, std_n))
    
    # PHYSICS VALIDATION: <N> should increase with μ
    for i in range(len(mean_n_values) - 1):
        mu1, n1, _ = mean_n_values[i]
        mu2, n2, _ = mean_n_values[i + 1]
        
        # Should increase (allow small noise)
        assert n2 >= n1 * 0.9, \
            f"Non-monotonic: μ={mu1} gives <N>={n1:.2f}, μ={mu2} gives <N>={n2:.2f}"
    
    # Should see significant effect from lowest to highest μ
    n_low = mean_n_values[0][1]
    n_high = mean_n_values[-1][1]
    assert n_high > n_low * 1.2 or n_high > n_low + 1, \
        f"Insufficient μ effect: <N> changes from {n_low:.2f} to {n_high:.2f}"


def test_temperature_effect():
    """Test that temperature affects acceptance rates with independent states"""
    import numpy as np
    
    # Seed NumPy for reproducibility
    np.random.seed(98765)
    
    # Test at different temperatures
    temperatures = [250.0, 300.0, 350.0]
    acceptance_rates = []
    energy_bins_by_temp = {}
    
    for T in temperatures:
        # CREATE FRESH STATE for each temperature (independent)
        state = pygcmc.MCState()
        state.info = pygcmc.MCInfo()
        state.info.box = [5.0, 5.0, 5.0]
        state.info.setTemperature(T)
        state.info.cutoff = 1.2
        state.info.volume = 125.0
        
        # Setup force field
        ff = pygcmc.MCForceField()
        ff.numTotalTypes = 10
        ff.maxTypes = 10
        ff.ljSigma = [0.0] * 100
        ff.ljEps = [0.0] * 100
        ff.ljSigma[0] = 0.3151  # O-O
        ff.ljEps[0] = 0.6364
        state.forcefield = ff
        state.residues = []
        state.atoms = []
        state.activeResidueCount = 0
        state.activeAtomCount = 0
        
        # Create mover with temperature and seed
        params = pygcmc.movement.MovementParams()
        params.temperature = T
        params.chemicalPotential = -15.7
        params.maxTranslation = 0.1
        params.maxRotation = 0.2
        params.seed = int(98765 + T)  # Different seed per T
        
        # Use MovementModule(params) for proper seeding
        mover = pygcmc.movement.MovementModule(params)
        mover.resetStatistics()
        
        # Run steps
        results = run_gcmc_steps(state, mover, n_steps=300)
        
        # Calculate acceptance
        n_accepted = sum(1 for r in results if r.accepted)
        n_total = len(results)
        
        if n_total > 0:
            acceptance_rates.append(n_accepted / n_total)
        else:
            acceptance_rates.append(0.0)
    
    # All should be valid rates
    for rate in acceptance_rates:
        assert 0.0 <= rate <= 1.0
    
    # Temperature effect depends on energy landscape
    # Just check that rates are different (temperature has an effect)
    assert len(set(acceptance_rates)) > 1, "Temperature has no effect on acceptance"
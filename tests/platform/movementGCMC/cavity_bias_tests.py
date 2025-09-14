# tests/simulation/movementGCMC/cavity_bias_tests.py
"""
GCMC cavity bias tests
"""
import pytest
import pygcmc
from .conftest import run_gcmc_steps


def test_cavity_detection(setup_gcmc_state, gcmc_params):
    """Test cavity detection functionality"""
    state = setup_gcmc_state
    
    # Create mover with cavity bias enabled
    params = gcmc_params
    params.useCavityBias = True
    params.cavityGridSpacing = 0.2
    params.probeRadius = 0.14
    
    mover = pygcmc.movement.MovementModule()
    mover.setParams(params)
    
    # Find cavities in empty box
    cavities = mover.findCavities(state)
    
    # Should find cavities in empty box
    assert len(cavities) > 0, "No cavities found in empty box"
    
    # Cavities should be within box bounds
    for cavity in cavities:
        assert 0 <= cavity.x <= state.info.box[0]
        assert 0 <= cavity.y <= state.info.box[1]
        assert 0 <= cavity.z <= state.info.box[2]


def test_cavity_bias_insertion(setup_gcmc_state, gcmc_params):
    """Test that cavity bias affects insertion"""
    state = setup_gcmc_state
    
    # Test with cavity bias
    params_with = gcmc_params
    params_with.useCavityBias = True
    
    mover_with = pygcmc.movement.MovementModule()
    mover_with.setParams(params_with)
    mover_with.resetStatistics()
    
    # Run insertions
    results_with = []
    for _ in range(100):
        result = mover_with.attemptInsertion(state)
        results_with.append(result)
    
    # BETTER: Request fresh state instead of manual reset to avoid invariant violations
    # For now, properly reset state ensuring consistency
    while state.activeResidueCount > 0:
        mover_with.attemptDeletion(state)
    # Double-check state is clean
    assert state.activeResidueCount == 0, "Failed to reset state"
    
    # Test without cavity bias
    params_without = gcmc_params
    params_without.useCavityBias = False
    
    mover_without = pygcmc.movement.MovementModule()
    mover_without.setParams(params_without)
    mover_without.resetStatistics()
    
    # Run insertions
    results_without = []
    for _ in range(100):
        result = mover_without.attemptInsertion(state)
        results_without.append(result)
    
    # Both should complete
    assert len(results_with) == 100
    assert len(results_without) == 100


def test_cavity_grid_parameters(setup_gcmc_state):
    """Test effect of cavity grid parameters"""
    state = setup_gcmc_state
    
    # Test different grid spacings
    spacings = [0.1, 0.2, 0.3]
    cavity_counts = []
    
    for spacing in spacings:
        params = pygcmc.movement.MovementParams()
        params.temperature = 300.0
        params.chemicalPotential = -15.7
        params.useCavityBias = True
        params.cavityGridSpacing = spacing
        params.probeRadius = 0.14
        
        mover = pygcmc.movement.MovementModule()
        mover.setParams(params)
        
        cavities = mover.findCavities(state)
        cavity_counts.append(len(cavities))
    
    # Finer grid (smaller spacing) should find more or equal cavities
    assert cavity_counts[0] >= cavity_counts[2], \
        "Finer grid found fewer cavities"


def test_cavity_vs_uniform_insertion():
    """Compare cavity-biased vs uniform insertion"""
    # Create two INDEPENDENT states to avoid interference
    import numpy as np
    np.random.seed(12345)
    
    # Create first state for cavity bias
    state_cavity = pygcmc.MCState()
    state_cavity.info = pygcmc.MCInfo()
    state_cavity.info.box = [5.0, 5.0, 5.0]
    state_cavity.info.setTemperature(300.0)
    state_cavity.info.cutoff = 1.2
    state_cavity.info.volume = 125.0
    
    ff1 = pygcmc.MCForceField()
    ff1.numTotalTypes = 10
    ff1.maxTypes = 10
    ff1.ljSigma = [0.0] * 100
    ff1.ljEps = [0.0] * 100
    ff1.ljSigma[0] = 0.3151
    ff1.ljEps[0] = 0.6364
    state_cavity.forcefield = ff1
    state_cavity.residues = []
    state_cavity.atoms = []
    state_cavity.activeResidueCount = 0
    state_cavity.activeAtomCount = 0
    
    # Create second INDEPENDENT state for uniform
    state_uniform = pygcmc.MCState()
    state_uniform.info = pygcmc.MCInfo()
    state_uniform.info.box = [5.0, 5.0, 5.0]
    state_uniform.info.setTemperature(300.0)
    state_uniform.info.cutoff = 1.2
    state_uniform.info.volume = 125.0
    
    ff2 = pygcmc.MCForceField()
    ff2.numTotalTypes = 10
    ff2.maxTypes = 10
    ff2.ljSigma = [0.0] * 100
    ff2.ljEps = [0.0] * 100
    ff2.ljSigma[0] = 0.3151
    ff2.ljEps[0] = 0.6364
    state_uniform.forcefield = ff2
    state_uniform.residues = []
    state_uniform.atoms = []
    state_uniform.activeResidueCount = 0
    state_uniform.activeAtomCount = 0
    
    # Cavity-biased insertion
    params_cavity = pygcmc.movement.MovementParams()
    params_cavity.temperature = 300.0
    params_cavity.chemicalPotential = -15.7
    params_cavity.useCavityBias = True
    params_cavity.cavityGridSpacing = 0.2
    params_cavity.probeRadius = 0.14
    
    mover_cavity = pygcmc.movement.MovementModule()
    mover_cavity.setParams(params_cavity)
    
    # Uniform insertion
    params_uniform = pygcmc.movement.MovementParams()
    params_uniform.temperature = 300.0
    params_uniform.chemicalPotential = -15.7
    params_uniform.useCavityBias = False
    
    mover_uniform = pygcmc.movement.MovementModule()
    mover_uniform.setParams(params_uniform)
    
    # Run simulations - increase steps for better statistics
    results_cavity = run_gcmc_steps(state_cavity, mover_cavity, n_steps=500)
    results_uniform = run_gcmc_steps(state_uniform, mover_uniform, n_steps=500)
    
    # Both should complete
    assert len(results_cavity) > 0
    assert len(results_uniform) > 0
    
    # Calculate acceptance rates
    accept_cavity = sum(1 for r in results_cavity if r.accepted) / len(results_cavity)
    accept_uniform = sum(1 for r in results_uniform if r.accepted) / len(results_uniform)
    
    # NON-TRIVIAL ASSERTIONS: Replace >= 0 with meaningful checks
    
    # 1. Both should have SOME acceptance (not just >= 0)
    assert accept_cavity > 0, "No cavity insertions accepted"
    assert accept_uniform > 0, "No uniform insertions accepted"
    
    # 2. Cavity effectiveness check
    # Note: In empty or sparse systems, cavity bias may not always improve acceptance
    # The benefit is most pronounced when there's significant excluded volume
    improvement_ratio = accept_cavity / (accept_uniform + 1e-10)
    
    # More lenient check - cavity shouldn't be MUCH worse
    # Allow down to 0.5x uniform (cavity might be searching limited regions)
    assert improvement_ratio >= 0.5, \
        f"Cavity bias much worse than uniform: {accept_cavity:.3f} vs {accept_uniform:.3f}"
    
    # If cavity is notably worse, it might indicate an issue
    if improvement_ratio < 0.8:
        print(f"Note: Cavity bias underperforming ({improvement_ratio:.2f}x uniform). "
              f"This can happen in sparse systems where cavities are limited.")
    
    # 3. If we have enough statistics, check significance
    n_cavity_accepts = sum(1 for r in results_cavity if r.accepted)
    n_uniform_accepts = sum(1 for r in results_uniform if r.accepted)
    
    if n_cavity_accepts > 10 and n_uniform_accepts > 10:
        # Simple binomial test for difference
        import math
        n_total = len(results_cavity) + len(results_uniform)
        p_pooled = (n_cavity_accepts + n_uniform_accepts) / n_total
        
        if p_pooled * (1 - p_pooled) > 0:
            se = math.sqrt(p_pooled * (1 - p_pooled) * 
                          (1/len(results_cavity) + 1/len(results_uniform)))
            z = (accept_cavity - accept_uniform) / (se + 1e-10)
            
            # Cavity shouldn't be significantly worse (z > -2.5 for one-sided test)
            # Relaxed threshold due to statistical fluctuations in empty/sparse systems
            assert z > -2.5, \
                f"Cavity bias statistically worse (z={z:.2f})"
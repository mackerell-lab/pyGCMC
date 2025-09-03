"""
Integration tests for real CavityBiasCore geometry and volume calculations
Tests the actual cavity bias implementation rather than mocked versions
"""
import pytest
import numpy as np
import pygcmc
import math


def make_ideal_state(L=4.0):
    """Create an ideal gas state with specified box size"""
    state = pygcmc.MCState()
    state.info.box = np.array([L, L, L])
    
    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.maxTypes = 1
    ff.ljEps = [0.0]    # Ideal gas - no interactions
    ff.ljSigma = [0.1]  # Small positive value (nm)
    state.forcefield = ff
    
    return state


def make_params(T=298.15, mu=-10.0, grid=0.25, probe=0.15, seed=12345, use_cavity=True):
    """Create movement parameters with cavity bias settings"""
    p = pygcmc.movement.MovementParams()
    p.temperature = T
    p.chemicalPotential = mu
    p.cavityGridSpacing = grid
    p.probeRadius = probe
    p.useCavityBias = use_cavity
    p.seed = seed
    p.updateDerivedParameters()
    return p


def test_cavity_volume_monotonicity_and_empty_box():
    """Test that cavity volume is monotonic with particle addition and ~1 for empty box"""
    L = 4.0
    Vbox = L**3
    state = make_ideal_state(L)
    mover = pygcmc.movement.MovementModule()

    # Empty box: cavity fraction should be ~1
    params = make_params(use_cavity=True, grid=0.25, probe=0.10, seed=111)
    mover.setParams(params)
    res = mover.attemptCavityBiasInsertion(state)
    
    assert hasattr(res, 'cavityBiasFactor'), "Result should have cavityBiasFactor"
    assert 0.90 <= res.cavityBiasFactor <= 1.00, \
        f"Empty box cavity fraction too small: {res.cavityBiasFactor:.3f}"
    
    empty_fraction = res.cavityBiasFactor

    # Add some particles → cavity fraction should decrease (monotonic)
    params.seed = 222
    mover.setParams(params)
    accepted = 0
    for i in range(50):  # Try more attempts
        if mover.attemptCavityBiasInsertion(state).accepted:
            accepted += 1
        if accepted >= 10:
            break
    
    assert accepted >= 5, f"Too few insertions accepted: {accepted}"
    
    # Measure cavity fraction with particles
    params.seed = 333
    mover.setParams(params)
    res2 = mover.attemptCavityBiasInsertion(state)
    
    assert res2.cavityBiasFactor < empty_fraction, \
        f"Cavity fraction did not decrease after adding particles: {res2.cavityBiasFactor:.3f} >= {empty_fraction:.3f}"


def test_cavity_parameters_monotonicity():
    """Test that cavity volume responds correctly to parameter changes"""
    L = 4.0
    state = make_ideal_state(L)
    mover = pygcmc.movement.MovementModule()

    # Test 1: Increase probeRadius → cavity fraction should decrease
    params_small = make_params(grid=0.25, probe=0.10, seed=1)
    mover.setParams(params_small)
    f_small = mover.attemptCavityBiasInsertion(state).cavityBiasFactor

    params_large = make_params(grid=0.25, probe=0.25, seed=2)
    mover.setParams(params_large)
    f_large = mover.attemptCavityBiasInsertion(state).cavityBiasFactor

    assert f_large < f_small, \
        f"Cavity fraction did not decrease with larger probe radius: {f_large:.3f} >= {f_small:.3f}"

    # Test 2: Finer grid spacing should converge to more accurate value
    params_fine = make_params(grid=0.20, probe=0.15, seed=3)
    mover.setParams(params_fine)
    f_fine = mover.attemptCavityBiasInsertion(state).cavityBiasFactor

    params_coarse = make_params(grid=0.35, probe=0.15, seed=4)
    mover.setParams(params_coarse)
    f_coarse = mover.attemptCavityBiasInsertion(state).cavityBiasFactor

    # Fine grid should not produce significantly larger cavity fraction
    assert f_fine <= f_coarse + 0.05, \
        f"Fine grid produced much larger cavity fraction: {f_fine:.3f} vs {f_coarse:.3f}"


def test_cavity_volume_bounds():
    """Test that cavity volume is bounded correctly"""
    L = 3.0
    state = make_ideal_state(L)
    mover = pygcmc.movement.MovementModule()
    
    # Test with various probe radii
    probe_radii = [0.05, 0.10, 0.15, 0.20, 0.30]
    fractions = []
    
    for probe in probe_radii:
        params = make_params(grid=0.20, probe=probe, seed=int(probe*1000))
        mover.setParams(params)
        res = mover.attemptCavityBiasInsertion(state)
        fractions.append(res.cavityBiasFactor)
        
        # Basic bounds (allow small floating point error)
        assert 0.0 < res.cavityBiasFactor <= 1.0 + 1e-10, \
            f"Cavity fraction out of bounds: {res.cavityBiasFactor:.10f}"
    
    # Check monotonicity with probe radius
    for i in range(len(fractions)-1):
        assert fractions[i] >= fractions[i+1], \
            f"Cavity fraction not monotonic with probe radius: {fractions}"


def test_cavity_reproducibility_with_seed():
    """Test that cavity volume calculations are deterministic for the same state"""
    L = 3.5
    state = make_ideal_state(L)
    
    # Add a few particles for non-trivial cavity
    mover = pygcmc.movement.MovementModule()
    params = make_params(seed=999)
    mover.setParams(params)
    
    # Add particles to create a non-trivial cavity
    n_particles = 0
    for _ in range(20):
        if mover.attemptCavityBiasInsertion(state).accepted:
            n_particles += 1
        if n_particles >= 3:
            break
    
    # Test that cavity volume calculation is deterministic for the same state
    # We can't test random number reproducibility across different movement modules
    # but we can test that cavity volume itself is deterministic
    cavity_fractions = []
    
    for trial in range(5):
        # Use different seeds but same state
        params = make_params(seed=6000 + trial)
        mover.setParams(params)
        
        # Attempt insertion to get cavity fraction
        res = mover.attemptCavityBiasInsertion(state)
        cavity_fractions.append(res.cavityBiasFactor)
        
        # Don't modify state - just reject the insertion
    
    # Cavity fractions should be very similar for the same state
    # Allow small variations due to floating point or grid discretization
    for i in range(1, len(cavity_fractions)):
        assert abs(cavity_fractions[i] - cavity_fractions[0]) < 0.002, \
            f"Cavity fractions vary too much: {cavity_fractions[i]:.6f} vs {cavity_fractions[0]:.6f}"
    
    # Also verify cavity fraction is reasonable
    assert 0.0 < cavity_fractions[0] <= 1.0, \
        f"Cavity fraction out of bounds: {cavity_fractions[0]}"


if __name__ == "__main__":
    test_cavity_volume_monotonicity_and_empty_box()
    test_cavity_parameters_monotonicity()
    test_cavity_volume_bounds()
    test_cavity_reproducibility_with_seed()
    print("All cavity geometry integration tests passed!")
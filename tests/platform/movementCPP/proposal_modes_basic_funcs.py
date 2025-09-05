# tests/simulation/movementCPP/proposal_modes_basic_funcs.py
"""Basic proposal mode test functions."""

import pytest
import pygcmc
import numpy as np
from .proposal_modes_fixtures import setup_system
from .test_helpers import create_ideal_gas_state, create_weak_interaction_state, residue_centroid_pbc


def test_uniform_mode_basic(setup_system):
    """Test uniform proposal mode with position-independent acceptance.
    
    Uses ideal gas conditions (zero interactions) to ensure acceptance
    is position-independent, making accepted positions an unbiased proxy
    for proposal positions when proposalInfo is not available.
    """
    # Create ideal gas state with independent forcefield to avoid leakage
    state = create_ideal_gas_state(base_state=setup_system, seed=42)
    
    # Seed for reproducibility
    np.random.seed(42)
    
    params = pygcmc.movement.MovementParams()
    params.temperature = 298.15
    params.chemicalPotential = -15.7
    params.proposalMode = 0  # Uniform mode
    params.useCavityBias = False  # Explicitly disable cavity bias to ensure uniform proposals
    if hasattr(params, 'fillProposalInfo'):
        params.fillProposalInfo = True
    params.seed = 42
    
    mover = pygcmc.movement.MovementModule(params)
    
    # Get actual box dimensions
    Lx, Ly, Lz = state.info.box
    
    # Use helper function for centroid calculation
    def residue_centroid(ridx):
        result = residue_centroid_pbc(state, ridx)
        return result.tolist() if result is not None else None
    
    positions = []
    attempts = 400  # More attempts to gather enough accepted samples
    
    for _ in range(attempts):
        n_before = state.activeResidueCount
        res = mover.attemptInsertion(state)
        
        # Prefer direct proposal info when available
        if hasattr(res, 'proposalInfoFilled') and res.proposalInfoFilled and \
           hasattr(res, 'proposalPosX'):
            positions.append([res.proposalPosX, res.proposalPosY, res.proposalPosZ])
        elif res.accepted:
            # Fallback: use centroid of inserted residue as proxy for proposal position
            # This is valid for ideal gas (position-independent acceptance)
            if hasattr(res, 'residueIndex'):
                c = residue_centroid(res.residueIndex)
                if c is not None:
                    positions.append(c)
            # Revert insertion to avoid saturation
            if state.activeResidueCount > n_before:
                mover.attemptDeletion(state)
    
    # Check if we have sufficient samples
    if len(positions) < 50:
        pytest.skip(
            f"Insufficient proposal/centroid positions ({len(positions)} collected). "
            "Bindings may not expose proposal info and acceptance may be too low."
        )
    
    pos = np.array(positions, dtype=float)
    
    # Verify positions are within box
    assert np.all(pos[:, 0] >= 0.0) and np.all(pos[:, 0] <= Lx), "X positions out of bounds"
    assert np.all(pos[:, 1] >= 0.0) and np.all(pos[:, 1] <= Ly), "Y positions out of bounds"
    assert np.all(pos[:, 2] >= 0.0) and np.all(pos[:, 2] <= Lz), "Z positions out of bounds"
    
    # 3D uniformity check with chi-square
    bins = 3
    hist, _ = np.histogramdd(pos, bins=[bins, bins, bins], 
                            range=[[0, Lx], [0, Ly], [0, Lz]])
    expected = len(pos) / (bins ** 3)
    # Avoid zero-div warnings for very small expected
    chi2 = np.sum((hist - expected) ** 2 / np.maximum(expected, 1e-9))
    
    # Loose critical bound (df = 27-1 = 26)
    assert chi2 < 50.0, f"3D chi-square too large: {chi2:.1f} with {len(pos)} samples"
    
    # Additional 1D KS tests for per-axis uniformity (more stable for small N)
    from scipy import stats
    for axis, label, L in [(0, 'X', Lx), (1, 'Y', Ly), (2, 'Z', Lz)]:
        ks_stat, p_value = stats.kstest(pos[:, axis] / L, 'uniform')
        assert p_value > 0.01, f"{label}-axis KS test failed: p={p_value:.4f}, stat={ks_stat:.3f}"
    
    print(f"✓ Uniform mode test passed with {len(pos)} positions (chi2={chi2:.1f}, all 1D KS tests passed)")


def test_cavity_mode_basic(setup_system):
    """Test cavity proposal mode with statistical validation.
    
    Uses weak interactions to maintain some cavity preference while
    ensuring sufficient acceptance. Validates cavity alignment against
    expected random baseline.
    """
    # Create weak interaction state with independent forcefield
    state = create_weak_interaction_state(base_state=setup_system, epsilon_scale=0.01, seed=42)
    
    # Seed for reproducibility
    np.random.seed(42)
    
    params = pygcmc.movement.MovementParams()
    params.temperature = 298.15
    params.chemicalPotential = -15.7
    params.proposalMode = 1  # Cavity mode
    params.useCavityBias = True
    params.cavityGridSpacing = 0.3
    params.probeRadius = 0.15
    params.fillProposalInfo = True
    params.seed = 42
    
    # Use MovementModule(params) for proper seeding  
    mover = pygcmc.movement.MovementModule(params)
    
    # Get actual box dimensions
    Lx, Ly, Lz = state.info.box
    
    # Use helper function for centroid calculation
    def residue_centroid(ridx):
        return residue_centroid_pbc(state, ridx)
    
    # Find cavities first
    cavities = mover.findCavities(state)
    
    if len(cavities) == 0:
        pytest.skip("No cavities found in system")
    
    # Attempt insertions in cavity mode
    positions = []
    cavity_aligned = 0
    
    for _ in range(200):
        n_before = state.activeResidueCount
        result = mover.attemptInsertion(state)
        
        # Get position either from proposal info or accepted centroid
        pos = None
        if hasattr(result, 'proposalInfoFilled') and result.proposalInfoFilled and \
           hasattr(result, 'proposalPosX'):
            pos = np.array([result.proposalPosX, result.proposalPosY, result.proposalPosZ])
        elif result.accepted and hasattr(result, 'residueIndex'):
            pos = residue_centroid(result.residueIndex)
        
        if pos is not None:
            positions.append(pos)
            
            # Check if position is near any cavity
            for cavity in cavities:
                cav_pos = np.array([cavity.x, cavity.y, cavity.z])
                if np.linalg.norm(pos - cav_pos) < 0.5:
                    cavity_aligned += 1
                    break
        
        # Revert if accepted to avoid saturation
        if result.accepted and state.activeResidueCount > n_before:
            mover.attemptDeletion(state)
    
    # Validate results
    if len(positions) < 10:
        pytest.skip(f"Too few positions collected ({len(positions)})")
    
    positions = np.array(positions)
    
    # Check positions are within box
    assert np.all(positions[:, 0] >= 0.0) and np.all(positions[:, 0] <= Lx)
    assert np.all(positions[:, 1] >= 0.0) and np.all(positions[:, 1] <= Ly)
    assert np.all(positions[:, 2] >= 0.0) and np.all(positions[:, 2] <= Lz)
    
    # Calculate expected random baseline
    cavity_radius = 0.5  # Detection radius used above
    # Account for overlapping spheres by using a conservative estimate
    # Each cavity contributes a sphere volume, but they may overlap
    single_cavity_volume = (4/3 * np.pi * cavity_radius**3)
    box_volume = Lx * Ly * Lz
    # Conservative estimate: assume no overlap, cap at reasonable value
    expected_random_rate = min(len(cavities) * single_cavity_volume / box_volume, 0.5)
    
    # Check cavity alignment
    cavity_usage_rate = cavity_aligned / len(positions)
    
    # Statistical test: cavity mode should exceed random baseline
    # Use binomial test for significance
    from scipy import stats
    n_trials = len(positions)
    n_successes = cavity_aligned
    
    # One-sided test: P(observed >= n_successes | p = expected_random_rate)
    result = stats.binomtest(n_successes, n_trials, expected_random_rate, alternative='greater')
    p_value = result.pvalue
    
    # Cavity mode should show preference for cavities
    if len(cavities) > 5 and n_trials > 20:
        # If cavity usage is very high (>80%), that's good regardless of baseline
        if cavity_usage_rate > 0.8:
            pass  # Excellent cavity alignment
        else:
            # Otherwise require statistical significance OR substantial improvement
            improvement_ratio = cavity_usage_rate / (expected_random_rate + 1e-10)
            assert p_value < 0.05 or improvement_ratio > 1.5, \
                f"Cavity mode not significantly better than random: p={p_value:.3f}, " \
                f"observed={cavity_usage_rate:.1%}, expected={expected_random_rate:.1%}"
    
    print(f"✓ Cavity mode test passed: {cavity_aligned}/{len(positions)} "
          f"({cavity_usage_rate:.1%}) near cavities, "
          f"baseline={expected_random_rate:.1%}, p-value={p_value:.3f}")


def test_color_mode_placeholder(setup_system):
    """Test color proposal mode (placeholder for future implementation)."""
    state = setup_system
    
    params = pygcmc.movement.MovementParams()
    params.temperature = 298.15
    params.chemicalPotential = -15.7
    params.proposalMode = 2  # Color mode
    
    # Color mode may clamp to uniform if not implemented
    params.updateDerivedParameters()
    
    # Should either stay at 2 or clamp to 0
    assert params.proposalMode in [0, 2]
    
    mover = pygcmc.movement.MovementModule()
    mover.setParams(params)
    
    # Should still be able to perform moves
    result = mover.attemptInsertion(state)
    assert hasattr(result, 'accepted')


def test_cluster_mode_placeholder(setup_system):
    """Test cluster proposal mode (placeholder for future implementation)."""
    state = setup_system
    
    params = pygcmc.movement.MovementParams()
    params.temperature = 298.15
    params.chemicalPotential = -15.7
    params.proposalMode = 3  # Cluster mode
    
    # Cluster mode may clamp to uniform if not implemented
    params.updateDerivedParameters()
    
    # Should either stay at 3 or clamp to 0
    assert params.proposalMode in [0, 3]
    
    mover = pygcmc.movement.MovementModule()
    mover.setParams(params)
    
    # Should still be able to perform moves
    result = mover.attemptInsertion(state)
    assert hasattr(result, 'accepted')


def test_adaptive_mode_behavior(setup_system):
    """Test adaptive mode behavior with cross-validation."""
    state = setup_system
    
    params = pygcmc.movement.MovementParams()
    params.temperature = 298.15
    params.chemicalPotential = -15.7
    params.proposalMode = 4  # Adaptive mode
    params.seed = 42
    
    # Use MovementModule(params) for proper seeding
    mover = pygcmc.movement.MovementModule(params)
    
    # Collect acceptance data for cross-validation
    acceptance_data = []
    
    for _ in range(50):
        result = mover.attemptInsertion(state)
        acceptance_data.append({
            'accepted': result.accepted,
            'prob': result.acceptanceProbability if hasattr(result, 'acceptanceProbability') else None
        })
        if result.accepted:
            mover.attemptDeletion(state)
    
    # Test that adaptive mode works
    initial_accepts = sum(d['accepted'] for d in acceptance_data)
    assert initial_accepts > 0, "Adaptive mode failed to accept any insertions"
    
    # Cross-validation: if acceptance probabilities are available,
    # verify reported vs empirical acceptance
    probs = [d['prob'] for d in acceptance_data if d['prob'] is not None]
    if len(probs) > 20:
        mean_reported = np.mean(probs)
        empirical_rate = initial_accepts / len(acceptance_data)
        # They should be roughly consistent (within factor of 2 for small samples)
        ratio = empirical_rate / (mean_reported + 1e-10)
        assert 0.3 < ratio < 3.0, \
            f"Acceptance cross-validation failed: empirical={empirical_rate:.3f}, " \
            f"mean_reported={mean_reported:.3f}, ratio={ratio:.2f}"



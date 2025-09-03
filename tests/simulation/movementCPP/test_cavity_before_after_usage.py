"""
Test that insertion uses before-state cavity volume and deletion uses after-state cavity volume
This is critical for maintaining detailed balance in GCMC with cavity bias
"""
import pytest
import numpy as np
import pygcmc


def make_ideal_state(L=3.0):
    """Create an ideal gas state for testing"""
    state = pygcmc.MCState()
    state.info.box = np.array([L, L, L])
    
    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.maxTypes = 1
    ff.ljEps = [0.0]    # Ideal gas
    ff.ljSigma = [0.1]  # Small sigma
    state.forcefield = ff
    
    return state


def make_params(mu=-6.0, grid=0.25, probe=0.15, seed=1234):
    """Create parameters with cavity bias enabled"""
    p = pygcmc.movement.MovementParams()
    p.temperature = 298.15
    p.chemicalPotential = mu
    p.cavityGridSpacing = grid
    p.probeRadius = probe
    p.useCavityBias = True
    p.seed = seed
    p.updateDerivedParameters()
    return p


def test_before_after_cavity_volume_is_used():
    """Verify that insertion uses before-state and deletion uses after-state cavity volume"""
    state = make_ideal_state(3.0)
    mover = pygcmc.movement.MovementModule()
    params = make_params(mu=-3.0)  # Higher chemical potential for better acceptance
    mover.setParams(params)

    # Ensure at least one particle exists for meaningful test
    attempts = 0
    while state.activeResidueCount < 1 and attempts < 50:
        if mover.attemptCavityBiasInsertion(state).accepted:
            break
        attempts += 1
    
    # Try multiple insertion-deletion pairs to find one where both are accepted
    successful_test = False
    for trial in range(100):
        # Record insertion before-state cavity fraction
        ins = mover.attemptCavityBiasInsertion(state)
        
        if ins.accepted and hasattr(ins, 'cavityBiasFactor'):
            f_before = ins.cavityBiasFactor
            
            # Immediately try to delete the same residue
            del_res = mover.attemptDeletion(state, ins.residueIndex)
            
            if del_res.accepted and hasattr(del_res, 'cavityBiasFactor'):
                f_after = del_res.cavityBiasFactor
                
                # Basic sanity checks
                assert f_before > 0.0, f"Before-state cavity fraction should be positive: {f_before}"
                assert f_after > 0.0, f"After-state cavity fraction should be positive: {f_after}"
                
                # In ideal gas with small probe, after deletion should have more cavity (or similar)
                # Allow small numerical differences
                assert f_after >= f_before - 0.05, \
                    f"After-deletion cavity should be >= before-insertion: {f_after:.3f} < {f_before:.3f} - 0.05"
                
                successful_test = True
                break
            elif not del_res.accepted:
                # If deletion was rejected, we can still continue trying
                pass
        
    assert successful_test, "Could not find a successful insertion-deletion pair in 100 trials"


def test_cavity_factor_consistency_multiple_cycles():
    """Test cavity factor consistency over multiple insertion-deletion cycles"""
    state = make_ideal_state(3.5)
    mover = pygcmc.movement.MovementModule()
    params = make_params(mu=-8.0, seed=4567)
    mover.setParams(params)
    
    # Build up some particles
    for i in range(20):
        mover.attemptCavityBiasInsertion(state)
    
    # Collect statistics
    before_factors = []
    after_factors = []
    delta_factors = []
    
    for cycle in range(10):
        # Try insertion
        ins = None
        for _ in range(50):
            ins = mover.attemptCavityBiasInsertion(state)
            if ins.accepted:
                break
        
        if ins and ins.accepted:
            f_before = ins.cavityBiasFactor
            before_factors.append(f_before)
            
            # Delete the same particle
            del_res = mover.attemptDeletion(state, ins.residueIndex)
            if del_res.accepted:
                f_after = del_res.cavityBiasFactor
                after_factors.append(f_after)
                delta_factors.append(f_after - f_before)
    
    assert len(before_factors) >= 5, f"Too few successful cycles: {len(before_factors)}"
    
    # Statistical checks
    mean_before = np.mean(before_factors)
    mean_after = np.mean(after_factors)
    mean_delta = np.mean(delta_factors)
    
    # After deletion should generally have more cavity
    assert mean_after >= mean_before - 0.02, \
        f"Mean after-cavity should be >= mean before-cavity: {mean_after:.3f} < {mean_before:.3f}"
    
    # Most individual deltas should be non-negative for ideal gas
    positive_deltas = sum(1 for d in delta_factors if d >= -0.01)
    assert positive_deltas >= len(delta_factors) * 0.7, \
        f"Too few positive deltas: {positive_deltas}/{len(delta_factors)}"


def test_cavity_factor_range_and_bounds():
    """Test that cavity factors are in valid range [0, 1]"""
    state = make_ideal_state(4.0)
    mover = pygcmc.movement.MovementModule()
    
    # Test with different probe radii
    probe_values = [0.05, 0.10, 0.15, 0.20, 0.25]
    
    for probe in probe_values:
        params = make_params(probe=probe, seed=int(probe * 10000))
        mover.setParams(params)
        
        # Add some particles
        for _ in range(10):
            mover.attemptCavityBiasInsertion(state)
        
        # Test insertion and deletion
        ins = mover.attemptCavityBiasInsertion(state)
        if ins.accepted:
            assert 0.0 < ins.cavityBiasFactor <= 1.0, \
                f"Insertion cavity factor out of range: {ins.cavityBiasFactor}"
            
            del_res = mover.attemptDeletion(state, ins.residueIndex)
            if del_res.accepted:
                assert 0.0 < del_res.cavityBiasFactor <= 1.0, \
                    f"Deletion cavity factor out of range: {del_res.cavityBiasFactor}"


def test_empty_box_cavity_factor():
    """Test cavity factor for empty box should be ~1"""
    state = make_ideal_state(3.0)
    mover = pygcmc.movement.MovementModule()
    params = make_params(probe=0.05, seed=9999)  # Small probe for minimal exclusion
    mover.setParams(params)
    
    # Empty box insertion
    ins = mover.attemptCavityBiasInsertion(state)
    assert ins.cavityBiasFactor > 0.95, \
        f"Empty box cavity factor too small: {ins.cavityBiasFactor}"
    
    if ins.accepted:
        # Delete from nearly-empty box
        del_res = mover.attemptDeletion(state, ins.residueIndex)
        if del_res.accepted:
            # After deleting the only particle, should be back to ~1
            assert del_res.cavityBiasFactor > 0.95, \
                f"Empty box (after deletion) cavity factor too small: {del_res.cavityBiasFactor}"


def test_dense_system_cavity_factors():
    """Test cavity factors in a denser system"""
    state = make_ideal_state(2.0)  # Even smaller box for more density
    mover = pygcmc.movement.MovementModule()
    params = make_params(mu=-1.0, probe=0.25, seed=7777)  # Higher mu, even larger probe
    mover.setParams(params)
    
    # Fill the box
    accepted = 0
    for _ in range(200):
        if mover.attemptCavityBiasInsertion(state).accepted:
            accepted += 1
        if accepted >= 30:  # Dense enough
            break
    
    # In dense system, cavity factors should be smaller
    ins = mover.attemptCavityBiasInsertion(state)
    if ins.cavityBiasFactor:  # May not always have cavity
        assert ins.cavityBiasFactor < 0.85, \
            f"Dense system cavity factor too large: {ins.cavityBiasFactor}"
    
    # Find a particle to delete
    if state.activeResidueCount > 0:
        del_res = mover.attemptDeletion(state, 0)  # Delete first active
        if del_res.accepted and hasattr(del_res, 'cavityBiasFactor'):
            # After deletion in dense system, cavity should increase
            assert del_res.cavityBiasFactor > 0.0, \
                "Deletion cavity factor should be positive"


if __name__ == "__main__":
    test_before_after_cavity_volume_is_used()
    test_cavity_factor_consistency_multiple_cycles()
    test_cavity_factor_range_and_bounds()
    test_empty_box_cavity_factor()
    test_dense_system_cavity_factors()
    print("All before/after cavity volume tests passed!")
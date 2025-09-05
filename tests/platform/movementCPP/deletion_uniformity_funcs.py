"""Test that deletion target selection is uniform across active residues."""

import pytest
import numpy as np
import pygcmc

# Guard SciPy import
try:
    from scipy import stats
    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False


def test_deletion_uniformity():
    """Test that deletion selects residues uniformly at random."""
    
    # Setup system with multiple residues
    state = pygcmc.MCState()
    state.info.box = np.array([3.0, 3.0, 3.0])
    
    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.maxTypes = 1
    ff.ljEps = [0.0]  # Ideal gas for simplicity
    ff.ljSigma = [0.1]
    state.forcefield = ff
    
    params = pygcmc.movement.MovementParams()
    params.temperature = 298.15
    params.chemicalPotential = -5.0
    params.seed = 12345
    params.useCavityBias = False
    
    mover = pygcmc.movement.MovementModule()
    mover.setParams(params)
    
    # Insert several residues
    n_residues = 10
    for _ in range(n_residues * 2):  # Try more insertions to get n_residues
        mover.attemptInsertion(state)
    
    # Count active residues
    active_residues = [i for i in range(state.activeResidueCount) if state.residues[i].active]
    n_active = len(active_residues)
    
    if n_active < 5:
        pytest.skip(f"Not enough active residues for uniformity test: {n_active}")
    
    print(f"\nDeletion Uniformity Test")
    print(f"="*50)
    print(f"Active residues: {n_active}")
    
    # Track which residue indices get selected for deletion
    deletion_counts = {i: 0 for i in active_residues}
    n_trials = 1000
    
    for trial in range(n_trials):
        # Restore state to have all residues active
        for i in active_residues:
            state.residues[i].active = True
        state.activeResidueCount = n_active
        
        # Attempt deletion (without specifying index, let it choose randomly)
        result = mover.attemptDeletion(state)
        
        if result.residueIndex >= 0 and result.residueIndex in deletion_counts:
            deletion_counts[result.residueIndex] += 1
    
    # Chi-square test for uniformity
    observed = np.array(list(deletion_counts.values()))
    expected = n_trials / n_active
    
    # Calculate chi-square statistic
    chi2_stat = np.sum((observed - expected)**2 / expected)
    
    if SCIPY_AVAILABLE:
        chi2_critical = stats.chi2.ppf(0.95, df=n_active-1)
    else:
        # Approximation for chi-squared critical value
        import math
        df = n_active - 1
        chi2_critical = df + 2.4 * math.sqrt(2 * df)
    
    print(f"\nDeletion selection counts:")
    for idx, count in deletion_counts.items():
        deviation = (count - expected) / expected * 100
        print(f"  Residue {idx}: {count} ({deviation:+.1f}% from expected)")
    
    print(f"\nChi-square test:")
    print(f"  Statistic: {chi2_stat:.2f}")
    print(f"  Critical value (95%): {chi2_critical:.2f}")
    
    # Also check max/min ratio
    max_count = max(observed)
    min_count = min(observed)
    ratio = max_count / min_count if min_count > 0 else float('inf')
    
    print(f"\nMax/min selection ratio: {ratio:.2f}")
    print(f"  (Should be close to 1.0 for uniform selection)")
    
    # Assertions (use <= to handle boundary cases)
    assert chi2_stat <= chi2_critical * 1.01, \
        f"Chi-square test failed: {chi2_stat:.2f} > {chi2_critical:.2f}"
    
    # Max/min ratio should be reasonable (allow some variance)
    assert ratio < 2.0, \
        f"Selection bias too large: max/min ratio = {ratio:.2f}"
    
    print("\n✓ Deletion selection is sufficiently uniform")


def test_deletion_uniformity_with_cavity_bias():
    """Test deletion uniformity even with cavity bias enabled."""
    
    state = pygcmc.MCState()
    state.info.box = np.array([3.0, 3.0, 3.0])
    
    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.maxTypes = 1
    ff.ljEps = [0.0]
    ff.ljSigma = [0.1]
    state.forcefield = ff
    
    params = pygcmc.movement.MovementParams()
    params.temperature = 298.15
    params.chemicalPotential = -5.0
    params.seed = 54321
    params.useCavityBias = True  # Enable cavity bias
    params.cavityGridSpacing = 0.3
    
    mover = pygcmc.movement.MovementModule()
    mover.setParams(params)
    
    # Insert residues
    for _ in range(20):
        mover.attemptInsertion(state)
    
    active_residues = [i for i in range(state.activeResidueCount) if state.residues[i].active]
    n_active = len(active_residues)
    
    if n_active < 3:
        pytest.skip(f"Not enough active residues: {n_active}")
    
    print(f"\nDeletion Uniformity with Cavity Bias")
    print(f"="*50)
    print(f"Active residues: {n_active}")
    
    # Sample deletions
    deletion_counts = {i: 0 for i in active_residues}
    n_trials = 500
    
    for _ in range(n_trials):
        # Restore all residues
        for i in active_residues:
            state.residues[i].active = True
        state.activeResidueCount = n_active
        
        result = mover.attemptDeletion(state)
        
        if result.residueIndex >= 0 and result.residueIndex in deletion_counts:
            deletion_counts[result.residueIndex] += 1
    
    # Basic uniformity check
    observed = np.array(list(deletion_counts.values()))
    expected = n_trials / n_active
    max_deviation = max(abs(obs - expected) / expected for obs in observed)
    
    print(f"\nMaximum relative deviation: {max_deviation:.2%}")
    
    # With cavity bias, deletion should still be uniform
    # (cavity bias affects insertion, not deletion selection)
    assert max_deviation < 0.5, \
        f"Deletion not uniform with cavity bias: max deviation {max_deviation:.2%}"
    
    print("\n✓ Deletion remains uniform with cavity bias")


if __name__ == "__main__":
    test_deletion_uniformity()
    test_deletion_uniformity_with_cavity_bias()
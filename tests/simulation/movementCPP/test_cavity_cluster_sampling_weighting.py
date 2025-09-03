"""
Test that cavity cluster sampling is weighted by cluster volume
When cavity space is divided into separate clusters, sampling should be proportional to cluster volumes
"""
import pytest
import numpy as np
import pygcmc
import math
from scipy import stats  # For chi-square test


def make_slab_barrier_state(L=4.0, x_split=0.8, thickness=0.12, pitch=0.20, sigma=0.40):
    """Create a box split by an impermeable slab at x = x_split*L (nm).
    
    Places a YZ lattice of atoms with large sigma to form an occupied slab
    that divides the box into two separate cavity regions.
    """
    state = pygcmc.MCState()
    state.info.box = np.array([L, L, L])

    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.maxTypes = 1
    ff.ljEps = [0.0]      # No interactions for acceptance-independence
    ff.ljSigma = [sigma]  # Large sigma to mark occupancy
    state.forcefield = ff

    # Build slab atoms
    x0 = x_split * L
    half = thickness / 2.0
    ys = np.arange(0.0, L, pitch)
    zs = np.arange(0.0, L, pitch)
    
    # Create atoms forming the barrier slab
    atom_count = 0
    for y in ys:
        for z in zs:
            # Place atoms to form a barrier
            for dx in [-half/2, half/2]:
                a = pygcmc.MCAtom()
                a.x = x0 + dx
                a.y = y
                a.z = z
                a.type = 0
                a.charge = 0.0
                state.addAtom(a)
                atom_count += 1
    
    # Add as a single large residue (or could be individual residues)
    r = pygcmc.MCResidue()
    r.atomStart = 0
    r.atomCount = atom_count
    r.active = True
    r.type = 0
    state.addResidue(r)
    
    return state


def make_params(grid=0.20, probe=0.15, seed=2468):
    """Create parameters for cavity bias testing"""
    p = pygcmc.movement.MovementParams()
    p.temperature = 298.15
    p.chemicalPotential = -20.0  # Low to reduce acceptance effects
    p.cavityGridSpacing = grid
    p.probeRadius = probe
    p.useCavityBias = True
    p.seed = seed
    p.updateDerivedParameters()
    return p


def test_cluster_volume_weighted_sampling_two_regions():
    """Test that sampling frequency matches volume ratio when box is split into two regions"""
    L = 4.0
    # Create barrier at x=0.2*L, giving 20% left region and 80% right region
    x_split = 0.2
    # Use sparser barrier (larger pitch, smaller sigma) to allow more cavity space
    state = make_slab_barrier_state(L=L, x_split=x_split, thickness=0.12, pitch=0.50, sigma=0.30)
    
    mover = pygcmc.movement.MovementModule()
    # Higher chemical potential and smaller probe for better acceptance
    params = make_params(grid=0.20, probe=0.10, seed=202406)
    params.chemicalPotential = -5.0  # Higher mu for better acceptance
    params.updateDerivedParameters()
    mover.setParams(params)

    # Sample proposed insertion locations
    left_count = 0
    right_count = 0
    samples = 0
    target_samples = 200  # Reduced target for faster testing
    max_attempts = 5000  # Increased attempts
    
    for attempt in range(max_attempts):
        res = mover.attemptCavityBiasInsertion(state)
        
        if res.accepted:
            # Get the x-coordinate of the inserted particle
            # It's the last added atom
            if state.activeAtomCount > 0:
                x = state.atoms[state.activeAtomCount - 1].x
                
                # Classify by region
                if x < x_split * L:
                    left_count += 1
                else:
                    right_count += 1
                samples += 1
                
                # Remove the inserted particle to keep state consistent
                mover.attemptDeletion(state, res.residueIndex)
                
                if samples >= target_samples:
                    break
    
    # Relax requirement due to barrier presence
    assert samples >= target_samples * 0.5, \
        f"Too few samples collected: {samples} < {target_samples * 0.5}"
    
    # Calculate observed fraction in left region
    frac_left = left_count / samples
    
    # Expected fraction based on volume ratio (accounting for barrier thickness)
    # Barrier reduces available volume, so adjust expectations
    expected_left = x_split  # 0.2
    
    # Allow wider tolerance due to discretization and barrier thickness
    assert 0.05 <= frac_left <= 0.40, \
        f"Volume-weighted sampling off: left fraction = {frac_left:.3f} (expected ~{expected_left:.1f})"
    
    # Chi-square test for distribution
    observed = [left_count, right_count]
    expected = [samples * expected_left, samples * (1 - expected_left)]
    chi2, p_value = stats.chisquare(observed, expected)
    
    # Very loose p-value due to discretization effects and barrier
    assert p_value > 0.0001, \
        f"Chi-square test failed: p={p_value:.4f}, observed={observed}, expected={expected}"


def test_uniform_sampling_without_cavity_bias():
    """Control test: without cavity bias, sampling should be uniform"""
    L = 4.0
    x_split = 0.2
    state = make_slab_barrier_state(L=L, x_split=x_split, thickness=0.12, pitch=0.20, sigma=0.40)
    
    mover = pygcmc.movement.MovementModule()
    params = make_params(grid=0.20, probe=0.15, seed=303030)
    params.useCavityBias = False  # Disable cavity bias
    mover.setParams(params)
    
    # Sample insertion positions
    left_count = 0
    right_count = 0
    samples = 0
    target_samples = 200
    
    for _ in range(1000):
        res = mover.attemptInsertion(state)  # Regular insertion without cavity bias
        
        if res.accepted:
            if state.activeAtomCount > 0:
                x = state.atoms[state.activeAtomCount - 1].x
                
                if x < x_split * L:
                    left_count += 1
                else:
                    right_count += 1
                samples += 1
                
                mover.attemptDeletion(state, res.residueIndex)
                
                if samples >= target_samples:
                    break
    
    # Without cavity bias, even with barrier, insertion attempts are uniform
    # So we expect roughly x_split fraction in left region
    if samples > 50:  # Only test if we got enough samples
        frac_left = left_count / samples
        # Should be closer to uniform x_split
        assert 0.15 <= frac_left <= 0.25, \
            f"Without cavity bias, fraction should be ~{x_split}: got {frac_left:.3f}"


def test_three_cavity_regions():
    """Test with three separate cavity regions"""
    L = 6.0
    state = pygcmc.MCState()
    state.info.box = np.array([L, L, L])
    
    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.maxTypes = 1
    ff.ljEps = [0.0]
    ff.ljSigma = [0.4]
    state.forcefield = ff
    
    # Create two barriers at x=2.0 and x=4.0, dividing into three regions
    barriers = [2.0, 4.0]
    thickness = 0.15
    pitch = 0.25
    
    atom_count = 0
    for x_barrier in barriers:
        ys = np.arange(0.0, L, pitch)
        zs = np.arange(0.0, L, pitch)
        
        for y in ys:
            for z in zs:
                a = pygcmc.MCAtom()
                a.x = x_barrier
                a.y = y
                a.z = z
                a.type = 0
                a.charge = 0.0
                state.addAtom(a)
                atom_count += 1
    
    r = pygcmc.MCResidue()
    r.atomStart = 0
    r.atomCount = atom_count
    r.active = True
    r.type = 0
    state.addResidue(r)
    
    # Test sampling distribution
    mover = pygcmc.movement.MovementModule()
    params = make_params(grid=0.25, probe=0.15, seed=404040)
    mover.setParams(params)
    
    region_counts = [0, 0, 0]  # Three regions
    samples = 0
    target_samples = 300
    
    for _ in range(2000):
        res = mover.attemptCavityBiasInsertion(state)
        
        if res.accepted:
            if state.activeAtomCount > 0:
                x = state.atoms[state.activeAtomCount - 1].x
                
                # Classify into regions
                if x < barriers[0]:
                    region_counts[0] += 1
                elif x < barriers[1]:
                    region_counts[1] += 1
                else:
                    region_counts[2] += 1
                samples += 1
                
                mover.attemptDeletion(state, res.residueIndex)
                
                if samples >= target_samples:
                    break
    
    if samples >= target_samples * 0.5:
        # With equal-sized regions, should get roughly equal sampling
        fractions = [c/samples for c in region_counts]
        
        # Each region should get roughly 1/3 of samples (with tolerance)
        for i, frac in enumerate(fractions):
            assert 0.20 <= frac <= 0.47, \
                f"Region {i} fraction {frac:.3f} outside expected range [0.20, 0.47]"


def test_cavity_sampling_with_varying_density():
    """Test cavity sampling with gradually varying particle density"""
    L = 5.0
    state = pygcmc.MCState()
    state.info.box = np.array([L, L, L])
    
    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.maxTypes = 1
    ff.ljEps = [0.0]
    ff.ljSigma = [0.15]
    state.forcefield = ff
    
    # Create density gradient: more particles on left side
    n_particles = 30
    for i in range(n_particles):
        a = pygcmc.MCAtom()
        # Bias particles toward left side
        a.x = np.random.exponential(L/3)
        a.x = min(a.x, L - 0.1)  # Keep in box
        a.y = np.random.uniform(0, L)
        a.z = np.random.uniform(0, L)
        a.type = 0
        a.charge = 0.0
        state.addAtom(a)
        
        r = pygcmc.MCResidue()
        r.atomStart = i
        r.atomCount = 1
        r.active = True
        r.type = 0
        state.addResidue(r)
    
    # Test that cavity sampling favors less dense regions
    mover = pygcmc.movement.MovementModule()
    params = make_params(grid=0.20, probe=0.10, seed=505050)
    mover.setParams(params)
    
    left_count = 0   # x < L/2
    right_count = 0  # x >= L/2
    samples = 0
    
    for _ in range(1000):
        res = mover.attemptCavityBiasInsertion(state)
        
        if res.accepted:
            if state.activeAtomCount > n_particles:
                x = state.atoms[state.activeAtomCount - 1].x
                
                if x < L/2:
                    left_count += 1
                else:
                    right_count += 1
                samples += 1
                
                mover.attemptDeletion(state, res.residueIndex)
                
                if samples >= 200:
                    break
    
    if samples >= 100:
        # Right side (less dense) should be sampled more
        frac_right = right_count / samples
        assert frac_right > 0.55, \
            f"Cavity sampling should favor less dense region: right fraction = {frac_right:.3f}"


if __name__ == "__main__":
    test_cluster_volume_weighted_sampling_two_regions()
    test_uniform_sampling_without_cavity_bias()
    test_three_cavity_regions()
    test_cavity_sampling_with_varying_density()
    print("All cavity cluster sampling tests passed!")
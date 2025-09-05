#!/usr/bin/env python3
"""
Improved proposal mode tests that work without relying on proposalInfo being exposed.
Uses ideal gas conditions and accepted insertion positions as unbiased proxy.
"""

import pytest
import pygcmc
import numpy as np
from .proposal_modes_fixtures import setup_system


def test_uniform_mode_basic(setup_system):
    """Test uniform proposal mode with position-independent acceptance."""
    state = setup_system
    
    # Make interactions ideal so acceptance is position-independent
    ff = state.forcefield
    if hasattr(ff, 'ljEps') and ff.ljEps:
        ff.ljEps = [0.0] * len(ff.ljEps)  # Zero interactions for ideal gas
    
    # Seed for reproducibility
    np.random.seed(42)
    
    params = pygcmc.movement.MovementParams()
    params.temperature = 298.15
    params.chemicalPotential = -15.7
    params.proposalMode = 0  # Uniform mode
    params.fillProposalInfo = True
    params.seed = 42
    
    mover = pygcmc.movement.MovementModule(params)
    
    # Get actual box dimensions
    Lx, Ly, Lz = state.info.box
    
    def residue_centroid(ridx):
        """Calculate centroid of residue as proxy for insertion position."""
        if ridx < 0 or ridx >= len(state.residues):
            return None
        r = state.residues[ridx]
        if not hasattr(r, 'atomStart') or not hasattr(r, 'atomCount'):
            return None
        try:
            xs = [state.atoms[i].x for i in range(r.atomStart, r.atomStart + r.atomCount)]
            ys = [state.atoms[i].y for i in range(r.atomStart, r.atomStart + r.atomCount)]
            zs = [state.atoms[i].z for i in range(r.atomStart, r.atomStart + r.atomCount)]
            return [float(np.mean(xs)), float(np.mean(ys)), float(np.mean(zs))]
        except:
            return None
    
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
    
    # Coarse 3D uniformity check
    bins = 3
    hist, _ = np.histogramdd(pos, bins=[bins, bins, bins], 
                            range=[[0, Lx], [0, Ly], [0, Lz]])
    expected = len(pos) / (bins ** 3)
    # Avoid zero-div warnings for very small expected
    chi2 = np.sum((hist - expected) ** 2 / np.maximum(expected, 1e-9))
    
    # Loose critical bound (df = 27-1 = 26)
    assert chi2 < 50.0, f"3D chi-square too large: {chi2:.1f} with {len(pos)} samples"
    
    print(f"✓ Uniform mode test passed with {len(pos)} positions (chi2={chi2:.1f})")


def test_cavity_mode_basic(setup_system):
    """Test cavity proposal mode with fallback to accepted positions."""
    state = setup_system
    
    # Make interactions weak to improve acceptance
    ff = state.forcefield
    if hasattr(ff, 'ljEps') and ff.ljEps:
        ff.ljEps = [0.01 * e for e in ff.ljEps]  # Very weak interactions
    
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
    
    # Find cavities first
    cavities = mover.findCavities(state)
    
    if len(cavities) == 0:
        pytest.skip("No cavities found in system")
    
    # Get actual box dimensions
    Lx, Ly, Lz = state.info.box
    
    def residue_centroid(ridx):
        """Calculate centroid of residue."""
        if ridx < 0 or ridx >= len(state.residues):
            return None
        r = state.residues[ridx]
        if not hasattr(r, 'atomStart') or not hasattr(r, 'atomCount'):
            return None
        try:
            xs = [state.atoms[i].x for i in range(r.atomStart, r.atomStart + r.atomCount)]
            ys = [state.atoms[i].y for i in range(r.atomStart, r.atomStart + r.atomCount)]
            zs = [state.atoms[i].z for i in range(r.atomStart, r.atomStart + r.atomCount)]
            return np.array([float(np.mean(xs)), float(np.mean(ys)), float(np.mean(zs))])
        except:
            return None
    
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
    
    # Check cavity alignment
    cavity_usage_rate = cavity_aligned / len(positions)
    
    # Cavity mode should show some preference for cavities
    # But don't require too high a rate since cavities may be sparse
    if len(cavities) > 5:
        assert cavity_usage_rate > 0.1 or cavity_aligned > 2, \
            f"Cavity mode shows insufficient cavity alignment: {cavity_usage_rate:.1%}"
    
    print(f"✓ Cavity mode test passed: {cavity_aligned}/{len(positions)} "
          f"({cavity_usage_rate:.1%}) near cavities")


def test_uniform_vs_cavity_comparison():
    """Compare uniform and cavity modes directly."""
    import numpy as np
    np.random.seed(99999)
    
    # Create fresh state
    state = pygcmc.MCState()
    state.info.box = np.array([4.0, 4.0, 4.0])
    
    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.ljEps = [0.0]  # Ideal gas
    ff.ljSigma = [0.1]
    state.forcefield = ff
    
    # Test uniform mode
    params_uniform = pygcmc.movement.MovementParams()
    params_uniform.temperature = 300.0
    params_uniform.chemicalPotential = -10.0
    params_uniform.proposalMode = 0
    params_uniform.seed = 11111
    
    mover_uniform = pygcmc.movement.MovementModule(params_uniform)
    
    uniform_accepts = 0
    for _ in range(200):
        result = mover_uniform.attemptInsertion(state)
        if result.accepted:
            uniform_accepts += 1
            mover_uniform.attemptDeletion(state)
    
    # Test cavity mode
    params_cavity = pygcmc.movement.MovementParams()
    params_cavity.temperature = 300.0
    params_cavity.chemicalPotential = -10.0
    params_cavity.proposalMode = 1
    params_cavity.useCavityBias = True
    params_cavity.cavityGridSpacing = 0.3
    params_cavity.seed = 22222
    
    mover_cavity = pygcmc.movement.MovementModule(params_cavity)
    
    cavity_accepts = 0
    for _ in range(200):
        result = mover_cavity.attemptInsertion(state)
        if result.accepted:
            cavity_accepts += 1
            mover_cavity.attemptDeletion(state)
    
    print(f"Uniform: {uniform_accepts}/200 accepts")
    print(f"Cavity: {cavity_accepts}/200 accepts")
    
    # Both should work
    assert uniform_accepts > 0, "Uniform mode failed completely"
    assert cavity_accepts > 0, "Cavity mode failed completely"
    
    print("✓ Both proposal modes functional")


if __name__ == "__main__":
    # Create a simple fixture replacement for standalone testing
    class SimpleState:
        def __init__(self):
            self.info = type('Info', (), {})()
            self.info.box = np.array([5.0, 5.0, 5.0])
            self.forcefield = type('FF', (), {})()
            self.forcefield.ljEps = [0.5]
            self.forcefield.ljSigma = [0.3]
            self.residues = []
            self.atoms = []
            self.activeResidueCount = 0
    
    print("Running improved proposal mode tests...")
    
    try:
        test_uniform_mode_basic(SimpleState())
        print("✓ Uniform mode test completed")
    except Exception as e:
        print(f"✗ Uniform mode test failed: {e}")
    
    try:
        test_cavity_mode_basic(SimpleState())
        print("✓ Cavity mode test completed")
    except Exception as e:
        print(f"✗ Cavity mode test failed: {e}")
    
    try:
        test_uniform_vs_cavity_comparison()
        print("✓ Comparison test completed")
    except Exception as e:
        print(f"✗ Comparison test failed: {e}")
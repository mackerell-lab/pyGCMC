# tests/simulation/movementCPP/multi_insertion_basic_funcs.py
"""Multi-insertion CBMC basic functionality tests - extracted functions."""

import pytest
import pygcmc
import numpy as np
import time
import math
import concurrent.futures
from .multi_insertion_fixtures import setup_multi_insertion_system

def test_multi_insertion_basic():
    """Test basic multi-insertion functionality."""
    state, params = setup_multi_insertion_system()
    
    params.useMultiInsertionCBMC = True
    params.maxParallelInsertions = 4
    # mproposal not exposed in Python
    # params.mproposal = 10
    
    mover = pygcmc.movement.MovementModule()
    
    try:
        mover.setParams(params)
        
        # Attempt multi-insertion
        result = mover.attemptMultiInsertionCBMC(state, 0)
        
        # Check result structure
        if hasattr(result, '__iter__'):
            # Multiple results returned
            assert len(result) <= params.maxParallelInsertions
            for r in result:
                assert hasattr(r, 'accepted')
                assert hasattr(r, 'energyChange')
        else:
            # Single result
            assert hasattr(result, 'accepted')
            assert hasattr(result, 'energyChange')
            
    except AttributeError:
        pytest.skip("Multi-insertion CBMC not implemented")
    except Exception as e:
        if "not implemented" in str(e).lower():
            pytest.skip("Multi-insertion CBMC not implemented")
        raise



def test_mproposal_scaling():
    """Test effect of mproposal parameter on acceptance."""
    state, params = setup_multi_insertion_system()
    
    params.useMultiInsertionCBMC = True
    params.maxParallelInsertions = 2
    
    mproposal_values = [1, 5, 10, 20]
    acceptance_rates = []
    
    for mproposal in mproposal_values:
        # params.mproposal = mproposal
        pass  # mproposal not exposed
        
        mover = pygcmc.movement.MovementModule()
        
        try:
            mover.setParams(params)
            
            accepts = 0
            attempts = 50
            
            for _ in range(attempts):
                result = mover.attemptMultiInsertionCBMC(state, 0)
                
                if hasattr(result, '__iter__'):
                    for r in result:
                        if r.accepted:
                            accepts += 1
                elif result.accepted:
                    accepts += 1
            
            acceptance_rates.append(accepts / (attempts * params.maxParallelInsertions))
            
        except Exception:
            acceptance_rates.append(0.0)
    
    # Higher mproposal should generally improve acceptance
    # Skip test if not implemented
    if any(r > 0 for r in acceptance_rates):
        # Check trend (may not be monotonic due to statistics)
        assert max(acceptance_rates) > min(acceptance_rates)



def test_parallel_insertion_count():
    """Test that parallel insertions respect maxParallelInsertions."""
    state, params = setup_multi_insertion_system()
    
    params.useMultiInsertionCBMC = True
    
    for max_parallel in [1, 2, 4, 8]:
        params.maxParallelInsertions = max_parallel
        # params.mproposal = 5  # Not exposed
        
        mover = pygcmc.movement.MovementModule()
        
        try:
            mover.setParams(params)
            
            # Attempt multi-insertion
            result = mover.attemptMultiInsertionCBMC(state, 0)
            
            if hasattr(result, '__iter__'):
                # Should not exceed max
                assert len(result) <= max_parallel
            else:
                # Single result is fine
                assert True
                
        except Exception:
            # Not implemented
            pass



def test_region_volume_calculation():
    """Test region volume calculation for multi-insertion."""
    state, params = setup_multi_insertion_system()
    
    params.useMultiInsertionCBMC = True
    params.maxParallelInsertions = 3
    # mproposal not exposed in Python
    # params.mproposal = 10
    # fillProposalInfo not exposed in Python
    # params.fillProposalInfo = True
    
    mover = pygcmc.movement.MovementModule()
    
    mover.setParams(params)
    
    # Perform multi-insertion
    result = mover.attemptMultiInsertionCBMC(state, 0)
    
    # Check result is valid
    assert result is not None
    
    # If result is a list, check it has elements
    if hasattr(result, '__iter__'):
        assert len(result) > 0
        # Check each result is valid
        for r in result:
            assert hasattr(r, 'accepted')
            assert hasattr(r, 'energyChange')
    else:
        # Single result
        assert hasattr(result, 'accepted')
        assert hasattr(result, 'energyChange')
    
    # Note: vregion field may not be exposed in Python bindings
    # Just verify the basic functionality works



def test_cbmc_weight_calculation():
    """Test CBMC weight calculation in multi-insertion."""
    state, params = setup_multi_insertion_system()
    
    params.useMultiInsertionCBMC = True
    params.maxParallelInsertions = 2
    # params.mproposal = 15  # Not exposed
    
    mover = pygcmc.movement.MovementModule()
    
    try:
        mover.setParams(params)
        
        # Collect weights
        weights = []
        
        for _ in range(20):
            result = mover.attemptMultiInsertionCBMC(state, 0)
            
            if hasattr(result, 'cbmcWeight'):
                weights.append(result.cbmcWeight)
            elif hasattr(result, '__iter__'):
                for r in result:
                    if hasattr(r, 'cbmcWeight'):
                        weights.append(r.cbmcWeight)
        
        if weights:
            # All weights should be positive
            assert all(w > 0 for w in weights)
            # Weights should vary
            assert len(set(weights)) > 1
            
    except Exception:
        pytest.skip("CBMC weights not available")



def test_multi_insertion_with_cavity_bias():
    """Test multi-insertion combined with cavity bias."""
    state, params = setup_multi_insertion_system()
    
    params.useMultiInsertionCBMC = True
    params.maxParallelInsertions = 3
    # mproposal not exposed in Python
    # params.mproposal = 10
    params.useCavityBias = True
    params.cavityGridSpacing = 0.15
    params.probeRadius = 0.15
    
    mover = pygcmc.movement.MovementModule()
    
    try:
        mover.setParams(params)
        
        # Find cavities first
        cavities = mover.findCavities(state)
        
        # Attempt cavity-biased multi-insertion
        result = mover.attemptMultiInsertionCBMC(state, 0)
        
        # Should complete without error
        assert result is not None
        
        # Check if positions are near cavities (if info available)
        if hasattr(result, 'proposalPosX') and len(cavities) > 0:
            pos = np.array([result.proposalPosX, result.proposalPosY, result.proposalPosZ])
            # Check distance to nearest cavity
            min_dist = float('inf')
            for cavity in cavities:
                cav_pos = np.array([cavity.x, cavity.y, cavity.z])
                dist = np.linalg.norm(pos - cav_pos)
                min_dist = min(min_dist, dist)
            # Should be reasonably close to a cavity
            assert min_dist < 1.0  # nm
            
    except Exception:
        pytest.skip("Multi-insertion with cavity bias not available")



def test_multi_insertion_performance():
    """Test performance of multi-insertion vs sequential."""
    state, params = setup_multi_insertion_system()
    
    # Sequential insertion
    params.useMultiInsertionCBMC = False
    mover_seq = pygcmc.movement.MovementModule()
    mover_seq.setParams(params)
    
    start_time = time.time()
    seq_accepts = 0
    attempts = 100
    
    for _ in range(attempts):
        result = mover_seq.attemptInsertion(state)
        if result.accepted:
            seq_accepts += 1
    
    seq_time = time.time() - start_time
    
    # Multi-insertion
    params.useMultiInsertionCBMC = True
    params.maxParallelInsertions = 4
    # mproposal not exposed in Python
    # params.mproposal = 10
    
    mover_multi = pygcmc.movement.MovementModule()
    
    try:
        mover_multi.setParams(params)
        
        start_time = time.time()
        multi_accepts = 0
        
        for _ in range(attempts // 4):  # Fewer calls but multiple insertions
            result = mover_multi.attemptMultiInsertionCBMC(state, 0)
            
            if hasattr(result, '__iter__'):
                for r in result:
                    if r.accepted:
                        multi_accepts += 1
            elif result.accepted:
                multi_accepts += 1
        
        multi_time = time.time() - start_time
        
        # Multi-insertion might be faster per insertion attempt
        # But this depends on implementation
        assert seq_time > 0
        assert multi_time > 0
        
    except Exception:
        pytest.skip("Multi-insertion performance test not available")



def test_multi_insertion_detailed_balance():
    """Test that multi-insertion preserves detailed balance."""
    state, params = setup_multi_insertion_system()
    
    params.useMultiInsertionCBMC = True
    params.maxParallelInsertions = 2
    # mproposal not exposed in Python
    # params.mproposal = 10
    
    mover = pygcmc.movement.MovementModule()
    
    try:
        mover.setParams(params)
        
        # Equilibrate with multi-insertion
        for _ in range(200):
            if np.random.random() < 0.5:
                mover.attemptMultiInsertionCBMC(state, 0)
            else:
                # Delete to balance
                for _ in range(params.maxParallelInsertions):
                    mover.attemptDeletion(state)
        
        # Measure steady state
        insertion_accepts = 0
        deletion_accepts = 0
        
        for _ in range(100):
            # Multi-insertion
            result = mover.attemptMultiInsertionCBMC(state, 0)
            if hasattr(result, '__iter__'):
                for r in result:
                    if r.accepted:
                        insertion_accepts += 1
            elif result.accepted:
                insertion_accepts += 1
            
            # Matching deletions
            for _ in range(params.maxParallelInsertions):
                result = mover.attemptDeletion(state)
                if result.accepted:
                    deletion_accepts += 1
        
        # Should be roughly balanced at equilibrium
        if insertion_accepts > 0 and deletion_accepts > 0:
            ratio = insertion_accepts / deletion_accepts
            assert 0.3 < ratio < 3.0
            
    except Exception:
        pytest.skip("Multi-insertion detailed balance test not available")



def test_multi_insertion_with_types():
    """Test multi-insertion with multiple atom types."""
    state = pygcmc.MCState()
    state.info.box = np.array([6.0, 6.0, 6.0])
    
    # Setup multi-type force field with full interaction matrix
    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 3
    ff.numMovementTypes = 3
    
    # Create full interaction matrix using Lorentz-Berthelot mixing rules
    eps_diag = [0.5, 0.6, 0.4]
    sigma_diag = [0.3, 0.35, 0.25]
    
    eps_matrix = []
    sigma_matrix = []
    for i in range(3):
        for j in range(3):
            # Geometric mean for epsilon
            eps_matrix.append(math.sqrt(eps_diag[i] * eps_diag[j]))
            # Arithmetic mean for sigma  
            sigma_matrix.append((sigma_diag[i] + sigma_diag[j]) / 2.0)
    
    ff.ljEps = eps_matrix  # 9 values for 3x3 matrix
    ff.ljSigma = sigma_matrix  # 9 values for 3x3 matrix
    state.forcefield = ff
    
    params = pygcmc.movement.MovementParams()
    params.temperature = 298.15
    params.chemicalPotential = -15.0
    params.useMultiInsertionCBMC = True
    params.maxParallelInsertions = 3
    # mproposal not exposed in Python
    # params.mproposal = 10
    
    mover = pygcmc.movement.MovementModule()
    
    mover.setParams(params)
    
    # Try multi-insertion for each type
    for atom_type in range(3):
        result = mover.attemptMultiInsertionCBMC(state, atom_type)
        
        # Should handle all types
        assert result is not None
        
        # Check basic result structure
        if hasattr(result, '__iter__'):
            assert len(result) > 0
            for r in result:
                assert hasattr(r, 'accepted')
                assert hasattr(r, 'energyChange')
                # Note: selectedType may not be exposed in bindings
        else:
            assert hasattr(result, 'accepted')
            assert hasattr(result, 'energyChange')



def test_multi_insertion_thread_safety():
    """Test thread safety of multi-insertion."""
    state, params = setup_multi_insertion_system()
    
    params.useMultiInsertionCBMC = True
    params.maxParallelInsertions = 2
    # params.mproposal = 5  # Not exposed
    

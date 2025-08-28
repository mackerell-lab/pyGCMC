# tests/simulation/movementCPP/proposal_reproducibility_funcs.py
"""Reproducibility test functions."""

import pytest
import pygcmc
import numpy as np


def test_seed_reproducibility():
    """Test that same seed produces same results."""
    seeds_to_test = [42, 12345, 99999]
    
    for seed in seeds_to_test:
        results1 = _run_with_seed(seed)
        results2 = _run_with_seed(seed)
        
        # Should produce identical results
        assert results1['positions'] == results2['positions'], \
            f"Different positions with seed {seed}"
        assert results1['accepts'] == results2['accepts'], \
            f"Different acceptance with seed {seed}"


def _run_with_seed(seed):
    """Helper to run insertion with specific seed."""
    state = pygcmc.MCState()
    state.info.box = np.array([5.0, 5.0, 5.0])
    
    # Need forcefield
    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.ljEps = [0.65]
    ff.ljSigma = [0.3]
    state.forcefield = ff
    
    params = pygcmc.movement.MovementParams()
    params.temperature = 298.15
    params.chemicalPotential = -15.0
    params.seed = seed
    params.fillProposalInfo = True
    params.useCavityBias = False  # Uniform for simplicity
    
    mover = pygcmc.movement.MovementModule()
    mover.setParams(params)
    
    positions = []
    accepts = 0
    
    for _ in range(20):
        result = mover.attemptInsertion(state)
        if result.proposalInfoFilled:
            positions.append((
                round(result.proposalPosX, 6),
                round(result.proposalPosY, 6),
                round(result.proposalPosZ, 6)
            ))
        if result.accepted:
            accepts += 1
    
    return {'positions': positions, 'accepts': accepts}


# Remove pytest.main call - not needed in module



"""
Improved Metropolis criterion validation test following the provided framework.
This version shows how to properly test p = min(1, exp(-βΔE)) for translation moves.
"""

import math
import numpy as np
import pytest
import pygcmc
from .test_helpers import save_forcefield, restore_forcefield


def create_two_particle_lj_system(box_len=8.0, sigma=1.0, epsilon=1.0, rcut=2.5):
    """Create a simple 2-particle LJ system for controlled testing.
    
    Args:
        box_len: Box length (should be >= 6*sigma to avoid PBC effects)
        sigma: LJ sigma parameter
        epsilon: LJ epsilon parameter
        rcut: Cutoff radius (typically 2.5*sigma)
        
    Returns:
        MCState with 2 LJ particles
    """
    state = pygcmc.MCState()
    state.info = pygcmc.SimInfo()
    state.info.box = (box_len, box_len, box_len)
    state.info.temperature = 1.0  # Reduced units
    
    # Create forcefield with LJ parameters
    ff = pygcmc.ForceField()
    # Note: This assumes the API supports setting LJ params
    # Actual implementation depends on PyGCMC API
    if hasattr(ff, 'ljEps'):
        ff.ljEps = [epsilon]
        ff.ljSigma = [sigma]
    else:
        # Alternative API
        ff.add_lj_params(0, epsilon, sigma)
    
    if hasattr(ff, 'masses'):
        ff.masses = [1.0]
    else:
        ff.add_atom_mass(0, 1.0)
    
    state.forcefield = ff
    
    # Add two particles at specific positions
    # Particle 1 at origin (will be fixed conceptually)
    atom1 = pygcmc.Atom()
    atom1.x, atom1.y, atom1.z = 0.0, 0.0, 0.0
    atom1.type = 0
    
    # Particle 2 at distance ~1.0*sigma (steep potential region)
    atom2 = pygcmc.Atom()
    atom2.x, atom2.y, atom2.z = 1.0 * sigma, 0.0, 0.0
    atom2.type = 0
    
    state.atoms = [atom1, atom2]
    
    # Create residues
    res1 = pygcmc.Residue()
    res1.atomStart = 0
    res1.atomCount = 1
    res1.active = True
    
    res2 = pygcmc.Residue()
    res2.atomStart = 1
    res2.atomCount = 1
    res2.active = True
    
    state.residues = [res1, res2]
    state.activeResidueCount = 2
    
    return state


def test_metropolis_criterion_validation_improved():
    """Test that translation moves follow Metropolis criterion p = min(1, exp(-βΔE)).
    
    This test:
    1. Creates a 2-particle LJ system
    2. Performs only translation moves
    3. Validates acceptance probabilities statistically
    """
    
    # Skip if deltaE is not available
    # First check if the field is exposed
    dummy_state = create_two_particle_lj_system()
    params = pygcmc.movement.MovementParams()
    params.temperature = 1.0
    params.seed = 12345
    mover = pygcmc.movement.MovementModule(params)
    
    # Test if deltaE is available
    test_result = mover.attemptTranslation(dummy_state)
    if not hasattr(test_result, 'deltaE'):
        pytest.skip("deltaE field not exposed in MovementResult - cannot validate Metropolis criterion")
    
    # Setup the actual test system
    state = create_two_particle_lj_system(
        box_len=8.0,  # Large enough to avoid PBC
        sigma=1.0,
        epsilon=1.0,
        rcut=2.5
    )
    
    # Configure movement parameters
    params = pygcmc.movement.MovementParams()
    params.temperature = 1.0  # Reduced units, so β = 1/T = 1.0
    params.seed = 12345
    params.proposalMode = 0  # Uniform mode (no cavity bias)
    params.useCavityBias = False  # Explicitly disable
    params.maxTranslation = 0.08  # Small moves for good sampling
    
    mover = pygcmc.movement.MovementModule(params)
    
    # Key parameters
    beta = 1.0 / params.temperature
    N_attempts = 10000
    
    # Accumulators for statistical test
    exp_sum = 0.0  # Sum of expected acceptance probabilities
    var_sum = 0.0  # Sum of variances
    acc_sum = 0    # Actual acceptances
    uphill_count = 0
    downhill_rejects = 0
    
    # Collect all deltaE values for analysis
    deltaE_values = []
    
    for i in range(N_attempts):
        # Only attempt translation of the second particle
        # (conceptually fixing the first)
        result = mover.attemptTranslation(state)
        
        if not hasattr(result, 'deltaE'):
            continue
            
        deltaE = result.deltaE
        deltaE_values.append(deltaE)
        
        if deltaE <= 0.0:
            # Downhill moves should always be accepted
            if not result.accepted:
                downhill_rejects += 1
        else:
            # Uphill moves: check Metropolis criterion
            uphill_count += 1
            p_metropolis = math.exp(-beta * deltaE)
            exp_sum += p_metropolis
            var_sum += p_metropolis * (1.0 - p_metropolis)
            
            if result.accepted:
                acc_sum += 1
    
    # Validation checks
    
    # 1. Must have sufficient uphill samples
    assert uphill_count > 100, \
        f"Insufficient uphill moves ({uphill_count}). Adjust max_translation or increase N_attempts"
    
    # 2. Downhill moves should never be rejected
    assert downhill_rejects == 0, \
        f"Found {downhill_rejects} rejected downhill moves - violates Metropolis criterion"
    
    # 3. Statistical test for uphill acceptance rate
    # Using z-score with threshold of 4.5 sigma
    if var_sum > 0:
        z_score = abs(acc_sum - exp_sum) / math.sqrt(var_sum)
        assert z_score < 4.5, \
            f"Metropolis validation failed: z-score = {z_score:.2f} " \
            f"(actual={acc_sum}, expected={exp_sum:.1f}, std={math.sqrt(var_sum):.1f})"
    
    # 4. Additional diagnostics (optional)
    print(f"Metropolis test passed:")
    print(f"  Total attempts: {N_attempts}")
    print(f"  Uphill moves: {uphill_count}")
    print(f"  Uphill acceptances: {acc_sum}")
    print(f"  Expected acceptances: {exp_sum:.1f}")
    print(f"  Z-score: {z_score:.2f}")
    
    # 5. Optional: Binned analysis for finer validation
    if len(deltaE_values) > 1000 and uphill_count > 200:
        validate_metropolis_by_bins(deltaE_values, beta, results=[])


def validate_metropolis_by_bins(deltaE_list, beta, results, n_bins=5):
    """Additional validation by binning deltaE values.
    
    Args:
        deltaE_list: List of all deltaE values
        beta: Inverse temperature
        results: List of MovementResult objects (if available)
        n_bins: Number of bins for analysis
    """
    # Filter for uphill moves
    uphill_deltaE = [dE for dE in deltaE_list if dE > 0]
    if len(uphill_deltaE) < n_bins * 20:
        return  # Not enough data for binned analysis
    
    # Sort and create equal-count bins
    uphill_deltaE.sort()
    bin_size = len(uphill_deltaE) // n_bins
    
    for i in range(n_bins):
        start_idx = i * bin_size
        end_idx = (i + 1) * bin_size if i < n_bins - 1 else len(uphill_deltaE)
        
        bin_deltaE = uphill_deltaE[start_idx:end_idx]
        bin_mean_dE = np.mean(bin_deltaE)
        
        # Expected acceptance for this bin
        exp_acc = sum(math.exp(-beta * dE) for dE in bin_deltaE)
        
        # Would need actual acceptance data per bin
        # This is a template for complete implementation
        print(f"  Bin {i+1}: ΔE range [{min(bin_deltaE):.3f}, {max(bin_deltaE):.3f}], "
              f"mean={bin_mean_dE:.3f}, expected_acc={exp_acc/len(bin_deltaE):.3f}")


def test_metropolis_with_external_potential():
    """Alternative test using external potential to create controlled energy landscape.
    
    This approach uses a harmonic potential to create predictable uphill moves.
    """
    # Skip if external potential or deltaE not available
    pytest.skip("External potential test - implement when API supports it")
    
    # Conceptual implementation:
    # 1. Create single particle in harmonic well
    # 2. Move particle away from minimum
    # 3. Validate acceptance follows exp(-β*k*Δx²/2)


# Minimal binding modification suggestion (for documentation)
BINDING_MODIFICATION_SUGGESTION = """
Minimal C++ binding modification to expose deltaE:

1. In MovementResult struct (C++):
   struct MovementResult {
       bool accepted;
       double acceptanceProbability;
       double deltaE;  // ADD THIS
       // ... other fields
   };

2. In movement attempt function (C++):
   MovementResult attemptTranslation(MCState& state) {
       // ... calculate deltaE ...
       double deltaE = calculateEnergyDifference(oldConfig, newConfig);
       
       // Metropolis criterion
       bool accepted = (deltaE <= 0) || (random() < exp(-beta * deltaE));
       
       MovementResult result;
       result.accepted = accepted;
       result.deltaE = deltaE;  // ADD THIS
       result.acceptanceProbability = min(1.0, exp(-beta * deltaE));
       
       return result;
   }

3. In pybind11 bindings:
   py::class_<MovementResult>(m, "MovementResult")
       .def_readonly("accepted", &MovementResult::accepted)
       .def_readonly("acceptanceProbability", &MovementResult::acceptanceProbability)
       .def_readonly("deltaE", &MovementResult::deltaE);  // ADD THIS

This minimal change enables complete Metropolis validation without affecting existing code.
"""
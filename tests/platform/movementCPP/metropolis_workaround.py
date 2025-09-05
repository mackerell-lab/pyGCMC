"""
Workaround for Metropolis criterion validation without deltaE field.
This shows how to validate using energy calculations as a temporary solution.
"""

import math
import numpy as np
import pytest
import pygcmc


def calculate_lj_energy(state, sigma=1.0, epsilon=1.0, rcut=2.5):
    """Calculate total LJ energy of the system.
    
    This is a workaround when the energy calculation API is available
    but deltaE is not directly exposed.
    
    Args:
        state: MCState object
        sigma: LJ sigma parameter
        epsilon: LJ epsilon parameter
        rcut: Cutoff radius
        
    Returns:
        Total energy
    """
    energy = 0.0
    atoms = state.atoms
    n_atoms = len(atoms)
    
    for i in range(n_atoms):
        for j in range(i + 1, n_atoms):
            # Calculate distance
            dx = atoms[j].x - atoms[i].x
            dy = atoms[j].y - atoms[i].y
            dz = atoms[j].z - atoms[i].z
            
            # Apply minimum image convention if needed
            Lx, Ly, Lz = state.info.box
            dx = dx - Lx * round(dx / Lx)
            dy = dy - Ly * round(dy / Ly)
            dz = dz - Lz * round(dz / Lz)
            
            r2 = dx*dx + dy*dy + dz*dz
            r = math.sqrt(r2)
            
            if r < rcut:
                # LJ potential
                r6i = (sigma / r) ** 6
                energy += 4.0 * epsilon * r6i * (r6i - 1.0)
    
    return energy


def test_metropolis_via_energy_calculation():
    """Test Metropolis criterion by calculating energies before and after moves.
    
    This is a workaround when deltaE is not exposed but we can:
    1. Save state before move
    2. Calculate energy before
    3. Attempt move
    4. Calculate energy after
    5. Compute deltaE = E_after - E_before
    """
    
    # Create simple 2-particle system
    state = pygcmc.MCState()
    state.info = pygcmc.SimInfo()
    state.info.box = (10.0, 10.0, 10.0)  # Large box to avoid PBC
    state.info.temperature = 1.0
    
    # Setup forcefield
    ff = pygcmc.ForceField()
    if hasattr(ff, 'ljEps'):
        ff.ljEps = [1.0]
        ff.ljSigma = [1.0]
    state.forcefield = ff
    
    # Add two particles
    atom1 = pygcmc.Atom()
    atom1.x, atom1.y, atom1.z = 0.0, 0.0, 0.0
    atom1.type = 0
    
    atom2 = pygcmc.Atom()
    atom2.x, atom2.y, atom2.z = 1.2, 0.0, 0.0  # Near minimum
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
    
    # Setup movement
    params = pygcmc.movement.MovementParams()
    params.temperature = 1.0
    params.seed = 42
    params.proposalMode = 0  # Uniform
    params.useCavityBias = False
    params.maxTranslation = 0.05  # Small moves
    
    mover = pygcmc.movement.MovementModule(params)
    
    # Test parameters
    beta = 1.0 / params.temperature
    n_trials = 1000
    
    # Track statistics
    uphill_data = []
    downhill_rejected = 0
    
    for _ in range(n_trials):
        # Save current positions
        old_positions = [(a.x, a.y, a.z) for a in state.atoms]
        
        # Calculate energy before
        E_before = calculate_lj_energy(state)
        
        # Attempt translation
        result = mover.attemptTranslation(state)
        
        # Calculate energy after
        E_after = calculate_lj_energy(state)
        
        # Compute deltaE
        deltaE = E_after - E_before
        
        if result.accepted:
            # Move was accepted, energy changed
            if deltaE > 0:
                # Uphill move accepted
                uphill_data.append({
                    'deltaE': deltaE,
                    'accepted': True,
                    'prob': math.exp(-beta * deltaE)
                })
        else:
            # Move was rejected, restore positions
            for i, (x, y, z) in enumerate(old_positions):
                state.atoms[i].x = x
                state.atoms[i].y = y
                state.atoms[i].z = z
            
            if deltaE <= 0:
                downhill_rejected += 1
            else:
                uphill_data.append({
                    'deltaE': deltaE,
                    'accepted': False,
                    'prob': math.exp(-beta * deltaE)
                })
    
    # Basic validation
    assert downhill_rejected == 0, f"Rejected {downhill_rejected} downhill moves!"
    
    # Statistical validation for uphill moves
    if len(uphill_data) > 50:
        # Group by deltaE bins
        deltaE_values = [d['deltaE'] for d in uphill_data]
        min_dE = min(deltaE_values)
        max_dE = max(deltaE_values)
        
        # Create bins
        n_bins = min(5, len(uphill_data) // 20)
        bins = np.linspace(min_dE, max_dE, n_bins + 1)
        
        for i in range(n_bins):
            bin_data = [d for d in uphill_data 
                       if bins[i] <= d['deltaE'] < bins[i+1]]
            
            if len(bin_data) > 10:
                expected_acc = sum(d['prob'] for d in bin_data)
                actual_acc = sum(d['accepted'] for d in bin_data)
                
                # Binomial test
                variance = sum(d['prob'] * (1 - d['prob']) for d in bin_data)
                if variance > 0:
                    z = abs(actual_acc - expected_acc) / math.sqrt(variance)
                    
                    # More lenient threshold for workaround
                    assert z < 5.0, \
                        f"Bin {i}: z={z:.2f} exceeds threshold"
    
    print(f"Workaround test completed: {len(uphill_data)} uphill moves analyzed")


def test_metropolis_dry_run_approach():
    """Alternative approach using dry-run if API supports it.
    
    Some implementations allow a 'dry run' mode where you can:
    1. Propose a move without applying it
    2. Get the proposed configuration
    3. Calculate energies for both configurations
    4. Compare with actual acceptance
    """
    
    # This would require API support for:
    # - mover.proposeTranslation(state, dry_run=True) -> proposed_config
    # - energy.calculate(config) -> energy
    
    pytest.skip("Dry-run approach requires API support for proposed configurations")


# Additional helper for debugging
def analyze_metropolis_deviation(results_with_deltaE):
    """Analyze where Metropolis criterion might be violated.
    
    Args:
        results_with_deltaE: List of dicts with 'deltaE', 'accepted', 'prob_reported'
        
    Returns:
        Dictionary with analysis results
    """
    analysis = {
        'total': len(results_with_deltaE),
        'downhill_rejected': 0,
        'uphill_stats': {},
        'deviation_by_range': []
    }
    
    # Separate by deltaE sign
    downhill = [r for r in results_with_deltaE if r['deltaE'] <= 0]
    uphill = [r for r in results_with_deltaE if r['deltaE'] > 0]
    
    # Check downhill
    analysis['downhill_rejected'] = sum(not r['accepted'] for r in downhill)
    
    # Analyze uphill by deltaE magnitude
    if uphill:
        # Sort by deltaE
        uphill.sort(key=lambda x: x['deltaE'])
        
        # Create quartiles
        n = len(uphill)
        quartiles = [
            uphill[:n//4],
            uphill[n//4:n//2],
            uphill[n//2:3*n//4],
            uphill[3*n//4:]
        ]
        
        for i, q_data in enumerate(quartiles):
            if q_data:
                mean_dE = np.mean([r['deltaE'] for r in q_data])
                expected_acc_rate = np.mean([math.exp(-r['deltaE']) for r in q_data])
                actual_acc_rate = np.mean([r['accepted'] for r in q_data])
                
                analysis['deviation_by_range'].append({
                    'quartile': i + 1,
                    'mean_deltaE': mean_dE,
                    'expected_acc': expected_acc_rate,
                    'actual_acc': actual_acc_rate,
                    'ratio': actual_acc_rate / (expected_acc_rate + 1e-10)
                })
    
    return analysis


# Export the binding modification suggestion
def get_minimal_binding_change():
    """Return the minimal C++ binding change needed to expose deltaE."""
    return """
    // In C++ MovementResult struct:
    struct MovementResult {
        bool accepted;
        double acceptanceProbability;
        double deltaE;  // <-- ADD THIS FIELD
        int residueIndex;
        // ... other fields ...
    };
    
    // In movement implementation:
    MovementResult MovementModule::attemptTranslation(MCState& state) {
        // ... select particle, propose new position ...
        
        double deltaE = calculateEnergyDifference(oldPos, newPos);
        
        bool accepted = acceptMove(deltaE, temperature);
        
        MovementResult result;
        result.accepted = accepted;
        result.deltaE = deltaE;  // <-- SET THE VALUE HERE
        result.acceptanceProbability = std::min(1.0, std::exp(-deltaE/kT));
        
        if (accepted) {
            applyMove(state, newPos);
        }
        
        return result;
    }
    
    // In Python bindings (pybind11):
    py::class_<MovementResult>(m, "MovementResult")
        .def_readonly("accepted", &MovementResult::accepted)
        .def_readonly("acceptanceProbability", &MovementResult::acceptanceProbability)
        .def_readonly("deltaE", &MovementResult::deltaE)  // <-- EXPOSE TO PYTHON
        .def_readonly("residueIndex", &MovementResult::residueIndex);
    
    This minimal change (3 lines total) enables complete Metropolis validation.
    """
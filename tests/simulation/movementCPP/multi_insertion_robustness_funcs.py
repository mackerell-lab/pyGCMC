# tests/simulation/movementCPP/multi_insertion_robustness_funcs.py
"""Multi-insertion CBMC robustness tests - extracted functions."""

import pytest
import pygcmc
import numpy as np
import os
from .multi_insertion_robustness_fixtures import basic_robustness_system, make_test_state

def test_acceptance_rate_monotonicity_chemical_potential():
    """Test acceptance rate increases with chemical potential."""
    state, params = basic_robustness_system()
    
    # Test range of chemical potentials
    mu_values = [-20.0, -15.0, -10.0, -5.0, 0.0]
    acceptance_rates = []
    
    for mu in mu_values:
        params.chemicalPotential = mu
        mover = pygcmc.movement.MovementModule()
        mover.setParams(params)
        
        # Run trials
        n_trials = 100
        accepts = 0
        attempts = 0
        
        for _ in range(n_trials):
            # Fresh state each time
            test_state = make_test_state()
            
            result = mover.attemptMultiInsertionCBMC(test_state, 0)
            if hasattr(result, '__iter__'):
                attempts += len(result)
                for r in result:
                    if r.accepted:
                        accepts += 1
            else:
                attempts += 1
                if result.accepted:
                    accepts += 1
        
        rate = accepts / attempts if attempts > 0 else 0
        acceptance_rates.append(rate)
    
    # Check monotonicity (allowing small numerical noise)
    for i in range(1, len(acceptance_rates)):
        assert acceptance_rates[i] >= acceptance_rates[i-1] - 0.01, \
            f"Acceptance rate decreases: μ={mu_values[i-1]}→{mu_values[i]}, " \
            f"rate={acceptance_rates[i-1]:.3f}→{acceptance_rates[i]:.3f}"
    

def test_acceptance_rate_monotonicity_volume():
    """Test acceptance rate increases with box volume."""
    state, params = basic_robustness_system()
    params.chemicalPotential = -10.0
    
    # Test range of box sizes (keeping shape cubic)
    box_sizes = [2.5, 3.0, 3.5, 4.0]  # nm
    acceptance_rates = []
    
    for box_size in box_sizes:
        test_state = make_test_state(box_size, min(1.2, box_size * 0.4))
        
        mover = pygcmc.movement.MovementModule()
        mover.setParams(params)
        
        # Run trials
        n_trials = 100
        accepts = 0
        attempts = 0
        
        for _ in range(n_trials):
            # Fresh state
            local_state = make_test_state(box_size, min(1.2, box_size * 0.4))
            
            result = mover.attemptMultiInsertionCBMC(local_state, 0)
            if hasattr(result, '__iter__'):
                attempts += len(result)
                for r in result:
                    if r.accepted:
                        accepts += 1
            else:
                attempts += 1
                if result.accepted:
                    accepts += 1
        
        rate = accepts / attempts if attempts > 0 else 0
        acceptance_rates.append(rate)
    
    # Check monotonicity
    for i in range(1, len(acceptance_rates)):
        assert acceptance_rates[i] >= acceptance_rates[i-1] - 0.01, \
            f"Acceptance rate decreases: V={box_sizes[i-1]**3:.1f}→{box_sizes[i]**3:.1f} nm³, " \
            f"rate={acceptance_rates[i-1]:.3f}→{acceptance_rates[i]:.3f}"
    

def test_energy_consistency_single_vs_system():
    """Test energy calculation consistency between methods."""
    state, params = basic_robustness_system()
    params.chemicalPotential = -5.0
    
    # First add some molecules to have non-zero interactions
    mover = pygcmc.movement.MovementModule()
    mover.setParams(params)
    
    # Add a few molecules
    for _ in range(3):
        result = mover.attemptInsertion(state, 0)
        if not result.accepted:
            # Force add for testing
            atom1 = pygcmc.MCAtom()
            atom1.x, atom1.y, atom1.z = 1.0, 1.0, 1.0
            atom1.type = 0
            atom1.charge = -0.834
            
            atom2 = pygcmc.MCAtom()
            atom2.x, atom2.y, atom2.z = 1.1, 1.0, 1.0
            atom2.type = 1
            atom2.charge = 0.417
            
            atom3 = pygcmc.MCAtom()
            atom3.x, atom3.y, atom3.z = 1.0, 1.1, 1.0
            atom3.type = 1
            atom3.charge = 0.417
            
            state.addAtom(atom1)
            state.addAtom(atom2)
            state.addAtom(atom3)
            
            residue = pygcmc.MCResidue()
            residue.atomStart = state.activeAtomCount - 3
            residue.atomCount = 3
            residue.type = 0
            residue.active = True
            state.addResidue(residue)
    
    if state.activeResidueCount == 0:
        pytest.skip("Could not create non-empty system")
    
    # Calculate energy before insertion
    # Use the exposed C++ function directly
    pygcmc.platform.cpu.computeSystemEnergyPBCCutoff(state)
    energy_before = sum(state.residues[i].energy_vdw + state.residues[i].energy_elec 
                       for i in range(state.activeResidueCount)) * 0.5
    
    # Add a test molecule
    test_atoms = [
        (2.0, 2.0, 2.0, 0, -0.834),  # O
        (2.1, 2.0, 2.0, 1, 0.417),   # H1
        (2.0, 2.1, 2.0, 1, 0.417)    # H2
    ]
    
    start_idx = state.activeAtomCount
    for x, y, z, typ, charge in test_atoms:
        atom = pygcmc.MCAtom()
        atom.x, atom.y, atom.z = x, y, z
        atom.type = typ
        atom.charge = charge
        state.addAtom(atom)
    
    res = pygcmc.MCResidue()
    res.atomStart = start_idx
    res.atomCount = 3
    res.type = 0
    res.active = True
    res_idx = state.addResidue(res)
    
    # Method 1: Single residue energy
    pygcmc.platform.cpu.computeResidueEnergyCutoffPBC(state, res_idx)
    single_residue_energy = state.residues[res_idx].energy_vdw + \
                            state.residues[res_idx].energy_elec
    
    # Method 2: System energy difference
    pygcmc.platform.cpu.computeSystemEnergyPBCCutoff(state)
    energy_after = sum(state.residues[i].energy_vdw + state.residues[i].energy_elec 
                      for i in range(state.activeResidueCount)) * 0.5
    
    system_energy_diff = energy_after - energy_before
    
    # Remove test molecule
    state.removeResidue(res_idx)
    for _ in range(3):
        state.removeAtom(state.activeAtomCount - 1)
    
    # Check consistency (allowing numerical tolerance)
    assert abs(system_energy_diff - single_residue_energy) < 1e-6, \
        f"Energy inconsistency: system diff={system_energy_diff:.8f}, " \
        f"single residue={single_residue_energy:.8f}"
    

def test_region_independence_boundary():
    """Test region selection with extreme separations."""
    state, params = basic_robustness_system()
    
    # Case 1: Separation near box size
    params.minRegionSeparationNm = 2.9  # Near 3.0 nm box
    mover = pygcmc.movement.MovementModule()
    mover.setParams(params)
    
    # Should still work but select fewer regions
    result = mover.attemptMultiInsertionCBMC(state, 0)
    if hasattr(result, '__iter__'):
        # With 5nm box and 2.9nm separation, at most 2 regions can fit
        assert len(result) <= 2, \
            f"Too many regions ({len(result)}) selected with large separation"
    
    # Case 2: Very small cutoff (force Keff=0 scenario)
    state.info.cutoff = 0.01  # nm - extremely small
    state.forcefield.ljSigma = [10.0, 10.0, 10.0, 10.0]  # Huge sigma
    
    params.minRegionSeparationNm = 1.0
    mover.setParams(params)
    
    # Should reject all but not crash
    result = mover.attemptMultiInsertionCBMC(state, 0)
    if hasattr(result, '__iter__'):
        for r in result:
            assert not r.accepted, "Should reject with extreme parameters"
    else:
        assert not result.accepted, "Should reject with extreme parameters"
    

def test_deterministic_config_selection():
    """Test reproducibility of configuration selection."""
    state, params = basic_robustness_system()
    params.seed = 12345
    
    # Run twice with same seed
    results1 = []
    results2 = []
    
    for run in [results1, results2]:
        mover = pygcmc.movement.MovementModule()
        mover.setParams(params)
        
        # Fresh state
        test_state = make_test_state()
        
        result = mover.attemptMultiInsertionCBMC(test_state, 0)
        
        if hasattr(result, '__iter__'):
            for r in result:
                run.append({
                    'accepted': r.accepted,
                    'energy': r.energyChange if hasattr(r, 'energyChange') else 0
                })
        else:
            run.append({
                'accepted': result.accepted,
                'energy': result.energyChange if hasattr(result, 'energyChange') else 0
            })
    
    # Check reproducibility
    assert len(results1) == len(results2), \
        f"Different number of results: {len(results1)} vs {len(results2)}"
    
    for i, (r1, r2) in enumerate(zip(results1, results2)):
        assert r1['accepted'] == r2['accepted'], \
            f"Region {i}: acceptance differs"
        assert abs(r1['energy'] - r2['energy']) < 1e-10, \
            f"Region {i}: energy differs {r1['energy']} vs {r2['energy']}"
    

def test_fallback_path_unit_consistency():
    """Test unit consistency in fallback creation path."""
    state, params = basic_robustness_system()
    
    # Don't set params.useMultiInsertionCBMC explicitly
    # This forces the fallback path in attemptMultiInsertionCBMC
    params_fresh = pygcmc.movement.MovementParams()
    params_fresh.temperature = params.temperature
    params_fresh.beta = params.beta
    params_fresh.maxParallelInsertions = 2
    params_fresh.numConfigTrials = 10
    params_fresh.minRegionSeparationNm = 1.5  # nm
    params_fresh.chemicalPotential = -10.0
    params_fresh.seed = 42
    
    # Create mover without setting useMultiInsertionCBMC
    mover = pygcmc.movement.MovementModule()
    # Set params but without useMultiInsertionCBMC to test fallback
    mover.setParams(params_fresh)
    
    # This should trigger the fallback initialization
    result = mover.attemptMultiInsertionCBMC(state, 0)
    
    # Should work without unit conversion errors
    assert result is not None, "Fallback path failed"
    
    # Compare with properly initialized path
    mover2 = pygcmc.movement.MovementModule()
    params_fresh.useMultiInsertionCBMC = True
    mover2.setParams(params_fresh)
    
    # Run multiple trials to check statistical consistency
    n_trials = 50
    accepts1 = 0
    accepts2 = 0
    
    for _ in range(n_trials):
        # Fallback path
        test_state1 = make_test_state()
        
        result1 = mover.attemptMultiInsertionCBMC(test_state1, 0)
        if hasattr(result1, '__iter__'):
            for r in result1:
                if r.accepted:
                    accepts1 += 1
        elif result1.accepted:
            accepts1 += 1
        
        # Normal path
        test_state2 = make_test_state()
        
        result2 = mover2.attemptMultiInsertionCBMC(test_state2, 0)
        if hasattr(result2, '__iter__'):
            for r in result2:
                if r.accepted:
                    accepts2 += 1
        elif result2.accepted:
            accepts2 += 1
    
    # Acceptance rates should be statistically similar
    rate1 = accepts1 / (n_trials * params_fresh.maxParallelInsertions)
    rate2 = accepts2 / (n_trials * params_fresh.maxParallelInsertions)
    
    # Allow some statistical variation
    assert abs(rate1 - rate2) < 0.2, \
        f"Unit inconsistency detected: fallback rate={rate1:.3f}, " \
        f"normal rate={rate2:.3f}"
    

def test_energy_range_parameterization():
    """Test energy values stay within expected ranges."""
    state, params = basic_robustness_system()
    params.chemicalPotential = -10.0
    
    # Expected range for water at 298K (empirical)
    ENERGY_MIN = -1000.0  # kJ/mol
    ENERGY_MAX = 100.0    # kJ/mol
    
    mover = pygcmc.movement.MovementModule()
    mover.setParams(params)
    
    # Add some molecules first
    for _ in range(5):
        mover.attemptInsertion(state, 0)
    
    # Collect energy changes
    energies = []
    for _ in range(50):
        result = mover.attemptMultiInsertionCBMC(state, 0)
        if hasattr(result, '__iter__'):
            for r in result:
                if hasattr(r, 'energyChange'):
                    energies.append(r.energyChange)
        elif hasattr(result, 'energyChange'):
            energies.append(result.energyChange)
    
    # Check range
    if energies:
        assert all(ENERGY_MIN <= e <= ENERGY_MAX for e in energies), \
            f"Energy outside expected range [{ENERGY_MIN}, {ENERGY_MAX}]: " \
            f"min={min(energies):.2f}, max={max(energies):.2f}"
@pytest.mark.skipif(not hasattr(pygcmc, 'PYGCMC_USE_OPENMP'), 
                    reason="OpenMP not enabled")
def test_openmp_determinism():
    """Test determinism with OpenMP if enabled."""
    import os
    state, params = basic_robustness_system()
    params.seed = 99999
    
    # Run with single thread
    os.environ['OMP_NUM_THREADS'] = '1'
    mover1 = pygcmc.movement.MovementModule()
    mover1.setParams(params)
    
    result1 = mover1.attemptMultiInsertionCBMC(state, 0)
    
    # Run again with single thread  
    state2 = make_test_state()
    
    mover2 = pygcmc.movement.MovementModule()
    mover2.setParams(params)
    
    result2 = mover2.attemptMultiInsertionCBMC(state2, 0)
    
    # Should be identical
    if hasattr(result1, '__iter__') and hasattr(result2, '__iter__'):
        assert len(result1) == len(result2)
        for r1, r2 in zip(result1, result2):
            assert r1.accepted == r2.accepted
            if hasattr(r1, 'energyChange') and hasattr(r2, 'energyChange'):
                assert abs(r1.energyChange - r2.energyChange) < 1e-10
    else:
        assert result1.accepted == result2.accepted
        if hasattr(result1, 'energyChange') and hasattr(result2, 'energyChange'):
            assert abs(result1.energyChange - result2.energyChange) < 1e-10
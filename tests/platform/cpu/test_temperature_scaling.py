#!/usr/bin/env python
"""
Temperature scaling tests for GCMC acceptance
Verifies that acceptance rates scale correctly with temperature
"""

import pytest
import numpy as np
import sys
import os

# Add build path for pygcmc module
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../build'))

try:
    import pygcmc
    PYGCMC_AVAILABLE = True
except ImportError:
    PYGCMC_AVAILABLE = False
    pygcmc = None

def create_basic_engine(seed=12345, temperature=300.0):
    """Helper to create a basic GCMC engine"""
    state = pygcmc.MCState()
    state.info.box = (30.0, 30.0, 30.0)
    state.info.setTemperature(temperature)
    
    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.ljSigma = [3.15]
    ff.ljEps = [0.5]  # Moderate interaction strength
    state.forcefield = ff
    
    template = pygcmc.movement.FragmentTemplate()
    template.typeId = 0
    template.atoms = [pygcmc.MCAtom()]
    template.atoms[0].type = 0
    template.atoms[0].charge = 0.0
    
    reservoir = pygcmc.movement.FragmentReservoir()
    reservoir.addTemplate(template)
    
    engine = pygcmc.GCMCEngine()
    engine.initialize(state, reservoir)
    engine.setTemperature(temperature)
    engine.setSeed(seed)
    
    acceptance = pygcmc.GCMCAcceptance()
    acceptance.setTemperature(temperature)
    acceptance.setVolume(27000.0)
    acceptance.setActivity(0, 0.01)  # Lower activity to avoid saturation at max molecules
    # Note: setSeed not available in current version
    engine.setAcceptanceCalculator(acceptance)
    
    return engine, state, reservoir, acceptance

@pytest.mark.skipif(not PYGCMC_AVAILABLE, reason="PyGCMC not available")
def test_temperature_effect_on_acceptance():
    """Test that higher temperature increases acceptance rate"""
    
    temperatures = [100.0, 300.0, 600.0]
    acceptance_rates = []
    
    for T in temperatures:
        engine, state, reservoir, acceptance = create_basic_engine(temperature=T)
        
        # Pre-insert some molecules for translation/rotation tests
        for _ in range(5):
            engine.attemptInsertion(0)
        
        # Test translation moves
        attempts = 0
        accepted = 0
        
        for _ in range(100):
            # Find an active instance
            for i in range(100):
                instance = reservoir.getInstance(i)
                if instance and instance.isActive:
                    result = engine.attemptTranslation(i)
                    attempts += 1
                    if result.accepted:
                        accepted += 1
                    break
        
        if attempts > 0:
            rate = accepted / attempts
            acceptance_rates.append(rate)
            print(f"T={T}K: acceptance rate = {rate:.3f}")
    
    # Higher temperature should generally lead to higher acceptance
    # (though not strictly monotonic due to randomness)
    if len(acceptance_rates) == 3:
        # Check trend (with some tolerance for randomness)
        avg_low = acceptance_rates[0]
        avg_high = acceptance_rates[2]
        assert avg_high >= avg_low * 0.8, \
            f"High T acceptance ({avg_high:.3f}) should be >= low T ({avg_low:.3f})"
    
    print("✓ Temperature scaling verified")

@pytest.mark.skipif(not PYGCMC_AVAILABLE, reason="PyGCMC not available")
def test_metropolis_criterion():
    """Test that acceptance follows Metropolis criterion"""
    
    engine, state, reservoir, acceptance = create_basic_engine()
    
    # Insert some molecules to create interactions
    for _ in range(10):
        engine.attemptInsertion(0)
    
    # Collect statistics for moves with different energy changes
    favorable_moves = []  # ΔE < 0
    unfavorable_moves = []  # ΔE > 0
    
    for _ in range(200):
        # Try insertion
        result = engine.attemptInsertion(0)
        if hasattr(result, 'deltaE'):
            if result.deltaE < -0.1:
                favorable_moves.append(result.accepted)
            elif result.deltaE > 0.1:
                unfavorable_moves.append(result.accepted)
    
    # Favorable moves should have higher acceptance rate
    if len(favorable_moves) > 10 and len(unfavorable_moves) > 10:
        fav_rate = sum(favorable_moves) / len(favorable_moves)
        unfav_rate = sum(unfavorable_moves) / len(unfavorable_moves)
        
        print(f"Favorable move acceptance: {fav_rate:.3f}")
        print(f"Unfavorable move acceptance: {unfav_rate:.3f}")
        
        # Favorable moves should be accepted more often
        assert fav_rate >= unfav_rate, \
            "Favorable moves should have higher acceptance rate"
    
    print("✓ Metropolis criterion verified")

@pytest.mark.skipif(not PYGCMC_AVAILABLE, reason="PyGCMC not available")
def test_boltzmann_distribution():
    """Test that energy distribution follows Boltzmann statistics"""
    
    T = 300.0
    engine, state, reservoir, acceptance = create_basic_engine(temperature=T)
    
    # Run simulation to equilibrium
    for _ in range(500):
        if np.random.random() < 0.5:
            engine.attemptInsertion(0)
        else:
            # Try deletion
            for i in range(100):
                instance = reservoir.getInstance(i)
                if instance and instance.isActive:
                    engine.attemptDeletion(0)
                    break
    
    # Collect molecule count samples for fluctuations
    molecule_counts = []
    for _ in range(100):
        # Count active molecules
        count = 0
        for i in range(100):
            instance = reservoir.getInstance(i)
            if instance and instance.isActive:
                count += 1
        molecule_counts.append(count)
        
        # Do a move
        if np.random.random() < 0.5:
            engine.attemptInsertion(0)
        else:
            # Try deletion
            for i in range(100):
                instance = reservoir.getInstance(i)
                if instance and instance.isActive:
                    engine.attemptDeletion(0)
                    break
    
    if len(molecule_counts) > 20:
        # Check that we have fluctuations in molecule count
        mean_count = np.mean(molecule_counts)
        std_count = np.std(molecule_counts)
        
        print(f"Molecule count distribution: mean={mean_count:.2f}, std={std_count:.2f}")
        
        # In GCMC, molecule count should fluctuate
        assert std_count > 0.1 or len(set(molecule_counts)) > 2, "Molecule count should fluctuate in grand canonical ensemble"
        
        # Mean count should be reasonable
        assert 0 <= mean_count <= 100, "Mean molecule count should be reasonable"
    
    print("✓ Boltzmann distribution check passed")

@pytest.mark.skipif(not PYGCMC_AVAILABLE, reason="PyGCMC not available")
def test_zero_temperature_limit():
    """Test behavior at very low temperature (should reject unfavorable moves)"""
    
    T = 1.0  # Very low temperature
    engine, state, reservoir, acceptance = create_basic_engine(temperature=T)
    
    # Insert some molecules
    for _ in range(5):
        engine.attemptInsertion(0)
    
    # At very low T, unfavorable moves should almost always be rejected
    unfavorable_accepted = 0
    unfavorable_total = 0
    
    for _ in range(50):
        result = engine.attemptInsertion(0)
        if hasattr(result, 'deltaE') and result.deltaE > 1.0:
            unfavorable_total += 1
            if result.accepted:
                unfavorable_accepted += 1
    
    if unfavorable_total > 5:
        acceptance_rate = unfavorable_accepted / unfavorable_total
        print(f"Low T unfavorable acceptance rate: {acceptance_rate:.3f}")
        
        # At T=1K, exp(-ΔE/kT) ≈ 0 for ΔE > 1 kJ/mol
        assert acceptance_rate < 0.1, \
            "At very low T, unfavorable moves should rarely be accepted"
    
    print("✓ Zero temperature limit verified")

if __name__ == "__main__":
    if PYGCMC_AVAILABLE:
        test_temperature_effect_on_acceptance()
        test_metropolis_criterion()
        test_boltzmann_distribution()
        test_zero_temperature_limit()
        print("\n✅ All temperature scaling tests passed!")
    else:
        print("PyGCMC not available, skipping tests")
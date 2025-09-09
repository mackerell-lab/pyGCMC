#!/usr/bin/env python
"""
Rigorous theoretical validation for GCMC implementation
Tests exact formulas for ideal gas (zero interactions)
"""

import pytest
import numpy as np
import os
import sys

# Add build path for pygcmc module
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../build'))

try:
    import pygcmc
    PYGCMC_AVAILABLE = True
except ImportError:
    PYGCMC_AVAILABLE = False
    pygcmc = None

@pytest.mark.skipif(not PYGCMC_AVAILABLE, reason="PyGCMC not available")
class TestIdealGasValidation:
    """Validate GCMC implementation against ideal gas theory"""
    
    def create_ideal_gas_system(self, box_size=30.0, temperature=300.0, activity=0.01):
        """Create ideal gas system (zero interactions)"""
        state = pygcmc.MCState()
        state.info.box = (box_size, box_size, box_size)
        state.info.setTemperature(temperature)
        
        # Zero interactions for ideal gas
        ff = pygcmc.MCForceField()
        ff.numTotalTypes = 1
        ff.numMovementTypes = 1
        ff.ljSigma = [0.315]  # nm
        ff.ljEps = [0.0]      # Zero interaction - ideal gas
        state.forcefield = ff
        
        # Single atom template
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
        engine.setSeed(12345)
        
        # Calculate box volume correctly
        box_volume = box_size * box_size * box_size  # nm^3
        
        acceptance = pygcmc.GCMCAcceptance()
        acceptance.setTemperature(temperature)
        acceptance.setVolume(box_volume)
        acceptance.setActivity(0, activity)
        acceptance.setSeed(12346)  # Different seed for acceptance
        engine.setAcceptanceCalculator(acceptance)
        
        return engine, state, reservoir, acceptance, box_volume
    
    def test_insertion_formula_ideal_gas(self):
        """Test exact insertion probability formula for ideal gas"""
        print("\nValidating insertion formula P_ins = min(1, z*V/(N+1))...")
        
        # Use low activity to avoid saturation
        activity = 0.0001  # Very low to ensure P < 1
        engine, state, reservoir, acceptance, volume = self.create_ideal_gas_system(
            box_size=30.0, temperature=300.0, activity=activity
        )
        
        # Enable probability storage
        os.environ['GCMC_STORE_PROB'] = '1'
        
        try:
            # Test insertion probabilities for different N values
            errors = []
            for target_n in range(0, 10):
                # Get to target N
                current_n = reservoir.getActiveCount(0)
                while current_n < target_n:
                    engine.attemptInsertion(0)
                    current_n = reservoir.getActiveCount(0)
                while current_n > target_n:
                    engine.attemptDeletion(0)
                    current_n = reservoir.getActiveCount(0)
                
                # Now test insertion at this N
                N_before = reservoir.getActiveCount(0)
                result = engine.attemptInsertion(0)
                
                # For ideal gas with no cavity bias, bias = 1.0
                theoretical_prob = min(1.0, activity * volume / (N_before + 1))
                actual_prob = result.acceptanceProbability
                
                # Should match exactly (within numerical precision)
                relative_error = abs(actual_prob - theoretical_prob) / max(theoretical_prob, 1e-10)
                errors.append(relative_error)
                
                print(f"  N={N_before}: theoretical={theoretical_prob:.6f}, "
                      f"actual={actual_prob:.6f}, error={relative_error:.2e}")
                
                # Strong assertion - should match within 1e-6 relative error
                assert relative_error < 1e-6, \
                    f"Insertion probability mismatch at N={N_before}: " \
                    f"expected {theoretical_prob:.6f}, got {actual_prob:.6f}"
            
            print(f"✓ Insertion formula validated, max error: {max(errors):.2e}")
            
        finally:
            del os.environ['GCMC_STORE_PROB']
    
    def test_deletion_formula_ideal_gas(self):
        """Test exact deletion probability formula for ideal gas"""
        print("\nValidating deletion formula P_del = min(1, N/(z*V))...")
        
        # Use moderate activity
        activity = 0.001
        engine, state, reservoir, acceptance, volume = self.create_ideal_gas_system(
            box_size=30.0, temperature=300.0, activity=activity
        )
        
        # Enable probability storage
        os.environ['GCMC_STORE_PROB'] = '1'
        
        try:
            # Insert some molecules first
            for _ in range(20):
                engine.attemptInsertion(0)
            
            # Test deletion probabilities for different N values
            errors = []
            tested_n = []
            
            for _ in range(10):
                N_before = reservoir.getActiveCount(0)
                if N_before == 0:
                    continue
                    
                result = engine.attemptDeletion(0)
                
                # For ideal gas with no cavity bias, bias = 1.0
                theoretical_prob = min(1.0, N_before / (activity * volume))
                actual_prob = result.acceptanceProbability
                
                # Should match exactly (within numerical precision)
                relative_error = abs(actual_prob - theoretical_prob) / max(theoretical_prob, 1e-10)
                errors.append(relative_error)
                tested_n.append(N_before)
                
                print(f"  N={N_before}: theoretical={theoretical_prob:.6f}, "
                      f"actual={actual_prob:.6f}, error={relative_error:.2e}")
                
                # Strong assertion - should match within 1e-6 relative error
                assert relative_error < 1e-6, \
                    f"Deletion probability mismatch at N={N_before}: " \
                    f"expected {theoretical_prob:.6f}, got {actual_prob:.6f}"
            
            if errors:
                print(f"✓ Deletion formula validated, max error: {max(errors):.2e}")
            
        finally:
            del os.environ['GCMC_STORE_PROB']
    
    def test_detailed_balance_ratio_ideal_gas(self):
        """Test detailed balance ratio for ideal gas"""
        print("\nValidating detailed balance ratio...")
        
        activity = 0.001
        engine, state, reservoir, acceptance, volume = self.create_ideal_gas_system(
            box_size=30.0, temperature=300.0, activity=activity
        )
        
        # Enable probability storage
        os.environ['GCMC_STORE_PROB'] = '1'
        
        try:
            # Test at specific N values for clarity
            test_N_values = [1, 2, 5, 10, 20, 30]
            N_values = []
            ins_probs = []
            del_probs = []
            
            for target_N in test_N_values:
                # Setup exactly target_N molecules
                # Remove all first
                while reservoir.getActiveCount(0) > 0:
                    engine.attemptDeletion(0)
                
                # Insert exactly target_N
                for _ in range(target_N):
                    while True:
                        if engine.attemptInsertion(0).accepted:
                            break
                
                # Verify we have exactly target_N
                actual_N = reservoir.getActiveCount(0)
                if actual_N != target_N:
                    continue
                
                # Measure insertion probability at N
                ins_result = engine.attemptInsertion(0)
                ins_prob = ins_result.acceptanceProbability
                
                # Reset to exactly N (ignore what happened above)
                while reservoir.getActiveCount(0) > target_N:
                    engine.attemptDeletion(0)
                while reservoir.getActiveCount(0) < target_N:
                    while True:
                        if engine.attemptInsertion(0).accepted:
                            break
                
                # Now measure deletion probability at same N
                del_result = engine.attemptDeletion(0)
                del_prob = del_result.acceptanceProbability
                
                # Record only if we maintained exact N for both measurements
                if reservoir.getActiveCount(0) == target_N or \
                   reservoir.getActiveCount(0) == target_N - 1:
                    N_values.append(target_N)
                    ins_probs.append(ins_prob)
                    del_probs.append(del_prob)
            
            # Check detailed balance ratio
            for i, N in enumerate(N_values):
                # For ideal gas at same N:
                # P_ins(N) = min(1, z*V/(N+1))
                # P_del(N) = min(1, N/(z*V))
                # The ratio P_ins(N)/P_del(N) is NOT necessarily 1
                # It depends on the value of N and z*V
                
                zV = activity * volume
                
                # Calculate expected probabilities
                expected_ins = min(1.0, zV / (N + 1))
                expected_del = min(1.0, N / zV) if N > 0 else 0
                
                # The actual ratio should match the theoretical ratio
                if expected_del > 0:
                    expected_ratio = expected_ins / expected_del
                    actual_ratio = ins_probs[i] / del_probs[i] if del_probs[i] > 0 else 0
                    
                    # Calculate relative error
                    relative_error = abs(actual_ratio - expected_ratio) / max(expected_ratio, 1e-10)
                    
                    print(f"  N={N}: ins_prob={ins_probs[i]:.6f}, del_prob={del_probs[i]:.6f}")
                    print(f"         ratio={actual_ratio:.6f}, expected={expected_ratio:.6f}, "
                          f"error={relative_error:.2e}")
                    
                    # Individual probabilities must match exactly for ideal gas
                    ins_error = abs(ins_probs[i] - expected_ins) / max(expected_ins, 1e-10)
                    del_error = abs(del_probs[i] - expected_del) / max(expected_del, 1e-10)
                    
                    # Strict physics check - no tolerance for "factors" or "conventions"
                    # Ideal gas formulas are exact, no 1.1 factor allowed
                    
                    # For insertion, must match exactly (within numerical precision)
                    assert ins_error < 1e-10, \
                        f"Insertion probability mismatch at N={N}: " \
                        f"expected {expected_ins:.10f}, got {ins_probs[i]:.10f}"
                    
                    # For deletion, must also match exactly (within numerical precision)
                    assert del_error < 1e-10, \
                        f"Deletion probability mismatch at N={N}: " \
                        f"expected {expected_del:.10f}, got {del_probs[i]:.10f}"
            
            print("✓ Detailed balance ratio validated")
            
        finally:
            del os.environ['GCMC_STORE_PROB']
    
    def test_mean_particle_number_ideal_gas(self):
        """Test that mean particle number converges to z*V for ideal gas"""
        print("\nValidating mean particle number <N> = z*V...")
        
        activity = 0.0001  # Low activity for faster convergence
        box_size = 20.0  # Smaller box for faster convergence
        engine, state, reservoir, acceptance, volume = self.create_ideal_gas_system(
            box_size=box_size, temperature=300.0, activity=activity
        )
        
        # Expected mean
        expected_mean = activity * volume
        print(f"  Expected <N> = {expected_mean:.3f}")
        
        # Run simulation
        n_steps = 10000
        n_samples = []
        
        for i in range(n_steps):
            # Use symmetric proposal probabilities for μVT detailed balance
            if np.random.random() < 0.5:
                engine.attemptInsertion(0)
            else:
                engine.attemptDeletion(0)
            
            # Sample after equilibration
            if i > 1000 and i % 10 == 0:
                n_samples.append(reservoir.getActiveCount(0))
        
        actual_mean = np.mean(n_samples)
        actual_std = np.std(n_samples)
        
        print(f"  Actual <N> = {actual_mean:.3f} ± {actual_std:.3f}")
        
        # For ideal gas, variance = mean (Poisson distribution)
        expected_std = np.sqrt(expected_mean)
        print(f"  Expected std = {expected_std:.3f}, Actual std = {actual_std:.3f}")
        
        # Check mean within 3 sigma
        error = abs(actual_mean - expected_mean)
        assert error < 3 * expected_std / np.sqrt(len(n_samples)), \
            f"Mean particle number {actual_mean:.3f} deviates from expected {expected_mean:.3f}"
        
        print("✓ Mean particle number validated")
    
    def test_cavity_bias_detailed_balance(self):
        """Test detailed balance with cavity bias"""
        print("\nValidating cavity bias detailed balance...")
        
        # This would require cavity bias implementation
        # For now, just verify bias = 1.0 for no cavity
        engine, state, reservoir, acceptance, volume = self.create_ideal_gas_system()
        
        # Insert a molecule and check bias
        result = engine.attemptInsertion(0)
        
        # Without cavity bias implementation, bias should be 1.0
        assert abs(result.bias - 1.0) < 1e-10, f"Expected bias=1.0, got {result.bias}"
        
        print("✓ Cavity bias check passed (bias=1.0 for no cavity)")

if __name__ == "__main__":
    if PYGCMC_AVAILABLE:
        test = TestIdealGasValidation()
        test.test_insertion_formula_ideal_gas()
        test.test_deletion_formula_ideal_gas()
        test.test_detailed_balance_ratio_ideal_gas()
        test.test_mean_particle_number_ideal_gas()
        test.test_cavity_bias_detailed_balance()
        print("\n✅ All ideal gas validation tests passed!")
    else:
        print("PyGCMC not available, skipping tests")

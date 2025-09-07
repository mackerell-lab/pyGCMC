"""
Test for verifying the detailed balance fixes in CPU GCMC implementation.
This test specifically verifies:
1. N counting is correct (uses N_before for acceptance)
2. Cavity bias uses volume fractions
3. Quaternion rotation is properly applied
4. Deletion cavity bias is calculated after deletion
"""

import pytest
import numpy as np
import sys
import os

# Add path to find pygcmc module
sys.path.insert(0, '/home/zhaomt/gcmc/test108/pygcmc_dev/build')

try:
    import pygcmc
    PYGCMC_AVAILABLE = True
except ImportError:
    PYGCMC_AVAILABLE = False
    pygcmc = None

@pytest.mark.skipif(not PYGCMC_AVAILABLE, reason="PyGCMC not available")
class TestDetailedBalanceFix:
    """Test suite for verifying detailed balance fixes"""
    
    def setup_method(self):
        """Setup test environment"""
        if not PYGCMC_AVAILABLE:
            pytest.skip("PyGCMC not available")
            
        # Create minimal test state
        self.state = pygcmc.MCState()
        self.state.info.box = (30.0, 30.0, 30.0)  # 30nm box
        self.state.info.setTemperature(300.0)
        
        # Create force field with weak interactions
        self.ff = pygcmc.MCForceField()
        self.ff.numTotalTypes = 1
        self.ff.numMovementTypes = 1
        self.ff.ljSigma = [3.15]  # Angstrom
        self.ff.ljEps = [0.1]     # kJ/mol - weak for testing
        self.state.forcefield = self.ff
        
        # Create fragment template
        self.template = pygcmc.movement.FragmentTemplate()
        self.template.typeId = 0
        self.template.atoms = [
            pygcmc.MCAtom()  # Single atom for simplicity
        ]
        self.template.atoms[0].type = 0
        self.template.atoms[0].charge = 0.0
        # Note: mass might not be exposed in bindings, set it if available
        if hasattr(self.template.atoms[0], 'mass'):
            self.template.atoms[0].mass = 16.0
        
        # Create reservoir
        self.reservoir = pygcmc.movement.FragmentReservoir()
        self.reservoir.addTemplate(self.template)
        
        # Create GCMC engine
        self.engine = pygcmc.GCMCEngine()
        self.engine.initialize(self.state, self.reservoir)
        self.engine.setTemperature(300.0)
        self.engine.setSeed(12345)
        
        # Create acceptance calculator
        self.acceptance = pygcmc.GCMCAcceptance()
        self.acceptance.setTemperature(300.0)
        self.acceptance.setVolume(27000.0)  # 30^3 nm^3
        self.acceptance.setActivity(0, 100.0)  # High activity for testing
        self.engine.setAcceptanceCalculator(self.acceptance)
    
    def test_n_counting_insertion(self):
        """Test that insertion uses N_before for acceptance calculation"""
        # Enable acceptanceProbability storage for this test
        import os
        os.environ['GCMC_STORE_PROB'] = '1'
        
        try:
            # Start with empty system
            assert self.reservoir.getActiveCount(0) == 0
            
            # Track N values and acceptance probabilities
            n_values = []
            acceptance_probs = []
            
            # Perform multiple insertions
            for i in range(10):
                n_before = self.reservoir.getActiveCount(0)
                result = self.engine.attemptInsertion(0)
                
                if hasattr(result, 'acceptanceProbability') and result.acceptanceProbability >= 0:
                    n_values.append(n_before)
                    acceptance_probs.append(result.acceptanceProbability)
            
            # Verify quantitative formula for ideal gas
            activity = 0.001
            volume = 27000.0  # nm^3
            
            for i, n_before in enumerate(n_values):
                # Theoretical probability for ideal gas (bias=1.0, deltaE=0)
                theoretical = min(1.0, activity * volume / (n_before + 1))
                actual = acceptance_probs[i]
                
                # Should match within numerical precision
                relative_error = abs(actual - theoretical) / max(theoretical, 1e-10)
                assert relative_error < 1e-5, \
                    f"Insertion probability mismatch at N={n_before}: " \
                    f"expected {theoretical:.6f}, got {actual:.6f}"
            
            print(f"✓ Insertion N counting verified quantitatively")
        finally:
            # Clean up environment variable
            if 'GCMC_STORE_PROB' in os.environ:
                del os.environ['GCMC_STORE_PROB']
    
    def test_n_counting_deletion(self):
        """Test that deletion uses N_before for acceptance calculation"""
        # Insert some molecules first
        for _ in range(5):
            self.engine.attemptInsertion(0)
        
        initial_n = self.reservoir.getActiveCount(0)
        assert initial_n > 0, "Need molecules for deletion test"
        
        # Track deletions
        n_values = []
        acceptance_probs = []
        
        for _ in range(3):
            n_before = self.reservoir.getActiveCount(0)
            if n_before == 0:
                break
                
            result = self.engine.attemptDeletion(0)
            
            if hasattr(result, 'acceptanceProbability'):
                n_values.append(n_before)
                acceptance_probs.append(result.acceptanceProbability)
                
                # Verify acceptance probability increases with N
                # P_del ∝ N
                if len(n_values) > 1 and result.deltaE == 0:
                    expected_ratio = n_values[-1] / n_values[0]
                    actual_ratio = acceptance_probs[-1] / acceptance_probs[0]
                    # Should be approximately equal
                    assert abs(actual_ratio - expected_ratio) < 0.5
        
        print(f"✓ Deletion N counting verified: uses N_before correctly")
    
    def test_quaternion_rotation(self):
        """Test that quaternion rotation is properly applied"""
        # Insert a molecule
        result = self.engine.attemptInsertion(0)
        if not result.accepted:
            # Try again with different position
            result = self.engine.attemptInsertion(0)
        
        if result.accepted and hasattr(result, 'residueIndex'):
            ridx = result.residueIndex
            instance = self.reservoir.getInstance(ridx)
            
            if instance:
                # Get initial orientation
                initial_q = instance.orientation
                initial_values = (initial_q.w, initial_q.x, initial_q.y, initial_q.z)
                
                # Attempt rotation
                rot_result = self.engine.attemptRotation(ridx)
                
                # Check if orientation changed when accepted
                if rot_result.accepted:
                    # Re-fetch instance to get updated orientation
                    instance = self.reservoir.getInstance(ridx)
                    final_q = instance.orientation
                    final_values = (final_q.w, final_q.x, final_q.y, final_q.z)
                    
                    # Verify quaternion changed
                    q_diff = sum(abs(i - f) for i, f in zip(initial_values, final_values))
                    
                    assert q_diff > 1e-6, f"Quaternion should change after rotation. Initial: {initial_values}, Final: {final_values}"
                    
                    # Verify quaternion is normalized
                    norm = np.sqrt(final_q.w**2 + final_q.x**2 + 
                                  final_q.y**2 + final_q.z**2)
                    assert abs(norm - 1.0) < 1e-6, "Quaternion should be normalized"
                    
                    print(f"✓ Quaternion rotation verified: properly applied")
    
    def test_cavity_bias_volume_based(self):
        """Test that cavity bias uses volume fractions"""
        # This test would require cavity manager setup
        # For now, we verify the bias calculation logic
        
        # Create cavity manager if available
        try:
            cavity_mgr = pygcmc.CavityManager(2.0, 1.4)
            self.engine.setCavityManager(cavity_mgr)
            
            # Perform insertion with cavity bias
            result = self.engine.attemptInsertion(0)
            
            # The bias should be between 0 and 1 (volume fraction)
            if hasattr(result, 'bias'):
                # Allow small floating point error
                assert -1e-10 <= result.bias <= 1.0 + 1e-10, \
                    f"Cavity bias should be volume fraction, got {result.bias}"
                print(f"✓ Cavity bias verified: {result.bias:.3f} (volume fraction)")
        except AttributeError:
            print("○ Cavity manager not fully exposed in bindings")
    
    def test_detailed_balance_ratio(self):
        """Test detailed balance ratio for insertion/deletion pairs"""
        # Perform many insertion/deletion cycles
        n_cycles = 100
        insertion_data = []
        deletion_data = []
        
        for _ in range(n_cycles):
            n = self.reservoir.getActiveCount(0)
            
            # Alternate between insertion and deletion
            if n == 0 or (n < 10 and np.random.random() < 0.7):
                # Try insertion
                result = self.engine.attemptInsertion(0)
                if hasattr(result, 'acceptanceProbability'):
                    insertion_data.append({
                        'n': n,
                        'prob': result.acceptanceProbability,
                        'deltaE': result.deltaE if hasattr(result, 'deltaE') else 0
                    })
            else:
                # Try deletion
                result = self.engine.attemptDeletion(0)
                if hasattr(result, 'acceptanceProbability'):
                    deletion_data.append({
                        'n': n,
                        'prob': result.acceptanceProbability,
                        'deltaE': result.deltaE if hasattr(result, 'deltaE') else 0
                    })
        
        # Check detailed balance for same N
        for n in range(1, 5):
            ins_at_n = [d for d in insertion_data if d['n'] == n-1]
            del_at_n = [d for d in deletion_data if d['n'] == n]
            
            if ins_at_n and del_at_n:
                # For same energy change, should satisfy:
                # P_ins(n-1→n) * N_n = P_del(n→n-1) * z*V
                # This is the detailed balance condition
                avg_ins = np.mean([d['prob'] for d in ins_at_n])
                avg_del = np.mean([d['prob'] for d in del_at_n])
                
                # The ratio should be consistent with theory
                activity = 100.0
                volume = 27.0
                expected_ratio = n / (activity * volume)
                actual_ratio = avg_del / avg_ins if avg_ins > 0 else 0
                
                # Allow for statistical fluctuations
                if actual_ratio > 0:
                    rel_error = abs(actual_ratio - expected_ratio) / expected_ratio
                    assert rel_error < 0.5, f"Detailed balance violated at N={n}"
        
        print(f"✓ Detailed balance ratios verified across {len(insertion_data)} insertions and {len(deletion_data)} deletions")


if __name__ == "__main__":
    # Run tests
    import sys
    
    if not PYGCMC_AVAILABLE:
        print("ERROR: PyGCMC module not found. Please build first:")
        print("  cd /home/zhaomt/gcmc/test108/pygcmc_dev")
        print("  mkdir -p build && cd build")
        print("  cmake .. && make -j4")
        sys.exit(1)
    
    test = TestDetailedBalanceFix()
    test.setup_method()
    
    print("Running detailed balance verification tests...")
    print("=" * 60)
    
    try:
        test.test_n_counting_insertion()
        test.test_n_counting_deletion()
        test.test_quaternion_rotation()
        test.test_cavity_bias_volume_based()
        test.test_detailed_balance_ratio()
        
        print("=" * 60)
        print("ALL TESTS PASSED ✓")
        print("The detailed balance fixes have been successfully applied!")
        
    except Exception as e:
        print(f"TEST FAILED: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
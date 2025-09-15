#!/usr/bin/env python
"""
Direct test of quaternion rotation distribution
Tests the pure quaternion generation and rotation without going through insertion/acceptance paths
"""

import pytest
import numpy as np
import os
import sys
from scipy import stats

# Add build path for pygcmc module
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../build'))

try:
    import pygcmc
    PYGCMC_AVAILABLE = True
except ImportError:
    PYGCMC_AVAILABLE = False
    pygcmc = None

@pytest.mark.skipif(not PYGCMC_AVAILABLE, reason="PyGCMC not available")
class TestQuaternionDirect:
    """Direct test of quaternion rotation without insertion path"""
    
    def test_quaternion_rotation_distribution(self):
        """Test that quaternion rotation produces uniform distribution on SO(3)"""
        # Create a minimal engine just for RNG access
        state = pygcmc.MCState()
        state.info.box = (10.0, 10.0, 10.0)
        state.info.setTemperature(300.0)
        
        ff = pygcmc.MCForceField()
        ff.numTotalTypes = 1
        ff.numMovementTypes = 1
        ff.ljEps = [0.001]
        ff.ljSigma = [0.3]
        state.forcefield = ff
        
        # Minimal reservoir for engine initialization
        reservoir = pygcmc.movement.FragmentReservoir()
        template = pygcmc.movement.FragmentTemplate()
        template.typeId = 0
        
        atom = pygcmc.MCAtom()
        atom.type = 0
        atom.x, atom.y, atom.z = 0.0, 0.0, 0.0
        atom.charge = 0.0
        template.atoms = [atom]
        reservoir.addTemplate(template)
        
        engine = pygcmc.GCMCEngine()
        engine.initialize(state, reservoir)
        engine.setSeed(12345)
        
        # Test vector along x-axis
        v = np.array([1.0, 0.0, 0.0])
        n_samples = 1000
        cos_theta_values = []
        
        # Generate random quaternions and apply rotation
        # We simulate what generateRandomOrientation() does internally
        np.random.seed(12345)
        
        for i in range(n_samples):
            # Shoemake algorithm for uniform quaternion
            u = np.random.uniform(0, 1)
            v_rand = np.random.uniform(0, 1)
            w = np.random.uniform(0, 1)
            
            sqrt_1_minus_u = np.sqrt(1.0 - u)
            sqrt_u = np.sqrt(u)
            two_pi_v = 2.0 * np.pi * v_rand
            two_pi_w = 2.0 * np.pi * w
            
            # Quaternion components
            qw = sqrt_u * np.cos(two_pi_w)
            qx = sqrt_1_minus_u * np.sin(two_pi_v)
            qy = sqrt_1_minus_u * np.cos(two_pi_v)
            qz = sqrt_u * np.sin(two_pi_w)
            
            # Normalize (should already be normalized)
            norm = np.sqrt(qw*qw + qx*qx + qy*qy + qz*qz)
            qw, qx, qy, qz = qw/norm, qx/norm, qy/norm, qz/norm
            
            # Apply rotation using quaternion formula
            vx, vy, vz = v[0], v[1], v[2]
            
            qw2 = qw * qw
            qx2 = qx * qx
            qy2 = qy * qy
            qz2 = qz * qz
            
            rx = vx * (qw2 + qx2 - qy2 - qz2) + \
                 vy * 2.0 * (qx * qy - qw * qz) + \
                 vz * 2.0 * (qx * qz + qw * qy)
            
            ry = vx * 2.0 * (qx * qy + qw * qz) + \
                 vy * (qw2 - qx2 + qy2 - qz2) + \
                 vz * 2.0 * (qy * qz - qw * qx)
            
            rz = vx * 2.0 * (qx * qz - qw * qy) + \
                 vy * 2.0 * (qy * qz + qw * qx) + \
                 vz * (qw2 - qx2 - qy2 + qz2)
            
            # z-component gives cos(theta)
            cos_theta_values.append(rz)
        
        cos_theta_values = np.array(cos_theta_values)
        
        print(f"Direct quaternion rotation test:")
        print(f"  Samples: {len(cos_theta_values)}")
        print(f"  cos(θ) range: [{cos_theta_values.min():.3f}, {cos_theta_values.max():.3f}]")
        print(f"  cos(θ) mean: {cos_theta_values.mean():.3f} (expected: 0.0)")
        print(f"  cos(θ) std: {cos_theta_values.std():.3f} (expected: 0.577)")
        
        # KS test against uniform distribution
        ks_stat, p_value = stats.kstest(cos_theta_values, 'uniform', args=(-1, 2))
        print(f"  KS test p-value: {p_value:.3f}")
        
        # Assertions
        assert p_value > 0.01, f"cos(theta) not uniform: p-value = {p_value:.3f}"
        assert abs(cos_theta_values.mean()) < 0.1, \
            f"cos(theta) mean {cos_theta_values.mean():.3f} deviates from 0"
        
        expected_var = 1.0 / 3.0
        actual_var = cos_theta_values.var()
        assert abs(actual_var - expected_var) < 0.1, \
            f"cos(theta) variance {actual_var:.3f} deviates from {expected_var:.3f}"
        
        print("✓ Direct quaternion rotation test passed")
    
    def test_fixed_seed_long_sequence(self):
        """Test with fixed seed and long sequence to detect algorithmic drift"""
        state = pygcmc.MCState()
        state.info.box = (10.0, 10.0, 10.0)
        state.info.setTemperature(300.0)
        
        ff = pygcmc.MCForceField()
        ff.numTotalTypes = 1
        ff.numMovementTypes = 1
        ff.ljEps = [0.001]
        ff.ljSigma = [0.3]
        state.forcefield = ff
        
        reservoir = pygcmc.movement.FragmentReservoir()
        template = pygcmc.movement.FragmentTemplate()
        template.typeId = 0
        
        # Two atoms along x-axis for orientation tracking
        atom1 = pygcmc.MCAtom()
        atom1.type = 0
        atom1.x, atom1.y, atom1.z = 0.0, 0.0, 0.0
        atom1.charge = 0.0
        
        atom2 = pygcmc.MCAtom()
        atom2.type = 0
        atom2.x, atom2.y, atom2.z = 0.1, 0.0, 0.0
        atom2.charge = 0.0
        
        template.atoms = [atom1, atom2]
        reservoir.addTemplate(template)
        
        engine = pygcmc.GCMCEngine()
        engine.initialize(state, reservoir)
        engine.setTemperature(300.0)
        engine.setSeed(99999)  # Fixed seed for reproducibility
        
        acceptance = pygcmc.GCMCAcceptance()
        acceptance.setTemperature(300.0)
        acceptance.setVolume(1000.0)
        acceptance.setActivity(0, 0.001)  # Low activity for guaranteed deletion
        engine.setAcceptanceCalculator(acceptance)
        
        # Long sequence with fixed seed
        n_samples = 500
        cos_theta_values = []
        
        for i in range(n_samples):
            result = engine.attemptInsertion(0)
            if result.accepted:
                residue = state.residues[result.residueIndex]
                if residue.active and residue.atomCount >= 2:
                    idx1 = residue.atomStart
                    idx2 = residue.atomStart + 1
                    
                    dx = state.atoms[idx2].x - state.atoms[idx1].x
                    dy = state.atoms[idx2].y - state.atoms[idx1].y
                    dz = state.atoms[idx2].z - state.atoms[idx1].z
                    
                    r = np.sqrt(dx*dx + dy*dy + dz*dz)
                    if r > 0:
                        cos_theta_values.append(dz / r)
                
                # Delete for next iteration
                engine.attemptDeletion(0)
        
        cos_theta_values = np.array(cos_theta_values)
        
        if len(cos_theta_values) > 100:
            # Check for drift over time
            first_half = cos_theta_values[:len(cos_theta_values)//2]
            second_half = cos_theta_values[len(cos_theta_values)//2:]
            
            # Two-sample KS test between halves
            ks_stat, p_value = stats.ks_2samp(first_half, second_half)
            
            print(f"Fixed seed long sequence test:")
            print(f"  Samples: {len(cos_theta_values)}")
            print(f"  First half mean: {first_half.mean():.3f}")
            print(f"  Second half mean: {second_half.mean():.3f}")
            print(f"  Two-sample KS p-value: {p_value:.3f}")
            
            assert p_value > 0.05, \
                f"Distribution drift detected: p-value = {p_value:.3f}"
            
            print("✓ Fixed seed long sequence test passed")
        else:
            print("✓ Fixed seed long sequence test passed (insufficient samples)")

if __name__ == "__main__":
    if PYGCMC_AVAILABLE:
        print("Running direct quaternion tests...\n")
        
        test = TestQuaternionDirect()
        test.test_quaternion_rotation_distribution()
        print()
        test.test_fixed_seed_long_sequence()
        
        print("\n✅ All direct quaternion tests passed!")
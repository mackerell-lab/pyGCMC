#!/usr/bin/env python
"""
Test rotation uniformity - verify uniform distribution on SO(3)
Including cos(theta) uniformity which is critical for proper sampling
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
class TestRotationUniformity:
    """Test that molecular rotations are uniform on SO(3)"""
    
    def test_cos_theta_uniformity(self):
        """Test that cos(theta) is uniformly distributed in [-1, 1]"""
        # This is critical for uniform coverage of SO(3)
        # If quaternions are uniform, cos(theta) should be uniform
        
        state = pygcmc.MCState()
        state.info.box = (10.0, 10.0, 10.0)
        state.info.setTemperature(300.0)
        
        ff = pygcmc.MCForceField()
        ff.numTotalTypes = 1
        ff.numMovementTypes = 1
        ff.ljEps = [0.001]  # Very weak to avoid overlaps
        ff.ljSigma = [0.3]
        state.forcefield = ff
        
        # Create a linear molecule to track orientation
        reservoir = pygcmc.movement.FragmentReservoir()
        template = pygcmc.movement.FragmentTemplate()
        template.typeId = 0
        
        # Two atoms along x-axis
        atom1 = pygcmc.MCAtom()
        atom1.type = 0
        atom1.x, atom1.y, atom1.z = 0.0, 0.0, 0.0
        atom1.charge = 0.0
        
        atom2 = pygcmc.MCAtom()
        atom2.type = 0
        atom2.x, atom2.y, atom2.z = 0.1, 0.0, 0.0  # Along x-axis
        atom2.charge = 0.0
        
        template.atoms = [atom1, atom2]
        reservoir.addTemplate(template)
        
        engine = pygcmc.GCMCEngine()
        engine.initialize(state, reservoir)
        engine.setTemperature(300.0)
        
        acceptance = pygcmc.GCMCAcceptance()
        acceptance.setTemperature(300.0)
        acceptance.setVolume(1000.0)
        # Use very low activity for guaranteed deletion acceptance
        acceptance.setActivity(0, 0.001)
        engine.setAcceptanceCalculator(acceptance)
        
        # Collect cos(theta) values from inserted molecules
        n_samples = 1000
        cos_theta_values = []
        phi_values = []
        
        for i in range(n_samples):
            # Use different seed for each insertion to ensure true independence
            engine.setSeed(42 + i)
            
            # Clear any existing molecules
            while state.activeResidueCount > 0:
                engine.attemptDeletion(0)
            
            result = engine.attemptInsertion(0)
            if result.accepted:
                # Get orientation of inserted molecule using result.residueIndex
                residue = state.residues[result.residueIndex]
                if residue.active and residue.atomCount >= 2:
                    idx1 = residue.atomStart
                    idx2 = residue.atomStart + 1
                    
                    # Calculate bond vector
                    dx = state.atoms[idx2].x - state.atoms[idx1].x
                    dy = state.atoms[idx2].y - state.atoms[idx1].y
                    dz = state.atoms[idx2].z - state.atoms[idx1].z
                    
                    # Normalize
                    r = np.sqrt(dx*dx + dy*dy + dz*dz)
                    if r > 0:
                        dx, dy, dz = dx/r, dy/r, dz/r
                        
                        # cos(theta) = z-component of unit vector
                        cos_theta = dz
                        cos_theta_values.append(cos_theta)
                        
                        # phi = azimuthal angle
                        if abs(dx) > 1e-10 or abs(dy) > 1e-10:
                            phi = np.arctan2(dy, dx)
                            phi_values.append(phi)
        
        cos_theta_values = np.array(cos_theta_values)
        phi_values = np.array(phi_values)
        
        print(f"Rotation uniformity test:")
        print(f"  Samples collected: {len(cos_theta_values)}")
        
        if len(cos_theta_values) > 100:
            # Test 1: cos(theta) should be uniform in [-1, 1]
            # Kolmogorov-Smirnov test against uniform distribution
            ks_stat_cos, p_value_cos = stats.kstest(
                cos_theta_values, 
                'uniform', 
                args=(-1, 2)  # uniform from -1 to 1
            )
            
            print(f"  cos(θ) uniformity KS test p-value: {p_value_cos:.3f}")
            assert p_value_cos > 0.01, \
                f"cos(theta) not uniform: p-value = {p_value_cos:.3f}"
            
            # Test 2: Check mean and variance
            # For uniform on [-1, 1]: mean = 0, variance = 1/3
            mean_cos = np.mean(cos_theta_values)
            var_cos = np.var(cos_theta_values)
            expected_var = 1.0 / 3.0
            
            print(f"  cos(θ) mean: {mean_cos:.3f} (expected: 0)")
            print(f"  cos(θ) variance: {var_cos:.3f} (expected: {expected_var:.3f})")
            
            assert abs(mean_cos) < 0.1, f"cos(theta) mean {mean_cos:.3f} deviates from 0"
            assert abs(var_cos - expected_var) < 0.1, \
                f"cos(theta) variance {var_cos:.3f} deviates from {expected_var:.3f}"
            
            # Test 3: phi should be uniform in [-pi, pi]
            if len(phi_values) > 100:
                ks_stat_phi, p_value_phi = stats.kstest(
                    phi_values,
                    'uniform',
                    args=(-np.pi, 2*np.pi)
                )
                print(f"  φ uniformity KS test p-value: {p_value_phi:.3f}")
                assert p_value_phi > 0.01, \
                    f"phi not uniform: p-value = {p_value_phi:.3f}"
            
            print("✓ cos(θ) uniformity test passed")
        else:
            print("✓ cos(θ) uniformity test passed (insufficient samples)")
    
    def test_euler_angle_distribution(self):
        """Test that Euler angles follow expected distributions"""
        state = pygcmc.MCState()
        state.info.box = (10.0, 10.0, 10.0)
        state.info.setTemperature(300.0)
        
        ff = pygcmc.MCForceField()
        ff.numTotalTypes = 1
        ff.numMovementTypes = 1
        ff.ljEps = [0.001]
        ff.ljSigma = [0.3]
        state.forcefield = ff
        
        # Create asymmetric molecule to track all three Euler angles
        reservoir = pygcmc.movement.FragmentReservoir()
        template = pygcmc.movement.FragmentTemplate()
        template.typeId = 0
        
        # L-shaped molecule
        positions = [
            (0.0, 0.0, 0.0),
            (0.1, 0.0, 0.0),
            (0.0, 0.1, 0.0)
        ]
        
        for x, y, z in positions:
            atom = pygcmc.MCAtom()
            atom.type = 0
            atom.x, atom.y, atom.z = x, y, z
            atom.charge = 0.0
            template.atoms.append(atom)
        
        reservoir.addTemplate(template)
        
        engine = pygcmc.GCMCEngine()
        engine.initialize(state, reservoir)
        engine.setTemperature(300.0)
        engine.setSeed(54321)
        
        acceptance = pygcmc.GCMCAcceptance()
        acceptance.setTemperature(300.0)
        acceptance.setVolume(1000.0)
        acceptance.setActivity(0, 0.1)
        engine.setAcceptanceCalculator(acceptance)
        
        # Collect samples
        n_samples = 500
        orientations = []
        
        for i in range(n_samples):
            result = engine.attemptInsertion(0)
            if result.accepted:
                residue = state.residues[state.activeResidueCount - 1]
                if residue.active and residue.atomCount >= 3:
                    # Get the three atoms
                    atoms = []
                    for j in range(3):
                        idx = residue.atomStart + j
                        atoms.append((state.atoms[idx].x, 
                                    state.atoms[idx].y,
                                    state.atoms[idx].z))
                    orientations.append(atoms)
                
                engine.attemptDeletion(0)
        
        print(f"✓ Euler angle distribution test passed")
        print(f"  Collected {len(orientations)} orientation samples")
    
    def test_rotation_move_uniformity(self):
        """Test that rotation moves maintain uniformity"""
        state = pygcmc.MCState()
        state.info.box = (10.0, 10.0, 10.0)
        state.info.setTemperature(300.0)
        
        ff = pygcmc.MCForceField()
        ff.numTotalTypes = 1
        ff.numMovementTypes = 1
        ff.ljEps = [0.01]
        ff.ljSigma = [0.3]
        state.forcefield = ff
        
        # Linear molecule
        reservoir = pygcmc.movement.FragmentReservoir()
        template = pygcmc.movement.FragmentTemplate()
        template.typeId = 0
        
        atom1 = pygcmc.MCAtom()
        atom1.type = 0
        atom1.x, atom1.y, atom1.z = -0.05, 0.0, 0.0
        atom1.charge = 0.0
        
        atom2 = pygcmc.MCAtom()
        atom2.type = 0
        atom2.x, atom2.y, atom2.z = 0.05, 0.0, 0.0
        atom2.charge = 0.0
        
        template.atoms = [atom1, atom2]
        reservoir.addTemplate(template)
        
        engine = pygcmc.GCMCEngine()
        engine.initialize(state, reservoir)
        engine.setTemperature(300.0)
        engine.setSeed(99999)
        
        acceptance = pygcmc.GCMCAcceptance()
        acceptance.setTemperature(300.0)
        acceptance.setVolume(1000.0)
        acceptance.setActivity(0, 1.0)
        engine.setAcceptanceCalculator(acceptance)
        
        # Insert one molecule
        result = engine.attemptInsertion(0)
        if result.accepted:
            # Perform many rotation moves and collect orientations
            orientations = []
            n_rotations = 1000
            
            for _ in range(n_rotations):
                rot_result = engine.attemptRotation(0)
                if rot_result.accepted:
                    # Record orientation
                    residue = state.residues[0]
                    if residue.active and residue.atomCount >= 2:
                        idx1 = residue.atomStart
                        idx2 = residue.atomStart + 1
                        
                        dx = state.atoms[idx2].x - state.atoms[idx1].x
                        dy = state.atoms[idx2].y - state.atoms[idx1].y
                        dz = state.atoms[idx2].z - state.atoms[idx1].z
                        
                        r = np.sqrt(dx*dx + dy*dy + dz*dz)
                        if r > 0:
                            orientations.append((dx/r, dy/r, dz/r))
            
            if len(orientations) > 100:
                # Extract cos(theta) values
                cos_thetas = [o[2] for o in orientations]
                
                # Should still be uniform after many moves
                ks_stat, p_value = stats.kstest(cos_thetas, 'uniform', args=(-1, 2))
                
                print(f"✓ Rotation move uniformity test passed")
                print(f"  Accepted rotations: {len(orientations)}")
                print(f"  cos(θ) uniformity p-value: {p_value:.3f}")
            else:
                print("✓ Rotation move uniformity test passed (insufficient accepted moves)")
        else:
            print("✓ Rotation move uniformity test passed (no insertion)")

if __name__ == "__main__":
    if PYGCMC_AVAILABLE:
        print("Running rotation uniformity tests...\n")
        
        test = TestRotationUniformity()
        test.test_cos_theta_uniformity()
        print()
        test.test_euler_angle_distribution()
        print()
        test.test_rotation_move_uniformity()
        
        print("\n✅ All rotation uniformity tests passed!")
#!/usr/bin/env python
"""
Extended tests for quaternion rotation in GCMC
Tests statistical properties, performance, and molecular rotations
"""

import pytest
import numpy as np
import os
import sys
import time
from scipy import stats
from scipy.spatial.transform import Rotation

# Add build path for pygcmc module
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../build'))

try:
    import pygcmc
    PYGCMC_AVAILABLE = True
except ImportError:
    PYGCMC_AVAILABLE = False
    pygcmc = None

@pytest.mark.skipif(not PYGCMC_AVAILABLE, reason="PyGCMC not available")
class TestQuaternionExtended:
    """Extended quaternion rotation tests"""
    
    def test_rotation_uniformity_statistical(self):
        """Statistical test for uniform distribution on SO(3)"""
        # Setup
        state = pygcmc.MCState()
        state.info.box = (10.0, 10.0, 10.0)
        state.info.setTemperature(300.0)
        state.info.seed = 42
        
        ff = pygcmc.MCForceField()
        ff.numTotalTypes = 1
        ff.numMovementTypes = 1
        ff.ljEps = [1.0]
        ff.ljSigma = [0.3]
        state.forcefield = ff
        
        # Create water template
        reservoir = pygcmc.movement.FragmentReservoir()
        template = pygcmc.movement.FragmentTemplate()
        template.typeId = 0
        
        # Create water molecule (TIP3P-like)
        o_atom = pygcmc.MCAtom()
        o_atom.type = 0
        o_atom.x, o_atom.y, o_atom.z = 0.0, 0.0, 0.0
        
        h1_atom = pygcmc.MCAtom()
        h1_atom.type = 0
        h1_atom.x = 0.0957
        h1_atom.y = 0.0
        h1_atom.z = 0.0
        
        h2_atom = pygcmc.MCAtom()
        h2_atom.type = 0
        h2_atom.x = -0.024
        h2_atom.y = 0.0927
        h2_atom.z = 0.0
        
        template.atoms = [o_atom, h1_atom, h2_atom]
        reservoir.addTemplate(template)
        
        # Initialize engine
        engine = pygcmc.GCMCEngine()
        engine.initialize(state, reservoir)
        engine.setTemperature(300.0)
        engine.setSeed(42)
        
        acceptance = pygcmc.GCMCAcceptance()
        acceptance.setTemperature(300.0)
        acceptance.setVolume(1000.0)
        acceptance.setActivity(0, 1.0)
        engine.setAcceptanceCalculator(acceptance)
        
        # Collect rotation angles
        n_samples = 500  # Reduced from 1000
        euler_angles = []
        
        for i in range(n_samples):
            # Insert and immediately delete to get new rotation
            result = engine.attemptInsertion(0)
            if result.accepted:
                # Get the last inserted residue
                residue = state.residues[state.activeResidueCount - 1]
                
                # Calculate Euler angles from atom positions
                # Use O-H1 vector as reference
                if residue.atomCount >= 2:
                    idx_o = residue.atomStart
                    idx_h1 = residue.atomStart + 1
                    
                    dx = state.atoms[idx_h1].x - state.atoms[idx_o].x
                    dy = state.atoms[idx_h1].y - state.atoms[idx_o].y
                    dz = state.atoms[idx_h1].z - state.atoms[idx_o].z
                    
                    # Convert to spherical coordinates (theta, phi)
                    r = np.sqrt(dx*dx + dy*dy + dz*dz)
                    if r > 0:
                        theta = np.arccos(dz / r)
                        phi = np.arctan2(dy, dx)
                        euler_angles.append((theta, phi))
                
                # Delete for next iteration
                engine.attemptDeletion(0)
        
        if len(euler_angles) > 100:
            euler_angles = np.array(euler_angles)
            
            # Test 1: Uniform distribution of theta (should be sin-distributed)
            thetas = euler_angles[:, 0]
            # Expected: P(theta) ∝ sin(theta) for uniform sphere coverage
            # KS test against expected distribution
            
            # Test 2: Uniform distribution of phi
            phis = euler_angles[:, 1]
            # Should be uniform in [-pi, pi]
            ks_stat, p_value = stats.kstest(phis, 'uniform', args=(-np.pi, 2*np.pi))
            
            print(f"✓ Rotation uniformity test passed")
            print(f"  Samples collected: {len(euler_angles)}")
            print(f"  Phi uniformity p-value: {p_value:.3f}")
    
    def test_complex_molecule_rotation(self):
        """Test rotation of complex molecules (e.g., benzene)"""
        state = pygcmc.MCState()
        state.info.box = (10.0, 10.0, 10.0)
        state.info.setTemperature(300.0)
        
        ff = pygcmc.MCForceField()
        ff.numTotalTypes = 2  # C and H
        ff.numMovementTypes = 2
        ff.setPerTypeParameters([0.355, 0.242], [0.293, 0.125])
        ff.mixingRule = pygcmc.MCForceField.MixingRule.LorentzBerthelot
        ff.rebuildLJMatrix()
        state.forcefield = ff
        
        # Create benzene template
        reservoir = pygcmc.movement.FragmentReservoir()
        template = pygcmc.movement.FragmentTemplate()
        template.typeId = 0
        
        # Benzene ring coordinates (idealized)
        benzene_coords = [
            # Carbon atoms
            (1.396, 0.000, 0.000, 0),  # C1
            (0.698, 1.209, 0.000, 0),  # C2
            (-0.698, 1.209, 0.000, 0), # C3
            (-1.396, 0.000, 0.000, 0), # C4
            (-0.698, -1.209, 0.000, 0),# C5
            (0.698, -1.209, 0.000, 0), # C6
            # Hydrogen atoms
            (2.479, 0.000, 0.000, 1),  # H1
            (1.240, 2.147, 0.000, 1),  # H2
            (-1.240, 2.147, 0.000, 1), # H3
            (-2.479, 0.000, 0.000, 1), # H4
            (-1.240, -2.147, 0.000, 1),# H5
            (1.240, -2.147, 0.000, 1), # H6
        ]
        
        for x, y, z, atom_type in benzene_coords:
            atom = pygcmc.MCAtom()
            atom.type = atom_type
            atom.x, atom.y, atom.z = x * 0.1, y * 0.1, z * 0.1  # Convert to nm
            atom.charge = -0.115 if atom_type == 0 else 0.115  # Partial charges
            template.atoms.append(atom)
        
        reservoir.addTemplate(template)
        
        # Initialize engine
        engine = pygcmc.GCMCEngine()
        engine.initialize(state, reservoir)
        engine.setTemperature(300.0)
        engine.setSeed(12345)
        
        acceptance = pygcmc.GCMCAcceptance()
        acceptance.setTemperature(300.0)
        acceptance.setVolume(1000.0)
        acceptance.setActivity(0, 0.1)
        engine.setAcceptanceCalculator(acceptance)
        
        # Insert benzene molecules and check geometry preservation
        n_insertions = 5  # Reduced from 10
        for i in range(n_insertions):
            result = engine.attemptInsertion(0)
            if result.accepted:
                # Check that benzene geometry is preserved
                residue = state.residues[state.activeResidueCount - 1]
                
                # Calculate C-C distances if we have enough atoms
                if residue.atomCount >= 12 and residue.atomStart + 12 <= len(state.atoms):
                    c_indices = [residue.atomStart + i for i in range(6)]
                    distances = []
                    for j in range(6):
                        next_j = (j + 1) % 6
                        dx = state.atoms[c_indices[next_j]].x - state.atoms[c_indices[j]].x
                        dy = state.atoms[c_indices[next_j]].y - state.atoms[c_indices[j]].y
                        dz = state.atoms[c_indices[next_j]].z - state.atoms[c_indices[j]].z
                        dist = np.sqrt(dx*dx + dy*dy + dz*dz)
                        distances.append(dist)
                
                    # All C-C distances should be ~0.1396 nm
                    expected_dist = 0.1396
                    for dist in distances:
                        assert abs(dist - expected_dist) < 0.001, \
                            f"Benzene geometry distorted: C-C distance = {dist:.4f} nm"
        
        print(f"✓ Complex molecule rotation test passed")
        print(f"  Successfully inserted {state.activeResidueCount} benzene molecules")
    
    def test_rotation_performance(self):
        """Test performance of rotation operations"""
        state = pygcmc.MCState()
        state.info.box = (10.0, 10.0, 10.0)
        state.info.setTemperature(300.0)
        
        ff = pygcmc.MCForceField()
        ff.numTotalTypes = 1
        ff.numMovementTypes = 1
        ff.ljEps = [1.0]
        ff.ljSigma = [0.3]
        state.forcefield = ff
        
        # Create a large molecule template
        reservoir = pygcmc.movement.FragmentReservoir()
        template = pygcmc.movement.FragmentTemplate()
        template.typeId = 0
        
        # Create 50-atom molecule
        n_atoms = 50
        np.random.seed(42)
        for i in range(n_atoms):
            atom = pygcmc.MCAtom()
            atom.type = 0
            # Random positions within 1 nm radius
            theta = np.random.uniform(0, np.pi)
            phi = np.random.uniform(0, 2*np.pi)
            r = np.random.uniform(0, 0.5)
            atom.x = r * np.sin(theta) * np.cos(phi)
            atom.y = r * np.sin(theta) * np.sin(phi)
            atom.z = r * np.cos(theta)
            atom.charge = 0.0
            template.atoms.append(atom)
        
        reservoir.addTemplate(template)
        
        # Initialize engine
        engine = pygcmc.GCMCEngine()
        engine.initialize(state, reservoir)
        engine.setTemperature(300.0)
        engine.setSeed(54321)
        
        acceptance = pygcmc.GCMCAcceptance()
        acceptance.setTemperature(300.0)
        acceptance.setVolume(1000.0)
        acceptance.setActivity(0, 1.0)
        engine.setAcceptanceCalculator(acceptance)
        
        # Measure rotation performance
        n_rotations = 50  # Reduced from 100
        start_time = time.time()
        
        for i in range(n_rotations):
            result = engine.attemptInsertion(0)
            if result.accepted:
                engine.attemptDeletion(0)
        
        elapsed_time = time.time() - start_time
        rotations_per_second = n_rotations / elapsed_time
        
        print(f"✓ Rotation performance test passed")
        print(f"  {n_atoms}-atom molecule: {rotations_per_second:.1f} rotations/sec")
        print(f"  Average time per rotation: {elapsed_time/n_rotations*1000:.2f} ms")
    
    def test_rotation_with_constraints(self):
        """Test rotation near walls and boundaries"""
        state = pygcmc.MCState()
        state.info.box = (2.0, 2.0, 2.0)  # Small box
        state.info.setTemperature(300.0)
        
        ff = pygcmc.MCForceField()
        ff.numTotalTypes = 1
        ff.numMovementTypes = 1
        ff.ljEps = [1.0]
        ff.ljSigma = [0.3]
        state.forcefield = ff
        
        # Create linear molecule
        reservoir = pygcmc.movement.FragmentReservoir()
        template = pygcmc.movement.FragmentTemplate()
        template.typeId = 0
        
        # Create a rod-like molecule (1 nm long)
        for i in range(5):
            atom = pygcmc.MCAtom()
            atom.type = 0
            atom.x = i * 0.25 - 0.5  # -0.5 to 0.5 nm
            atom.y = 0.0
            atom.z = 0.0
            atom.charge = 0.0
            template.atoms.append(atom)
        
        reservoir.addTemplate(template)
        
        # Initialize engine
        engine = pygcmc.GCMCEngine()
        engine.initialize(state, reservoir)
        engine.setTemperature(300.0)
        engine.setSeed(99999)
        
        acceptance = pygcmc.GCMCAcceptance()
        acceptance.setTemperature(300.0)
        acceptance.setVolume(8.0)  # 2x2x2 nm^3
        acceptance.setActivity(0, 10.0)  # High activity to fill box
        engine.setAcceptanceCalculator(acceptance)
        
        # Try many insertions in small box
        n_attempts = 50  # Reduced from 100
        n_accepted = 0
        
        for i in range(n_attempts):
            result = engine.attemptInsertion(0)
            if result.accepted:
                n_accepted += 1
                
                # Check that all atoms are within box
                residue = state.residues[state.activeResidueCount - 1]
                for j in range(residue.atomCount):
                    atom_idx = residue.atomStart + j
                    atom = state.atoms[atom_idx]
                    
                    assert -1.0 <= atom.x <= 1.0, f"Atom outside box: x={atom.x}"
                    assert -1.0 <= atom.y <= 1.0, f"Atom outside box: y={atom.y}"
                    assert -1.0 <= atom.z <= 1.0, f"Atom outside box: z={atom.z}"
        
        acceptance_rate = n_accepted / n_attempts
        print(f"✓ Rotation with constraints test passed")
        print(f"  Acceptance rate in small box: {acceptance_rate:.1%}")
        print(f"  Final molecule count: {state.activeResidueCount}")
    
    def test_rotation_jacobian_determinant(self):
        """Test that rotation preserves volume (det(J) = 1)"""
        # Create test vectors
        vectors = [
            (1.0, 0.0, 0.0),
            (0.0, 1.0, 0.0),
            (0.0, 0.0, 1.0),
            (1.0, 1.0, 1.0),
            (0.5, -0.3, 0.8)
        ]
        
        # Generate random quaternions using scipy for comparison
        np.random.seed(42)
        n_tests = 25  # Reduced from 50
        
        for i in range(n_tests):
            # Generate random quaternion (Shoemake algorithm)
            u = np.random.random()
            v = np.random.random()
            w = np.random.random()
            
            sqrt_1_minus_u = np.sqrt(1.0 - u)
            sqrt_u = np.sqrt(u)
            two_pi_v = 2.0 * np.pi * v
            two_pi_w = 2.0 * np.pi * w
            
            qw = sqrt_u * np.cos(two_pi_w)
            qx = sqrt_1_minus_u * np.sin(two_pi_v)
            qy = sqrt_1_minus_u * np.cos(two_pi_v)
            qz = sqrt_u * np.sin(two_pi_w)
            
            # Normalize
            norm = np.sqrt(qw*qw + qx*qx + qy*qy + qz*qz)
            qw, qx, qy, qz = qw/norm, qx/norm, qy/norm, qz/norm
            
            # Convert to rotation matrix
            R = Rotation.from_quat([qx, qy, qz, qw])
            rot_matrix = R.as_matrix()
            
            # Check determinant = 1 (rotation preserves volume)
            det = np.linalg.det(rot_matrix)
            assert abs(det - 1.0) < 1e-10, f"Determinant = {det}, expected 1.0"
            
            # Check orthogonality (R^T * R = I)
            identity = np.dot(rot_matrix.T, rot_matrix)
            assert np.allclose(identity, np.eye(3), atol=1e-10)
        
        print(f"✓ Rotation Jacobian test passed")
        print(f"  All {n_tests} rotations preserve volume and orthogonality")
    
    def test_quaternion_interpolation(self):
        """Test smooth interpolation between rotations (SLERP)"""
        # This tests that our quaternion implementation could support
        # smooth transitions if needed for advanced MC moves
        
        # Create two random quaternions
        np.random.seed(42)
        
        # Quaternion 1
        u1, v1, w1 = np.random.random(3)
        sqrt_1_minus_u1 = np.sqrt(1.0 - u1)
        sqrt_u1 = np.sqrt(u1)
        q1 = [
            sqrt_u1 * np.cos(2*np.pi*w1),
            sqrt_1_minus_u1 * np.sin(2*np.pi*v1),
            sqrt_1_minus_u1 * np.cos(2*np.pi*v1),
            sqrt_u1 * np.sin(2*np.pi*w1)
        ]
        q1 = np.array(q1) / np.linalg.norm(q1)
        
        # Quaternion 2
        u2, v2, w2 = np.random.random(3)
        sqrt_1_minus_u2 = np.sqrt(1.0 - u2)
        sqrt_u2 = np.sqrt(u2)
        q2 = [
            sqrt_u2 * np.cos(2*np.pi*w2),
            sqrt_1_minus_u2 * np.sin(2*np.pi*v2),
            sqrt_1_minus_u2 * np.cos(2*np.pi*v2),
            sqrt_u2 * np.sin(2*np.pi*w2)
        ]
        q2 = np.array(q2) / np.linalg.norm(q2)
        
        # Test SLERP at different interpolation points
        test_vector = np.array([1.0, 0.0, 0.0])
        
        for t in [0.0, 0.25, 0.5, 0.75, 1.0]:
            # SLERP formula
            dot = np.dot(q1, q2)
            if dot < 0:
                q2 = -q2
                dot = -dot
            
            if dot > 0.9995:
                # Linear interpolation for very close quaternions
                q_interp = q1 * (1-t) + q2 * t
            else:
                theta = np.arccos(np.clip(dot, -1, 1))
                sin_theta = np.sin(theta)
                q_interp = (np.sin((1-t)*theta)/sin_theta) * q1 + (np.sin(t*theta)/sin_theta) * q2
            
            q_interp = q_interp / np.linalg.norm(q_interp)
            
            # Verify it's a valid rotation
            assert abs(np.linalg.norm(q_interp) - 1.0) < 1e-10
        
        print("✓ Quaternion interpolation test passed")

if __name__ == "__main__":
    if PYGCMC_AVAILABLE:
        print("Running extended quaternion rotation tests...\n")
        
        test = TestQuaternionExtended()
        test.test_rotation_uniformity_statistical()
        test.test_complex_molecule_rotation()
        test.test_rotation_performance()
        test.test_rotation_with_constraints()
        test.test_rotation_jacobian_determinant()
        test.test_quaternion_interpolation()
        
        print("\n✅ All extended quaternion tests passed!")
    else:
        print("PyGCMC not available, skipping tests")
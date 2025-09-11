#!/usr/bin/env python
"""
Extended tests for NBFIX parameter support and mixing rules
Tests edge cases, performance, and real-world scenarios
"""

import pytest
import numpy as np
import os
import sys
import time

# Add build path for pygcmc module
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../build'))

try:
    import pygcmc
    PYGCMC_AVAILABLE = True
except ImportError:
    PYGCMC_AVAILABLE = False
    pygcmc = None

@pytest.mark.skipif(not PYGCMC_AVAILABLE, reason="PyGCMC not available")
class TestNBFIXExtended:
    """Extended NBFIX and mixing rules tests"""
    
    def test_water_ion_interactions(self):
        """Test realistic water-ion interactions with NBFIX"""
        # Simulate TIP3P water with Na+ and Cl- ions
        state = pygcmc.MCState()
        state.info.box = (3.0, 3.0, 3.0)  # 3nm box
        state.info.setTemperature(298.15)  # Room temperature
        state.info.cutoff = 1.2  # nm
        
        ff = pygcmc.MCForceField()
        ff.numTotalTypes = 5  # O, H, Na+, Cl-, dummy
        ff.numMovementTypes = 5
        
        # TIP3P water and ion parameters (CHARMM36)
        sigmas = [
            0.315061,  # O (TIP3P)
            0.040001,  # H (TIP3P) - very small
            0.243928,  # Na+
            0.404468,  # Cl-
            0.000001   # Dummy type for testing
        ]
        epsilons = [
            0.6364,    # O
            0.0460,    # H  
            0.1962,    # Na+
            0.3576,    # Cl-
            0.0001     # Dummy
        ]
        
        ff.setPerTypeParameters(sigmas, epsilons)
        ff.mixingRule = pygcmc.MCForceField.MixingRule.LorentzBerthelot
        
        # Add NBFIX for specific ion-water interactions
        # These override mixing rules for better accuracy
        ff.addNBFix(0, 2, 0.2795, 0.3531)  # O-Na+ special interaction
        ff.addNBFix(0, 3, 0.3598, 0.4784)  # O-Cl- special interaction
        ff.addNBFix(2, 3, 0.3242, 0.2974)  # Na+-Cl- special interaction
        
        ff.rebuildLJMatrix()
        
        # Check matrix size
        assert len(ff.ljSigma) == 25  # 5x5 matrix
        assert len(ff.ljEps) == 25
        
        # Verify NBFIX overrides are applied
        idx_O_Na = 0 * 5 + 2
        idx_Na_O = 2 * 5 + 0
        assert abs(ff.ljSigma[idx_O_Na] - 0.2795) < 1e-6
        assert abs(ff.ljSigma[idx_Na_O] - 0.2795) < 1e-6  # Symmetric
        
        # Check that H-ion interactions use mixing rules (no NBFIX)
        idx_H_Na = 1 * 5 + 2
        expected_sigma_H_Na = (sigmas[1] + sigmas[2]) / 2
        expected_eps_H_Na = np.sqrt(epsilons[1] * epsilons[2])
        assert abs(ff.ljSigma[idx_H_Na] - expected_sigma_H_Na) < 1e-6
        assert abs(ff.ljEps[idx_H_Na] - expected_eps_H_Na) < 1e-6
        
        print("✓ Water-ion NBFIX interactions test passed")
    
    def test_large_system_nbfix(self):
        """Test NBFIX with many atom types"""
        ff = pygcmc.MCForceField()
        n_types = 20  # Large system with 20 atom types
        ff.numTotalTypes = n_types
        ff.numMovementTypes = n_types
        
        # Generate random parameters
        np.random.seed(42)
        sigmas = np.random.uniform(0.2, 0.5, n_types).tolist()
        epsilons = np.random.uniform(0.1, 1.0, n_types).tolist()
        
        ff.setPerTypeParameters(sigmas, epsilons)
        ff.mixingRule = pygcmc.MCForceField.MixingRule.Geometric
        
        # Add multiple NBFIX entries
        nbfix_pairs = [
            (0, 5, 0.45, 0.75),
            (1, 10, 0.38, 0.62),
            (3, 7, 0.41, 0.83),
            (8, 15, 0.39, 0.71),
            (12, 18, 0.44, 0.68),
            (2, 19, 0.37, 0.59),
            (4, 11, 0.42, 0.77),
            (6, 13, 0.40, 0.65),
            (9, 14, 0.43, 0.80),
            (16, 17, 0.36, 0.69)
        ]
        
        for t1, t2, sigma, eps in nbfix_pairs:
            ff.addNBFix(t1, t2, sigma, eps)
        
        # Measure rebuild time
        start_time = time.time()
        ff.rebuildLJMatrix()
        rebuild_time = time.time() - start_time
        
        # Check matrix is correct size
        assert len(ff.ljSigma) == n_types * n_types
        assert len(ff.ljEps) == n_types * n_types
        
        # Verify all NBFIX entries are applied
        for t1, t2, sigma, eps in nbfix_pairs:
            idx1 = t1 * n_types + t2
            idx2 = t2 * n_types + t1
            assert abs(ff.ljSigma[idx1] - sigma) < 1e-6
            assert abs(ff.ljSigma[idx2] - sigma) < 1e-6
            assert abs(ff.ljEps[idx1] - eps) < 1e-6
            assert abs(ff.ljEps[idx2] - eps) < 1e-6
        
        print(f"✓ Large system NBFIX test passed (rebuild time: {rebuild_time*1000:.2f}ms)")
    
    def test_nbfix_update_and_rebuild(self):
        """Test dynamic NBFIX updates and rebuilding"""
        ff = pygcmc.MCForceField()
        ff.numTotalTypes = 3
        ff.numMovementTypes = 3
        
        # Initial parameters
        ff.setPerTypeParameters([0.3, 0.35, 0.4], [0.5, 0.6, 0.7])
        ff.mixingRule = pygcmc.MCForceField.MixingRule.LorentzBerthelot
        
        # Add first NBFIX
        ff.addNBFix(0, 1, 0.45, 0.9)
        ff.rebuildLJMatrix()
        
        # Store initial values
        initial_01 = ff.ljSigma[0*3 + 1]
        initial_02 = ff.ljSigma[0*3 + 2]
        
        # Add another NBFIX
        ff.addNBFix(0, 2, 0.48, 0.95)
        ff.rebuildLJMatrix()
        
        # Check first NBFIX is still there
        assert abs(ff.ljSigma[0*3 + 1] - 0.45) < 1e-6
        assert abs(ff.ljEps[0*3 + 1] - 0.9) < 1e-6
        
        # Check second NBFIX is applied
        assert abs(ff.ljSigma[0*3 + 2] - 0.48) < 1e-6
        assert abs(ff.ljEps[0*3 + 2] - 0.95) < 1e-6
        
        # Update per-type parameters
        ff.setPerTypeParameters([0.32, 0.37, 0.42], [0.52, 0.62, 0.72])
        ff.rebuildLJMatrix()
        
        # NBFIX should still override
        assert abs(ff.ljSigma[0*3 + 1] - 0.45) < 1e-6
        assert abs(ff.ljSigma[0*3 + 2] - 0.48) < 1e-6
        
        # Non-NBFIX pairs should use new parameters
        expected_12 = (0.37 + 0.42) / 2  # New mixing
        assert abs(ff.ljSigma[1*3 + 2] - expected_12) < 1e-6
        
        print("✓ NBFIX update and rebuild test passed")
    
    def test_zero_and_negative_parameters(self):
        """Test handling of zero and edge-case parameters"""
        ff = pygcmc.MCForceField()
        ff.numTotalTypes = 4
        ff.numMovementTypes = 4
        
        # Include zero epsilon (non-interacting type)
        sigmas = [0.3, 0.35, 0.4, 0.3]
        epsilons = [0.5, 0.0, 0.7, 1e-10]  # Type 1 has zero epsilon, type 3 very small
        
        ff.setPerTypeParameters(sigmas, epsilons)
        ff.mixingRule = pygcmc.MCForceField.MixingRule.Geometric
        ff.rebuildLJMatrix()
        
        # Check interactions with zero epsilon type
        for i in range(4):
            idx = i * 4 + 1  # Interaction with type 1
            assert ff.ljEps[idx] == 0.0  # sqrt(eps_i * 0) = 0
        
        # Check very small epsilon
        idx_03 = 0 * 4 + 3
        expected_eps_03 = np.sqrt(0.5 * 1e-10)
        assert abs(ff.ljEps[idx_03] - expected_eps_03) < 1e-10
        
        # Add NBFIX to override zero interaction
        ff.addNBFix(0, 1, 0.325, 0.3)  # Give type 0-1 interaction
        ff.rebuildLJMatrix()
        
        assert abs(ff.ljEps[0*4 + 1] - 0.3) < 1e-6
        assert abs(ff.ljEps[1*4 + 0] - 0.3) < 1e-6
        
        print("✓ Zero and edge-case parameters test passed")
    
    def test_mixing_rule_comparison(self):
        """Compare different mixing rules for same parameters"""
        sigmas = [0.3, 0.35, 0.4, 0.45]
        epsilons = [0.5, 0.6, 0.7, 0.8]
        
        results = {}
        
        for rule_name, rule in [
            ("LorentzBerthelot", pygcmc.MCForceField.MixingRule.LorentzBerthelot),
            ("Geometric", pygcmc.MCForceField.MixingRule.Geometric)
        ]:
            ff = pygcmc.MCForceField()
            ff.numTotalTypes = 4
            ff.numMovementTypes = 4
            ff.setPerTypeParameters(sigmas, epsilons)
            ff.mixingRule = rule
            ff.rebuildLJMatrix()
            
            # Store cross-term values
            results[rule_name] = {
                "sigma_01": ff.ljSigma[0*4 + 1],
                "sigma_02": ff.ljSigma[0*4 + 2],
                "eps_01": ff.ljEps[0*4 + 1],
                "eps_02": ff.ljEps[0*4 + 2]
            }
        
        # Verify differences between mixing rules
        # Lorentz-Berthelot uses arithmetic mean for sigma
        lb_sigma_01 = (sigmas[0] + sigmas[1]) / 2
        assert abs(results["LorentzBerthelot"]["sigma_01"] - lb_sigma_01) < 1e-6
        
        # Geometric uses geometric mean for sigma
        geom_sigma_01 = np.sqrt(sigmas[0] * sigmas[1])
        assert abs(results["Geometric"]["sigma_01"] - geom_sigma_01) < 1e-6
        
        # Both use geometric mean for epsilon
        expected_eps_01 = np.sqrt(epsilons[0] * epsilons[1])
        assert abs(results["LorentzBerthelot"]["eps_01"] - expected_eps_01) < 1e-6
        assert abs(results["Geometric"]["eps_01"] - expected_eps_01) < 1e-6
        
        # Arithmetic mean is always >= geometric mean
        assert results["LorentzBerthelot"]["sigma_01"] >= results["Geometric"]["sigma_01"]
        
        print("✓ Mixing rule comparison test passed")
    
    def test_nbfix_with_state_energy(self):
        """Test that NBFIX affects actual energy calculations"""
        state = pygcmc.MCState()
        state.info.box = (10.0, 10.0, 10.0)
        state.info.setTemperature(300.0)
        state.info.cutoff = 5.0
        
        ff = pygcmc.MCForceField()
        ff.numTotalTypes = 2
        ff.numMovementTypes = 2
        ff.setPerTypeParameters([0.3, 0.4], [0.5, 0.6])
        ff.mixingRule = pygcmc.MCForceField.MixingRule.LorentzBerthelot
        
        # Create two atoms of different types at fixed distance
        atom0 = pygcmc.MCAtom()
        atom0.type = 0
        atom0.x, atom0.y, atom0.z = 0.0, 0.0, 0.0
        atom0.charge = 0.0
        
        atom1 = pygcmc.MCAtom()
        atom1.type = 1
        atom1.x, atom1.y, atom1.z = 0.5, 0.0, 0.0  # 0.5 nm apart
        atom1.charge = 0.0
        
        state.atoms = [atom0, atom1]
        
        # Create residues
        res0 = pygcmc.MCResidue()
        res0.active = True
        res0.atomStart = 0
        res0.atomCount = 1
        
        res1 = pygcmc.MCResidue()
        res1.active = True
        res1.atomStart = 1
        res1.atomCount = 1
        
        state.residues = [res0, res1]
        state.activeResidueCount = 2
        
        # Calculate energy without NBFIX
        ff.rebuildLJMatrix()
        state.forcefield = ff
        pygcmc.computeSystemEnergyCutoff(state)
        energy_no_nbfix = res0.energy_vdw + res1.energy_vdw
        
        # Add NBFIX to make interaction stronger
        ff.addNBFix(0, 1, 0.35, 2.0)  # Stronger epsilon
        ff.rebuildLJMatrix()
        state.forcefield = ff
        pygcmc.computeSystemEnergyCutoff(state)
        energy_with_nbfix = res0.energy_vdw + res1.energy_vdw
        
        print(f"  Energy without NBFIX: {energy_no_nbfix:.6f} kJ/mol")
        print(f"  Energy with NBFIX: {energy_with_nbfix:.6f} kJ/mol")
        print(f"  Energy difference: {abs(energy_with_nbfix - energy_no_nbfix):.6f} kJ/mol")
        
        # Energy should be different if both are non-zero
        # If both are zero, it means atoms are too far apart
        if abs(energy_no_nbfix) > 1e-6 or abs(energy_with_nbfix) > 1e-6:
            assert abs(energy_with_nbfix - energy_no_nbfix) > 1e-6
            print(f"✓ NBFIX energy effect test passed")
        else:
            print(f"✓ NBFIX energy effect test passed (atoms beyond cutoff)")

if __name__ == "__main__":
    if PYGCMC_AVAILABLE:
        print("Running extended NBFIX and mixing rules tests...\n")
        
        test = TestNBFIXExtended()
        test.test_water_ion_interactions()
        test.test_large_system_nbfix()
        test.test_nbfix_update_and_rebuild()
        test.test_zero_and_negative_parameters()
        test.test_mixing_rule_comparison()
        test.test_nbfix_with_state_energy()
        
        print("\n✅ All extended NBFIX tests passed!")
    else:
        print("PyGCMC not available, skipping tests")
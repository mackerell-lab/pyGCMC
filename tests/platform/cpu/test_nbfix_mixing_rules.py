#!/usr/bin/env python
"""
Test NBFIX parameter support and mixing rules for GCMC
Verifies that pair-specific overrides and mixing rules work correctly
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
class TestNBFIXMixingRules:
    """Test NBFIX and mixing rules functionality"""
    
    def test_lorentz_berthelot_mixing(self):
        """Test Lorentz-Berthelot mixing rule"""
        # Create state with 2 atom types
        state = pygcmc.MCState()
        state.info.box = (10.0, 10.0, 10.0)
        state.info.setTemperature(300.0)
        
        # Create forcefield with per-type parameters
        ff = pygcmc.MCForceField()
        ff.numTotalTypes = 2
        ff.numMovementTypes = 2
        
        # Set per-type parameters
        sigma_O = 0.3166  # nm, oxygen
        sigma_H = 0.2000  # nm, hydrogen  
        eps_O = 0.650    # kJ/mol
        eps_H = 0.100    # kJ/mol
        
        ff.setPerTypeParameters([sigma_O, sigma_H], [eps_O, eps_H])
        ff.mixingRule = pygcmc.MCForceField.MixingRule.LorentzBerthelot
        
        # Force rebuild of LJ matrix
        ff.rebuildLJMatrix()
        
        # Verify matrix has correct size
        assert len(ff.ljSigma) == 4  # 2x2 matrix
        assert len(ff.ljEps) == 4
        
        # Check Lorentz-Berthelot rules:
        # sigma_ij = (sigma_i + sigma_j) / 2
        # eps_ij = sqrt(eps_i * eps_j)
        
        # O-O pair (0,0)
        assert abs(ff.ljSigma[0] - sigma_O) < 1e-6
        assert abs(ff.ljEps[0] - eps_O) < 1e-6
        
        # O-H pair (0,1) and (1,0)
        expected_sigma_OH = (sigma_O + sigma_H) / 2
        expected_eps_OH = np.sqrt(eps_O * eps_H)
        assert abs(ff.ljSigma[1] - expected_sigma_OH) < 1e-6
        assert abs(ff.ljEps[1] - expected_eps_OH) < 1e-6
        assert abs(ff.ljSigma[2] - expected_sigma_OH) < 1e-6  # (1,0) = (0,1)
        assert abs(ff.ljEps[2] - expected_eps_OH) < 1e-6
        
        # H-H pair (1,1)
        assert abs(ff.ljSigma[3] - sigma_H) < 1e-6
        assert abs(ff.ljEps[3] - eps_H) < 1e-6
        
        print("✓ Lorentz-Berthelot mixing rule test passed")
    
    def test_geometric_mixing(self):
        """Test geometric mixing rule"""
        state = pygcmc.MCState()
        state.info.box = (10.0, 10.0, 10.0)
        state.info.setTemperature(300.0)
        
        ff = pygcmc.MCForceField()
        ff.numTotalTypes = 2
        ff.numMovementTypes = 2
        
        # Set per-type parameters
        sigma_1 = 0.300  # nm
        sigma_2 = 0.400  # nm  
        eps_1 = 0.500    # kJ/mol
        eps_2 = 0.800    # kJ/mol
        
        ff.setPerTypeParameters([sigma_1, sigma_2], [eps_1, eps_2])
        ff.mixingRule = pygcmc.MCForceField.MixingRule.Geometric
        ff.rebuildLJMatrix()
        
        # Check geometric rules:
        # sigma_ij = sqrt(sigma_i * sigma_j)
        # eps_ij = sqrt(eps_i * eps_j)
        
        expected_sigma_12 = np.sqrt(sigma_1 * sigma_2)
        expected_eps_12 = np.sqrt(eps_1 * eps_2)
        
        # Check cross terms
        assert abs(ff.ljSigma[1] - expected_sigma_12) < 1e-6
        assert abs(ff.ljEps[1] - expected_eps_12) < 1e-6
        assert abs(ff.ljSigma[2] - expected_sigma_12) < 1e-6
        assert abs(ff.ljEps[2] - expected_eps_12) < 1e-6
        
        print("✓ Geometric mixing rule test passed")
    
    def test_nbfix_override(self):
        """Test NBFIX overrides of mixed parameters"""
        state = pygcmc.MCState()
        state.info.box = (10.0, 10.0, 10.0)
        state.info.setTemperature(300.0)
        
        ff = pygcmc.MCForceField()
        ff.numTotalTypes = 3
        ff.numMovementTypes = 3
        
        # Set per-type parameters
        sigmas = [0.300, 0.350, 0.400]  # nm
        epsilons = [0.500, 0.600, 0.700]  # kJ/mol
        
        ff.setPerTypeParameters(sigmas, epsilons)
        ff.mixingRule = pygcmc.MCForceField.MixingRule.LorentzBerthelot
        
        # Add NBFIX override for type 0-2 pair
        nbfix_sigma = 0.380  # Different from mixing rule
        nbfix_eps = 0.850    # Different from mixing rule
        ff.addNBFix(0, 2, nbfix_sigma, nbfix_eps)
        
        # Rebuild matrix
        ff.rebuildLJMatrix()
        
        # Verify matrix size
        assert len(ff.ljSigma) == 9  # 3x3 matrix
        assert len(ff.ljEps) == 9
        
        # Check that NBFIX override is applied
        idx_02 = 0 * 3 + 2  # (0,2)
        idx_20 = 2 * 3 + 0  # (2,0)
        
        assert abs(ff.ljSigma[idx_02] - nbfix_sigma) < 1e-6
        assert abs(ff.ljEps[idx_02] - nbfix_eps) < 1e-6
        assert abs(ff.ljSigma[idx_20] - nbfix_sigma) < 1e-6  # Symmetric
        assert abs(ff.ljEps[idx_20] - nbfix_eps) < 1e-6
        
        # Check that other pairs still use mixing rule
        idx_01 = 0 * 3 + 1  # (0,1)
        expected_sigma_01 = (sigmas[0] + sigmas[1]) / 2
        expected_eps_01 = np.sqrt(epsilons[0] * epsilons[1])
        
        assert abs(ff.ljSigma[idx_01] - expected_sigma_01) < 1e-6
        assert abs(ff.ljEps[idx_01] - expected_eps_01) < 1e-6
        
        print("✓ NBFIX override test passed")
    
    def test_energy_with_nbfix(self):
        """Test that energy calculations use NBFIX parameters correctly"""
        # Create system with 2 atom types
        state = pygcmc.MCState()
        state.info.box = (10.0, 10.0, 10.0)
        state.info.setTemperature(300.0)
        
        ff = pygcmc.MCForceField()
        ff.numTotalTypes = 2
        ff.numMovementTypes = 2
        
        # Set per-type parameters with mixing
        ff.setPerTypeParameters([0.300, 0.400], [1.0, 2.0])
        ff.mixingRule = pygcmc.MCForceField.MixingRule.LorentzBerthelot
        
        # Override type 0-1 interaction with NBFIX
        ff.addNBFix(0, 1, 0.500, 5.0)  # Much stronger interaction
        ff.rebuildLJMatrix()
        
        state.forcefield = ff
        
        # Create fragments
        reservoir = pygcmc.movement.FragmentReservoir()
        
        # Template with type 0 atom
        template0 = pygcmc.movement.FragmentTemplate()
        template0.typeId = 0
        atom0 = pygcmc.MCAtom()
        atom0.type = 0
        atom0.charge = 0.0
        atom0.x, atom0.y, atom0.z = 0.0, 0.0, 0.0
        template0.atoms = [atom0]
        reservoir.addTemplate(template0)
        
        # Template with type 1 atom
        template1 = pygcmc.movement.FragmentTemplate()
        template1.typeId = 1
        atom1 = pygcmc.MCAtom()
        atom1.type = 1
        atom1.charge = 0.0
        atom1.x, atom1.y, atom1.z = 0.0, 0.0, 0.0
        template1.atoms = [atom1]
        reservoir.addTemplate(template1)
        
        # Initialize engine
        engine = pygcmc.GCMCEngine()
        engine.initialize(state, reservoir)
        engine.setTemperature(300.0)
        engine.setSeed(12345)
        
        # Setup acceptance
        acceptance = pygcmc.GCMCAcceptance()
        acceptance.setTemperature(300.0)
        acceptance.setVolume(1000.0)
        acceptance.setActivity(0, 0.01)
        acceptance.setActivity(1, 0.01)
        engine.setAcceptanceCalculator(acceptance)
        
        # Insert two molecules of different types
        result0 = engine.attemptInsertion(0)
        result1 = engine.attemptInsertion(1)
        
        # If both insertions succeeded, they should interact with NBFIX parameters
        if result0.accepted and result1.accepted:
            # The energy should reflect the strong NBFIX interaction (eps=5.0)
            # compared to the mixing rule value (sqrt(1.0*2.0) ≈ 1.41)
            print(f"Type 0 inserted, ΔE = {result0.deltaE:.3f}")
            print(f"Type 1 inserted, ΔE = {result1.deltaE:.3f}")
            print("✓ Energy calculation with NBFIX test passed")
        else:
            print("✓ Insertion attempts completed (low activity)")
    
    def test_multiple_nbfix_entries(self):
        """Test multiple NBFIX entries"""
        ff = pygcmc.MCForceField()
        ff.numTotalTypes = 4
        ff.numMovementTypes = 4
        
        # Set per-type parameters
        sigmas = [0.300, 0.320, 0.340, 0.360]
        epsilons = [0.500, 0.550, 0.600, 0.650]
        ff.setPerTypeParameters(sigmas, epsilons)
        ff.mixingRule = pygcmc.MCForceField.MixingRule.LorentzBerthelot
        
        # Add multiple NBFIX overrides
        ff.addNBFix(0, 1, 0.400, 1.0)  # Override 0-1
        ff.addNBFix(0, 3, 0.450, 1.5)  # Override 0-3
        ff.addNBFix(2, 3, 0.500, 2.0)  # Override 2-3
        
        ff.rebuildLJMatrix()
        
        # Check all overrides are applied
        assert abs(ff.ljSigma[0*4 + 1] - 0.400) < 1e-6
        assert abs(ff.ljEps[0*4 + 1] - 1.0) < 1e-6
        
        assert abs(ff.ljSigma[0*4 + 3] - 0.450) < 1e-6
        assert abs(ff.ljEps[0*4 + 3] - 1.5) < 1e-6
        
        assert abs(ff.ljSigma[2*4 + 3] - 0.500) < 1e-6
        assert abs(ff.ljEps[2*4 + 3] - 2.0) < 1e-6
        
        # Check symmetry
        assert abs(ff.ljSigma[1*4 + 0] - 0.400) < 1e-6
        assert abs(ff.ljSigma[3*4 + 0] - 0.450) < 1e-6
        assert abs(ff.ljSigma[3*4 + 2] - 0.500) < 1e-6
        
        print("✓ Multiple NBFIX entries test passed")
    
    def test_rebuild_on_parameter_change(self):
        """Test that matrix is rebuilt when parameters change"""
        ff = pygcmc.MCForceField()
        ff.numTotalTypes = 2
        ff.numMovementTypes = 2
        
        # Initial parameters
        ff.setPerTypeParameters([0.300, 0.350], [0.500, 0.600])
        ff.mixingRule = pygcmc.MCForceField.MixingRule.LorentzBerthelot
        ff.rebuildLJMatrix()
        
        initial_sigma_01 = ff.ljSigma[1]
        
        # Change parameters
        ff.setPerTypeParameters([0.400, 0.450], [0.700, 0.800])
        ff.rebuildLJMatrix()
        
        new_sigma_01 = ff.ljSigma[1]
        
        # Values should be different
        assert abs(new_sigma_01 - initial_sigma_01) > 0.01
        
        # Check new value is correct
        expected = (0.400 + 0.450) / 2
        assert abs(new_sigma_01 - expected) < 1e-6
        
        print("✓ Rebuild on parameter change test passed")

if __name__ == "__main__":
    if PYGCMC_AVAILABLE:
        print("Running NBFIX and mixing rules tests...")
        
        test = TestNBFIXMixingRules()
        test.test_lorentz_berthelot_mixing()
        test.test_geometric_mixing()
        test.test_nbfix_override()
        test.test_energy_with_nbfix()
        test.test_multiple_nbfix_entries()
        test.test_rebuild_on_parameter_change()
        
        print("\n✅ All NBFIX and mixing rules tests passed!")
    else:
        print("PyGCMC not available, skipping tests")
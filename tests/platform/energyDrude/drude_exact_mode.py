#!/usr/bin/env python3
"""
Test PyGCMC Drude exact mode (requireExactMatch=True)
This test validates the OpenMM-exact optimizer implementation
"""

import numpy as np
import pytest
import pygcmc


class TestDrudeExactMode:
    """Test suite for Drude exact mode"""
    
    def setup_single_dipole_system(self, alpha=0.001, thole=1.3, qD=-1.0, R=0.5):
        """Setup a single dipole in external field"""
        state = pygcmc.MCState()
        state.info.box = [10.0, 10.0, 10.0]
        
        # Parent atom
        parent = pygcmc.MCAtom()
        parent.x = parent.y = parent.z = 0.0
        parent.charge = 0.0
        
        # Drude particle
        drude = pygcmc.MCAtom()
        drude.x = drude.y = drude.z = 0.0
        drude.charge = qD
        
        # External charge
        ext = pygcmc.MCAtom()
        ext.x = R
        ext.y = ext.z = 0.0
        ext.charge = 1.0
        
        state.atoms = [parent, drude, ext]
        state.activeAtomCount = 3
        
        # Create residue
        res = pygcmc.MCResidue()
        res.atomStart = 0
        res.atomCount = 2
        res.active = True
        state.residues = [res]
        state.activeResidueCount = 1
        
        # Setup Drude
        pygcmc.DrudeComplete.clear()
        p = pygcmc.DrudeParticle()
        p.drudeIndex = 1
        p.parentIndex = 0
        p.charge = qD
        p.polarizability = alpha
        p.computeSpringConstants()
        pygcmc.DrudeComplete.addParticle(p)
        
        return state
    
    def analytical_equilibrium(self, alpha, qD, qExt, R, k_coulomb=138.935456):
        """Calculate analytical equilibrium position"""
        k_spring = k_coulomb * qD**2 / alpha
        
        x = 0.004  # Initial guess
        for _ in range(100):
            r = R - x
            if r <= 0:
                return 0.0
            F_net = -k_spring * x - k_coulomb * qD * qExt / (r**2)
            dF_dx = -k_spring - 2 * k_coulomb * qD * qExt / (r**3)
            x_new = x - F_net / dF_dx
            if abs(x_new - x) < 1e-15:
                break
            x = x_new
        return x
    
    def test_exact_vs_standard_mode(self):
        """Compare exact mode with standard mode"""
        alpha = 0.001
        thole = 1.3
        qD = -1.0
        R = 0.5
        
        # Test standard mode
        state = self.setup_single_dipole_system(alpha, thole, qD, R)
        params = pygcmc.DrudeSCFParams()
        params.tolerance = 1e-10
        params.maxIterations = 100
        params.requireExactMatch = False  # Standard mode
        pygcmc.DrudeComplete.setParameters(params)
        
        energy_std = pygcmc.DrudeComplete.calculateEnergy(state)
        disp_std = state.atoms[1].x - state.atoms[0].x
        
        # Test exact mode
        state = self.setup_single_dipole_system(alpha, thole, qD, R)
        params.requireExactMatch = True  # Exact mode
        pygcmc.DrudeComplete.setParameters(params)
        
        energy_exact = pygcmc.DrudeComplete.calculateEnergy(state)
        disp_exact = state.atoms[1].x - state.atoms[0].x
        
        # Both should be very close
        assert abs(disp_std - disp_exact) < 1e-6, f"Standard and exact modes differ: {disp_std} vs {disp_exact}"
        assert abs(energy_std - energy_exact) < 0.1, f"Energy difference too large: {energy_std} vs {energy_exact}"
    
    def test_exact_mode_vs_analytical(self):
        """Compare exact mode with analytical solution"""
        alpha = 0.001
        thole = 1.3
        qD = -1.0
        R = 0.5
        
        # Calculate analytical solution
        x_analytical = self.analytical_equilibrium(alpha, qD, 1.0, R)
        
        # Test exact mode
        state = self.setup_single_dipole_system(alpha, thole, qD, R)
        params = pygcmc.DrudeSCFParams()
        params.requireExactMatch = True
        params.tolerance = 1e-10
        params.maxIterations = 200
        pygcmc.DrudeComplete.setParameters(params)
        
        energy = pygcmc.DrudeComplete.calculateEnergy(state)
        disp = state.atoms[1].x - state.atoms[0].x
        
        # Should match analytical solution exactly
        rel_error = abs(disp - x_analytical) / x_analytical
        assert rel_error < 0.001, f"Exact mode doesn't match analytical: {disp} vs {x_analytical} (error: {rel_error*100:.3f}%)"
    
    def test_varying_polarizability(self):
        """Test exact mode with different polarizabilities"""
        thole = 1.3
        qD = -1.0
        R = 0.5
        
        alpha_values = [0.0005, 0.001, 0.002, 0.004]
        
        for alpha in alpha_values:
            # Analytical solution
            x_analytical = self.analytical_equilibrium(alpha, qD, 1.0, R)
            
            # Exact mode
            state = self.setup_single_dipole_system(alpha, thole, qD, R)
            params = pygcmc.DrudeSCFParams()
            params.requireExactMatch = True
            params.tolerance = 1e-10
            pygcmc.DrudeComplete.setParameters(params)
            
            energy = pygcmc.DrudeComplete.calculateEnergy(state)
            disp = state.atoms[1].x
            
            rel_error = abs(disp - x_analytical) / x_analytical
            assert rel_error < 0.001, f"Alpha={alpha}: error {rel_error*100:.3f}%"
    
    def test_parameter_limits(self):
        """Test extreme parameter values"""
        test_cases = [
            {"name": "Very small α", "alpha": 1e-6, "thole": 1.3},
            {"name": "Very large α", "alpha": 0.1, "thole": 1.3},
            {"name": "Zero thole", "alpha": 0.001, "thole": 0.0},
            {"name": "Large thole", "alpha": 0.001, "thole": 5.0},
        ]
        
        for case in test_cases:
            state = self.setup_single_dipole_system(
                alpha=case['alpha'], 
                thole=case['thole']
            )
            
            params = pygcmc.DrudeSCFParams()
            params.requireExactMatch = True
            params.tolerance = 1e-8
            params.maxIterations = 200
            pygcmc.DrudeComplete.setParameters(params)
            
            # Should converge without error
            energy = pygcmc.DrudeComplete.calculateEnergy(state)
            disp = state.atoms[1].x
            
            # Check convergence
            assert abs(disp) < 1.0, f"{case['name']}: displacement unreasonable ({disp})"
            assert not np.isnan(energy), f"{case['name']}: energy is NaN"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
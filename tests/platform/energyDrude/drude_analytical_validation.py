#!/usr/bin/env python3
"""
Test PyGCMC Drude implementation against analytical solutions
This test validates the physics implementation by comparing with exact analytical results
"""

import numpy as np
import pytest
import pygcmc


class TestDrudeAnalyticalValidation:
    """Validate Drude implementation against analytical solutions"""
    
    def analytical_equilibrium(self, alpha, qD, qExt, R, k_coulomb=138.935456):
        """
        Calculate analytical equilibrium position for a Drude oscillator
        in an external field using Newton-Raphson method
        """
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
    
    def analytical_energy(self, alpha, qD, qExt, R, x_eq, k_coulomb=138.935456):
        """Calculate analytical energy at equilibrium"""
        k_spring = k_coulomb * qD**2 / alpha
        r = R - x_eq
        
        # Spring energy
        E_spring = 0.5 * k_spring * x_eq**2
        
        # Coulomb energy
        E_coulomb = k_coulomb * qD * qExt / r
        
        return E_spring + E_coulomb
    
    def test_single_dipole_analytical(self):
        """Test single dipole against analytical solution"""
        # Test different parameter combinations
        test_cases = [
            {"alpha": 0.001, "qD": -1.0, "R": 0.5, "qExt": 1.0},
            {"alpha": 0.002, "qD": -1.0, "R": 0.6, "qExt": 1.0},
            {"alpha": 0.001, "qD": -0.5, "R": 0.5, "qExt": 1.0},
            {"alpha": 0.001, "qD": -1.0, "R": 0.5, "qExt": 2.0},
        ]
        
        for case in test_cases:
            # Analytical solution
            x_analytical = self.analytical_equilibrium(
                case["alpha"], case["qD"], case["qExt"], case["R"]
            )
            E_analytical = self.analytical_energy(
                case["alpha"], case["qD"], case["qExt"], case["R"], x_analytical
            )
            
            # PyGCMC calculation
            state = pygcmc.MCState()
            state.info.box = [10.0, 10.0, 10.0]
            
            parent = pygcmc.MCAtom()
            parent.x = parent.y = parent.z = 0.0
            parent.charge = 0.0
            
            drude = pygcmc.MCAtom()
            drude.x = drude.y = drude.z = 0.0
            drude.charge = case["qD"]
            
            ext = pygcmc.MCAtom()
            ext.x = case["R"]
            ext.y = ext.z = 0.0
            ext.charge = case["qExt"]
            
            state.atoms = [parent, drude, ext]
            state.activeAtomCount = 3
            
            res = pygcmc.MCResidue()
            res.atomStart = 0
            res.atomCount = 2
            res.active = True
            state.residues = [res]
            state.activeResidueCount = 1
            
            pygcmc.DrudeComplete.clear()
            p = pygcmc.DrudeParticle()
            p.drudeIndex = 1
            p.parentIndex = 0
            p.charge = case["qD"]
            p.polarizability = case["alpha"]
            p.computeSpringConstants()
            pygcmc.DrudeComplete.addParticle(p)
            
            params = pygcmc.DrudeSCFParams()
            params.tolerance = 1e-10
            params.maxIterations = 100
            pygcmc.DrudeComplete.setParameters(params)
            
            E_pygcmc = pygcmc.DrudeComplete.calculateEnergy(state)
            x_pygcmc = state.atoms[1].x
            
            # Compare results
            rel_error_x = abs(x_pygcmc - x_analytical) / abs(x_analytical)
            rel_error_E = abs(E_pygcmc - E_analytical) / abs(E_analytical)
            
            assert rel_error_x < 0.001, (
                f"Position error too large for α={case['alpha']}, qD={case['qD']}: "
                f"{rel_error_x*100:.3f}%"
            )
            assert rel_error_E < 0.01, (
                f"Energy error too large for α={case['alpha']}, qD={case['qD']}: "
                f"{rel_error_E*100:.3f}%"
            )
    
    def test_force_balance_at_equilibrium(self):
        """Verify force balance at equilibrium position"""
        alpha = 0.001
        qD = -1.0
        R = 0.5
        qExt = 1.0
        k_coulomb = 138.935456
        
        # Get equilibrium position
        x_eq = self.analytical_equilibrium(alpha, qD, qExt, R)
        
        # Calculate forces at equilibrium
        k_spring = k_coulomb * qD**2 / alpha
        r = R - x_eq
        
        F_spring = -k_spring * x_eq
        F_coulomb = -k_coulomb * qD * qExt / (r**2)
        F_total = F_spring + F_coulomb
        
        # Force should be zero at equilibrium
        assert abs(F_total) < 1e-10, f"Force not balanced at equilibrium: {F_total}"
        
        # Test PyGCMC reaches same equilibrium
        state = pygcmc.MCState()
        state.info.box = [10.0, 10.0, 10.0]
        
        parent = pygcmc.MCAtom()
        parent.x = parent.y = parent.z = 0.0
        parent.charge = 0.0
        
        drude = pygcmc.MCAtom()
        drude.x = 0.01  # Start away from equilibrium
        drude.y = drude.z = 0.0
        drude.charge = qD
        
        ext = pygcmc.MCAtom()
        ext.x = R
        ext.y = ext.z = 0.0
        ext.charge = qExt
        
        state.atoms = [parent, drude, ext]
        state.activeAtomCount = 3
        
        res = pygcmc.MCResidue()
        res.atomStart = 0
        res.atomCount = 2
        res.active = True
        state.residues = [res]
        state.activeResidueCount = 1
        
        pygcmc.DrudeComplete.clear()
        p = pygcmc.DrudeParticle()
        p.drudeIndex = 1
        p.parentIndex = 0
        p.charge = qD
        p.polarizability = alpha
        p.computeSpringConstants()
        pygcmc.DrudeComplete.addParticle(p)
        
        params = pygcmc.DrudeSCFParams()
        params.tolerance = 1e-12
        params.maxIterations = 200
        pygcmc.DrudeComplete.setParameters(params)
        
        pygcmc.DrudeComplete.calculateEnergy(state)
        x_final = state.atoms[1].x
        
        # Should converge to analytical equilibrium
        assert abs(x_final - x_eq) < 1e-8, (
            f"Did not converge to equilibrium: {x_final} vs {x_eq}"
        )
    
    def test_polarizability_scaling(self):
        """Test that response scales correctly with polarizability"""
        qD = -1.0
        R = 0.5
        qExt = 1.0
        
        alpha_values = [0.0005, 0.001, 0.002, 0.004]
        displacements = []
        
        for alpha in alpha_values:
            x_eq = self.analytical_equilibrium(alpha, qD, qExt, R)
            displacements.append(x_eq)
        
        # Displacement should increase with polarizability
        for i in range(1, len(displacements)):
            assert displacements[i] > displacements[i-1], (
                f"Displacement not increasing with α: "
                f"α={alpha_values[i-1]}→{alpha_values[i]}, "
                f"x={displacements[i-1]}→{displacements[i]}"
            )
        
        # Check approximate scaling (for small displacements)
        # x ≈ α * E_field / k_coulomb
        for i in range(len(alpha_values)):
            # Approximate field at origin
            E_field = 138.935456 * qExt / (R**2)
            x_approx = alpha_values[i] * E_field / 138.935456
            
            # Should be reasonable approximation for small α
            if alpha_values[i] < 0.001:
                rel_diff = abs(displacements[i] - x_approx) / displacements[i]
                assert rel_diff < 0.1, (
                    f"Scaling deviation too large for α={alpha_values[i]}: "
                    f"{rel_diff*100:.1f}%"
                )


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
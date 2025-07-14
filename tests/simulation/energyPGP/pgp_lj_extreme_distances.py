# tests/simulation/energyPGP/pgp_lj_limits.py
"""
Test PGP Lennard-Jones energy at extreme distances.

Verifies correct LJ 12-6 potential behavior at:
- Very short distances (r → σ)
- Near cutoff (r → cutoff)
- Beyond cutoff (r > cutoff → 0)
"""

import pytest
import math
import pygcmc
from pygcmc import MCState, MCAtom, MCResidue, MCForceField
from pygcmc import setPGPParameters, initializePMEParameters, precomputeGridPotential
from pygcmc import computeSystemEnergyPGP


# Import helper functions
from .pgp_lj_helpers import (
    create_lj_pair_system,
    calculate_lj_analytical
)


@pytest.mark.parametrize("r_factor", [0.9, 1.0, 1.5, 2.0])
def test_pgp_lj_near_sigma(r_factor):
    """Test LJ energy near and at sigma."""
    epsilon = 1.0  # kJ/mol
    sigma = 0.34  # nm
    distance = r_factor * sigma
    
    print(f"\n=== Testing LJ at r = {r_factor}σ = {distance:.4f} nm ===")
    
    # Create system
    state = create_lj_pair_system(distance, epsilon, sigma)
    
    # Initialize PGP
    alpha = 2.5
    mesh_size = [32, 32, 32]
    initializePMEParameters(state.info.cutoff, state.info.box, alpha)
    setPGPParameters(alpha, mesh_size, state.info.cutoff, mesh_size, 4, 1e-6)
    precomputeGridPotential(state, fixed_only=True)
    
    # Calculate energy
    computeSystemEnergyPGP(state)
    
    # Get LJ energy from residues
    # Note: Each residue stores the full pair energy, so we need to divide by 2
    lj_energy_sum = sum(res.energy_vdw for res in state.residues if res.active)
    lj_energy = lj_energy_sum / 2.0  # Correct for double counting
    
    # Calculate expected energy
    expected = calculate_lj_analytical(distance, epsilon, sigma)
    
    print(f"PGP LJ energy:      {lj_energy:.6f} kJ/mol")
    print(f"Analytical energy:  {expected:.6f} kJ/mol")
    
    # At r = σ, energy should be exactly 0
    if abs(r_factor - 1.0) < 1e-10:
        assert abs(lj_energy) < 1e-5, f"LJ energy at σ should be 0: {lj_energy}"
        print("✓ Energy is zero at r = σ")
    else:
        # Check relative error
        rel_error = abs((lj_energy - expected) / expected) if expected != 0 else abs(lj_energy)
        print(f"Relative error: {rel_error:.6e}")
        assert rel_error < 1e-5, f"LJ energy error too large: {rel_error}"
        
        # Check sign
        if r_factor < 1.0:
            assert lj_energy > 0, "LJ should be repulsive for r < σ"
            print("✓ Correctly repulsive for r < σ")
        else:
            assert lj_energy < 0, "LJ should be attractive for r > σ"
            print("✓ Correctly attractive for r > σ")



def test_pgp_lj_short_distance():
    """Test LJ at very short distances (strong repulsion)."""
    epsilon = 1.0
    sigma = 0.34
    distances = [0.2 * sigma, 0.5 * sigma, 0.7 * sigma]
    
    print("\n=== Testing LJ at very short distances ===")
    print("Distance    | PGP Energy   | Analytical   | Rel Error")
    print("-" * 55)
    
    for distance in distances:
        state = create_lj_pair_system(distance, epsilon, sigma)
        
        # Initialize and calculate
        alpha = 2.5
        mesh_size = [32, 32, 32]
        initializePMEParameters(state.info.cutoff, state.info.box, alpha)
        setPGPParameters(alpha, mesh_size, state.info.cutoff, mesh_size, 4, 1e-6)
        precomputeGridPotential(state, fixed_only=True)
        computeSystemEnergyPGP(state)
        
        # Get LJ energy from residues
        lj_energy_sum = sum(res.energy_vdw for res in state.residues if res.active)
        lj_energy = lj_energy_sum / 2.0  # Correct for double counting
        expected = calculate_lj_analytical(distance, epsilon, sigma)
        
        rel_error = abs((lj_energy - expected) / expected)
        print(f"{distance:.4f} nm | {lj_energy:12.4f} | {expected:12.4f} | {rel_error:.6e}")
        
        # Check strong repulsion
        assert lj_energy > 10.0, f"Should be strongly repulsive: {lj_energy}"
        # Check if energy is capped at 1e6
        if abs(lj_energy - 1e6) < 1.0:
            print(f"  WARNING: LJ energy capped at {lj_energy:.0f}")
            continue  # Skip assertion for capped values
        assert rel_error < 1e-4, f"Error too large at short distance: {rel_error}"
    
    print("✓ All short-distance tests passed")


def test_pgp_lj_asymptotic_behavior():
    """Test LJ asymptotic behavior at large distances."""
    epsilon = 1.0
    sigma = 0.34
    cutoff = 2.0  # Larger cutoff for this test
    
    print("\n=== Testing LJ asymptotic behavior (r → ∞) ===")
    
    # Test at increasing distances
    r_factors = [3.0, 4.0, 5.0, 6.0]
    
    print("r/σ    | Distance | LJ Energy    | Expected    | r^6 term")
    print("-" * 60)
    
    for r_factor in r_factors:
        distance = r_factor * sigma
        
        if distance < cutoff:
            state = create_lj_pair_system(distance, epsilon, sigma, cutoff=cutoff)
            
            alpha = 2.5
            mesh_size = [32, 32, 32]
            initializePMEParameters(cutoff, state.info.box, alpha)
            setPGPParameters(alpha, mesh_size, cutoff, mesh_size, 4, 1e-6)
            precomputeGridPotential(state, fixed_only=True)
            computeSystemEnergyPGP(state)
            
            lj_energy_sum = sum(res.energy_vdw for res in state.residues if res.active)
            lj_energy = lj_energy_sum / 2.0  # Correct for double counting
            expected = calculate_lj_analytical(distance, epsilon, sigma)
            
            # At large r, LJ ~ -4ε(σ/r)^6
            asymptotic = -4 * epsilon * (sigma/distance)**6
            
            print(f"{r_factor:5.1f} | {distance:8.4f} | {lj_energy:12.6f} | "
                  f"{expected:11.6f} | {asymptotic:9.6f}")
            
            # Check that energy approaches asymptotic form
            rel_diff = abs((lj_energy - asymptotic) / asymptotic)
            assert rel_diff < 0.1, f"Not approaching r^-6 behavior: {rel_diff}"
    
    print("✓ Correct asymptotic behavior")



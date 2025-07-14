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
from . import pgp_wrapper
from .pgp_wrapper import setPGPParameters, initializePMEParameters, precomputeGridPotential, computeSystemEnergyPGP

# Import helper functions
from pygcmc import MCState, MCAtom, MCResidue
from pygcmc import MCForceField
from .pgp_lj_helpers import (
    create_lj_pair_system,
    calculate_lj_analytical
)

def test_pgp_lj_minimum():
    """Test LJ energy at minimum (r = 2^(1/6) * σ)."""
    epsilon = 1.0
    sigma = 0.34
    r_min = 2**(1/6) * sigma  # Distance at minimum
    
    print(f"\n=== Testing LJ at minimum r = 2^(1/6)σ = {r_min:.4f} nm ===")
    
    state = create_lj_pair_system(r_min, epsilon, sigma)
    
    # Initialize and calculate
    alpha = 2.5
    mesh_size = [32, 32, 32]
    initializePMEParameters(state.info.cutoff, state.info.box, alpha)
    setPGPParameters(alpha, mesh_size, state.info.cutoff, mesh_size, 4, 1e-6)
    precomputeGridPotential(state)
    computeSystemEnergyPGP(state)
    
    lj_energy_sum = sum(res.energy_vdw for res in state.residues if res.active)
    lj_energy = lj_energy_sum  # PGPContext already distributes energy
    expected = -epsilon  # Minimum energy is -ε
    
    print(f"PGP LJ energy:     {lj_energy:.6f} kJ/mol")
    print(f"Expected minimum:  {expected:.6f} kJ/mol")
    
    rel_error = abs((lj_energy - expected) / expected)
    assert rel_error < 1e-5, f"LJ minimum error: {rel_error}"
    print("✓ Correct energy at LJ minimum")

@pytest.mark.parametrize("cutoff_fraction", [0.8, 0.9, 0.95, 0.99, 1.0, 1.01])
def test_pgp_lj_near_cutoff(cutoff_fraction):
    """Test LJ behavior approaching and beyond cutoff."""
    epsilon = 1.0
    sigma = 0.34
    cutoff = 1.2
    distance = cutoff_fraction * cutoff
    
    print(f"\n=== Testing LJ at {cutoff_fraction:.2f} × cutoff = {distance:.4f} nm ===")
    
    state = create_lj_pair_system(distance, epsilon, sigma, cutoff=cutoff)
    
    # Initialize and calculate
    alpha = 2.5
    mesh_size = [32, 32, 32]
    initializePMEParameters(cutoff, state.info.box, alpha)
    setPGPParameters(alpha, mesh_size, state.info.cutoff, mesh_size, 4, 1e-6)
    precomputeGridPotential(state)
    computeSystemEnergyPGP(state)
    
    lj_energy_sum = sum(res.energy_vdw for res in state.residues if res.active)
    lj_energy = lj_energy_sum  # PGPContext already distributes energy
    
    if distance < cutoff:
        expected = calculate_lj_analytical(distance, epsilon, sigma)
        print(f"PGP LJ energy:     {lj_energy:.6f} kJ/mol")
        print(f"Analytical energy: {expected:.6f} kJ/mol")
        
        rel_error = abs((lj_energy - expected) / expected) if expected != 0 else abs(lj_energy)
        assert rel_error < 1e-5, f"LJ error within cutoff: {rel_error}"
        print("✓ Correct energy within cutoff")
    else:
        print(f"PGP LJ energy: {lj_energy:.6f} kJ/mol")
        assert abs(lj_energy) < 1e-10, f"LJ should be zero beyond cutoff: {lj_energy}"
        print("✓ Correctly zero beyond cutoff")

def test_pgp_lj_mixed_distances():
    """Test system with multiple LJ pairs at different distances."""
    epsilon = 1.0
    sigma = 0.34
    
    print("\n=== Testing multi-particle LJ system ===")
    
    # Create 4-particle system
    state = MCState()
    state.info.box = [6.0, 6.0, 6.0]
    state.info.cutoff = 1.5
    
    ff = MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.ljEps = [epsilon]
    ff.ljSigma = [sigma]
    state.forcefield = ff
    
    # Place 4 atoms in a square
    positions = [
        [3.0, 3.0, 3.0],
        [3.4, 3.0, 3.0],  # 0.4 nm from first
        [3.0, 3.6, 3.0],  # 0.6 nm from first
        [3.4, 3.6, 3.0]   # Various distances
    ]
    
    atoms = []
    for i, pos in enumerate(positions):
        atom = MCAtom()
        atom.x, atom.y, atom.z = pos
        atom.charge = 0.0
        atom.type = 0
        atoms.append(atom)
    
    state.atoms = atoms
    state.activeAtomCount = 4
    
    # Create residues
    residues = []
    for i in range(4):
        res = MCResidue()
        res.active = True
        res.fixed = False
        res.atomStart = i
        res.atomCount = 1
        res.type = 0
        residues.append(res)
    
    state.residues = residues
    state.activeResidueCount = 4
    
    # Initialize and calculate
    alpha = 2.5
    mesh_size = [32, 32, 32]  # Power of 2
    initializePMEParameters(state.info.cutoff, state.info.box, alpha)
    setPGPParameters(alpha, mesh_size, state.info.cutoff, mesh_size, 4, 1e-6)
    precomputeGridPotential(state)
    computeSystemEnergyPGP(state)
    
    # For multi-particle system, sum the residue energies
    total_lj = sum(res.energy_vdw for res in state.residues if res.active)
    
    # Calculate expected energy manually
    expected_total = 0.0
    pair_count = 0
    for i in range(4):
        for j in range(i+1, 4):
            dx = positions[i][0] - positions[j][0]
            dy = positions[i][1] - positions[j][1]
            dz = positions[i][2] - positions[j][2]
            r = math.sqrt(dx*dx + dy*dy + dz*dz)
            if r < state.info.cutoff:
                expected_total += calculate_lj_analytical(r, epsilon, sigma)
                pair_count += 1
                print(f"Pair {i}-{j}: r = {r:.4f} nm, E = {calculate_lj_analytical(r, epsilon, sigma):.6f}")
    
    print(f"\nTotal PGP LJ energy: {total_lj:.6f} kJ/mol")
    print(f"Expected total:      {expected_total:.6f} kJ/mol")
    print(f"Number of pairs within cutoff: {pair_count}")
    
    rel_error = abs((total_lj - expected_total) / expected_total) if expected_total != 0 else abs(total_lj)
    assert rel_error < 1e-4, f"Multi-particle LJ error: {rel_error}"
    print("✓ Multi-particle system correct")

if __name__ == "__main__":
    # Test near sigma
    for r_factor in [0.9, 1.0, 1.5, 2.0]:
        test_pgp_lj_near_sigma(r_factor)
    
    # Test at minimum
    test_pgp_lj_minimum()
    
    # Test near cutoff
    for frac in [0.9, 0.99, 1.01]:
        test_pgp_lj_near_cutoff(frac)
    
    # Test extremes
    test_pgp_lj_short_distance()
    test_pgp_lj_asymptotic_behavior()
    test_pgp_lj_mixed_distances()

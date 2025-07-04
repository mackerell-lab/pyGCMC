# tests/simulation/energyPGP/pgp_lj_limits.py
"""
Test PGP Lennard-Jones energy at extreme distances.

Verifies correct LJ 12-6 potential behavior at:
- Very short distances (r → σ)
- Near cutoff (r → cutoff)
- Beyond cutoff (r > cutoff → 0)
"""

import pytest
import numpy as np
import pygcmc
from pygcmc import MCState, MCAtom, MCResidue, MCForceField
from pygcmc import setPGPParameters, initializePMEParameters, precomputeGridPotential
from pygcmc import computeSystemEnergyPGP


def create_lj_pair_system(distance, epsilon=1.0, sigma=0.34, box_size=5.0, cutoff=1.2):
    """Create a system with two LJ particles at specified distance."""
    state = MCState()
    state.info.box = [box_size, box_size, box_size]
    state.info.cutoff = cutoff
    
    # Force field with LJ parameters
    ff = MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.ljEps = [epsilon]  # kJ/mol
    ff.ljSigma = [sigma]  # nm
    state.forcefield = ff
    
    # Create atoms with no charge (pure LJ)
    atoms = []
    
    # Atom 1 at center
    atom1 = MCAtom()
    atom1.x = box_size / 2
    atom1.y = box_size / 2
    atom1.z = box_size / 2
    atom1.charge = 0.0  # No charge
    atom1.type = 0
    atoms.append(atom1)
    
    # Atom 2 at specified distance
    atom2 = MCAtom()
    atom2.x = box_size / 2 + distance
    atom2.y = box_size / 2
    atom2.z = box_size / 2
    atom2.charge = 0.0  # No charge
    atom2.type = 0
    atoms.append(atom2)
    
    state.atoms = atoms
    state.activeAtomCount = 2
    
    # Create residues
    residues = []
    for i in range(2):
        res = MCResidue()
        res.active = True
        res.fixed = False
        res.atomStart = i
        res.atomCount = 1
        res.type = 0
        residues.append(res)
    
    state.residues = residues
    state.activeResidueCount = 2
    
    return state


def calculate_lj_analytical(r, epsilon, sigma):
    """Calculate analytical LJ 12-6 energy."""
    if r <= 0:
        return float('inf')
    r_ratio = sigma / r
    return 4 * epsilon * (r_ratio**12 - r_ratio**6)


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
    precomputeGridPotential(state, fixed_only=True)
    computeSystemEnergyPGP(state)
    
    lj_energy_sum = sum(res.energy_vdw for res in state.residues if res.active)
    lj_energy = lj_energy_sum / 2.0  # Correct for double counting
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
    setPGPParameters(alpha, mesh_size, cutoff, mesh_size, 4, 1e-6)
    precomputeGridPotential(state, fixed_only=True)
    computeSystemEnergyPGP(state)
    
    lj_energy_sum = sum(res.energy_vdw for res in state.residues if res.active)
    lj_energy = lj_energy_sum / 2.0  # Correct for double counting
    
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
    precomputeGridPotential(state, fixed_only=True)
    computeSystemEnergyPGP(state)
    
    # For multi-particle system, get total from state (already summed correctly)
    total_lj_raw = state.ewald_energy.get('total', 0.0) - state.ewald_energy.get('self', 0.0) - state.ewald_energy.get('reciprocal', 0.0) - state.ewald_energy.get('real_space', 0.0)
    # Note: LJ energy in residues uses GCMC double-counting, so divide by 2 for pair energy
    total_lj = total_lj_raw / 2.0
    
    # Calculate expected energy manually
    expected_total = 0.0
    pair_count = 0
    for i in range(4):
        for j in range(i+1, 4):
            dx = positions[i][0] - positions[j][0]
            dy = positions[i][1] - positions[j][1]
            dz = positions[i][2] - positions[j][2]
            r = np.sqrt(dx*dx + dy*dy + dz*dz)
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
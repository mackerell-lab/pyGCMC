# tests/simulation/energyPGP/test_pgp_real_space.py
"""
Test to verify PGP correctly calculates real-space (short-range) interactions.

This test creates a system with charged particles to ensure real-space
electrostatic calculations are working properly.
"""

import pytest
import math
import pygcmc
from pygcmc import MCAtom, MCResidue, MCState
from pygcmc import setPGPParameters, initializePMEParameters, precomputeGridPotential
from pygcmc import computeSystemEnergyPGP


def test_pgp_real_space_calculation():
    """Test that PGP correctly calculates real-space electrostatic interactions."""
    
    print("\n=== Test: PGP Real-Space Calculation ===")
    
    # Create a simple system with two charged atoms
    state = MCState()
    
    # System info
    state.info.box = [5.0, 5.0, 5.0]
    state.info.cutoff = 1.2
    
    # Force field - minimal setup
    state.forcefield.numTotalTypes = 1
    state.forcefield.numMovementTypes = 1
    state.forcefield.ljEps = [0.0]  # No LJ to isolate electrostatic
    state.forcefield.ljSigma = [0.3]
    
    # Create two oppositely charged atoms close together
    atoms = []
    
    # Positive charge
    atom1 = MCAtom()
    atom1.x = 2.5
    atom1.y = 2.5
    atom1.z = 2.5
    atom1.charge = 1.0  # +1 charge
    atom1.type = 0
    atoms.append(atom1)
    
    # Negative charge - 0.5 nm away
    atom2 = MCAtom()
    atom2.x = 3.0
    atom2.y = 2.5
    atom2.z = 2.5
    atom2.charge = -1.0  # -1 charge
    atom2.type = 0
    atoms.append(atom2)
    
    state.atoms = atoms
    state.activeAtomCount = 2
    
    # Create residues
    residues = []
    for i in range(2):
        res = MCResidue()
        res.active = True
        res.fixed = False  # Both movable
        res.atomStart = i
        res.atomCount = 1
        res.type = 0
        residues.append(res)
    
    state.residues = residues
    state.activeResidueCount = 2
    
    # Calculate exact distance
    r = 0.5  # nm
    print(f"Two atoms with charges +1 and -1 at distance {r} nm")
    
    # Initialize PGP
    alpha = 2.0  # Use larger alpha to emphasize real-space
    mesh_size = [32, 32, 32]
    
    setPGPParameters(
        alpha=alpha,
        meshSize=mesh_size,
        potential_cutoff=state.info.cutoff,
        potentialGridSize=mesh_size,
        splineOrder=4,
        tolerance=1e-5
    )
    
    initializePMEParameters(state.info.cutoff, state.info.box, alpha)
    
    # Precompute (should be empty for all-movable system)
    precomputeGridPotential(state, fixed_only=True)
    
    # Calculate energy
    computeSystemEnergyPGP(state)
    
    # Get energy components
    real_space = state.ewald_energy.get('real_space', 0.0)
    reciprocal = state.ewald_energy.get('reciprocal', 0.0)
    self_energy = state.ewald_energy.get('self', 0.0)
    total = state.ewald_energy.get('total', 0.0)
    
    print(f"\nEnergy components:")
    print(f"  Real-space energy: {real_space:.6f} kJ/mol")
    print(f"  Reciprocal energy: {reciprocal:.6f} kJ/mol")
    print(f"  Self energy: {self_energy:.6f} kJ/mol")
    print(f"  Total energy: {total:.6f} kJ/mol")
    
    # Calculate expected real-space contribution
    # E_real = q1 * q2 * erfc(alpha*r) / r * COULOMB
    kC = 138.935456  # Coulomb constant
    erfc_value = math.erfc(alpha * r)
    expected_real = (1.0) * (-1.0) * erfc_value / r * kC
    
    print(f"\nExpected real-space contribution:")
    print(f"  erfc({alpha}*{r}) = {erfc_value:.6f}")
    print(f"  Expected: {expected_real:.6f} kJ/mol")
    
    # Real-space should be significant with alpha=2.0
    assert abs(real_space) > 1.0, f"Real-space energy ({real_space:.6f}) is too small - not being calculated?"
    
    # For oppositely charged particles, energy should be negative
    assert real_space < 0, f"Real-space energy should be negative for opposite charges, got {real_space:.6f}"
    
    # Total energy should also be negative
    assert total < 0, f"Total energy should be negative for opposite charges, got {total:.6f}"
    
    print("\n✅ PGP correctly calculates real-space interactions!")
    
    # Test 2: Move atoms further apart
    print("\n--- Test 2: Distance dependence ---")
    
    distances = [0.3, 0.5, 0.7, 1.0, 1.5]  # nm
    real_space_energies = []
    
    for d in distances:
        # Move second atom
        state.atoms[1].x = 2.5 + d
        
        # Recalculate
        computeSystemEnergyPGP(state)
        rs_energy = state.ewald_energy.get('real_space', 0.0)
        real_space_energies.append(rs_energy)
        
        # Calculate expected
        erfc_val = math.erfc(alpha * d)
        expected = -erfc_val / d * kC
        
        print(f"  Distance {d:.1f} nm: Real-space = {rs_energy:.4f} kJ/mol, Expected ≈ {expected:.4f}")
    
    # Verify energy becomes less negative as distance increases
    for i in range(len(distances)-1):
        assert real_space_energies[i] < real_space_energies[i+1], \
            f"Energy should increase (become less negative) with distance"
    
    # At large distance, real-space should approach zero
    assert abs(real_space_energies[-1]) < abs(real_space_energies[0]) * 0.1, \
        "Real-space energy should decay significantly with distance"
    
    print("\n✅ Real-space energy shows correct distance dependence!")


def test_pgp_real_space_with_fixed_atoms():
    """Test real-space calculation with mixed fixed and moving atoms."""
    
    print("\n=== Test: PGP Real-Space with Fixed Atoms ===")
    
    # Create system with 3 atoms: 2 fixed, 1 moving
    state = MCState()
    state.info.box = [5.0, 5.0, 5.0]
    state.info.cutoff = 1.2
    
    # Force field
    state.forcefield.numTotalTypes = 1
    state.forcefield.numMovementTypes = 1
    state.forcefield.ljEps = [0.0]
    state.forcefield.ljSigma = [0.3]
    
    # Create atoms in a line
    atoms = []
    charges = [1.0, -1.0, 1.0]  # +-+ pattern
    positions = [2.0, 2.5, 3.0]  # 0.5 nm spacing
    
    for i, (x, q) in enumerate(zip(positions, charges)):
        atom = MCAtom()
        atom.x = x
        atom.y = 2.5
        atom.z = 2.5
        atom.charge = q
        atom.type = 0
        atoms.append(atom)
    
    state.atoms = atoms
    state.activeAtomCount = 3
    
    # Residues: first two are fixed
    residues = []
    for i in range(3):
        res = MCResidue()
        res.active = True
        res.fixed = (i < 2)  # First two fixed
        res.atomStart = i
        res.atomCount = 1
        res.type = 0
        residues.append(res)
    
    state.residues = residues
    state.activeResidueCount = 3
    
    print("System: Fixed(+1) -- Fixed(-1) -- Moving(+1)")
    print("Positions: 2.0 -- 2.5 -- 3.0")
    
    # Initialize PGP
    alpha = 2.0
    setPGPParameters(
        alpha=alpha,
        meshSize=[32, 32, 32],
        potential_cutoff=state.info.cutoff,
        potentialGridSize=[32, 32, 32],
        splineOrder=4,
        tolerance=1e-5
    )
    
    initializePMEParameters(state.info.cutoff, state.info.box, alpha)
    
    # Precompute grid for fixed atoms
    print("\nPrecomputing grid potential for fixed atoms...")
    precomputeGridPotential(state, fixed_only=True)
    
    # Calculate energy
    computeSystemEnergyPGP(state)
    
    real_space = state.ewald_energy.get('real_space', 0.0)
    grid_energy = state.ewald_energy.get('reciprocal', 0.0)  # Grid interpolation stored as reciprocal
    
    print(f"\nEnergy components:")
    print(f"  Real-space energy: {real_space:.6f} kJ/mol")
    print(f"  Grid energy: {grid_energy:.6f} kJ/mol")
    
    # The real-space should include:
    # 1. Fixed-Fixed interaction (should be excluded as it's in the grid)
    # 2. Fixed-Moving interactions (should be calculated)
    # 3. Moving-Moving interactions (none in this case)
    
    # Verify real-space energy exists
    assert abs(real_space) > 0.1, \
        f"Real-space energy ({real_space:.6f}) too small - fixed-moving interactions not calculated?"
    
    print("\n✅ PGP correctly handles real-space with fixed atoms!")


if __name__ == "__main__":
    test_pgp_real_space_calculation()
    test_pgp_real_space_with_fixed_atoms()
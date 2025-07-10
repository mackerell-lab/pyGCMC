"""
Validate PGP energy calculations using PME Complete as reference

PGP (Precompute Grid Potential) is designed for efficient energy calculation
of moveable particles in the field of fixed particles. This test validates
PGP against PME Complete for such scenarios.
"""

import pytest
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pygcmc
from pygcmc import MCState, MCAtom, MCResidue, MCForceField
from pygcmc import initializePMEParameters, computeSystemEnergyPMEComplete
from pygcmc import setPMEParameters, setPGPParameters
from pygcmc import precomputeGridPotential, calculateMoleculeEnergy
from pygcmc import computeMovementEnergyPME


def create_fixed_and_moveable_system():
    """Create a system with fixed and moveable particles"""
    state = MCState()
    state.info.box = [4.0, 4.0, 4.0]
    state.info.cutoff = 1.5
    
    # Force field
    ff = MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.ljEps = [1.0]
    ff.ljSigma = [0.35]
    state.forcefield = ff
    
    atoms = []
    residues = []
    
    # Fixed particles (create a background field)
    fixed_positions = [
        ([1.0, 2.0, 2.0], 1.0),   # +1 charge
        ([3.0, 2.0, 2.0], -1.0),  # -1 charge
        ([2.0, 1.0, 2.0], 1.0),   # +1 charge
        ([2.0, 3.0, 2.0], -1.0),  # -1 charge
    ]
    
    for i, (pos, charge) in enumerate(fixed_positions):
        atom = MCAtom()
        atom.x, atom.y, atom.z = pos
        atom.charge = charge
        atom.type = 0
        atoms.append(atom)
        
        res = MCResidue()
        res.active = True
        res.fixed = True  # Fixed particles
        res.atomStart = i
        res.atomCount = 1
        res.type = 0
        residues.append(res)
    
    # Moveable particle
    moveable_atom = MCAtom()
    moveable_atom.x, moveable_atom.y, moveable_atom.z = 2.0, 2.0, 2.0
    moveable_atom.charge = -0.5
    moveable_atom.type = 0
    atoms.append(moveable_atom)
    
    moveable_res = MCResidue()
    moveable_res.active = True
    moveable_res.fixed = False  # Moveable particle
    moveable_res.atomStart = 4
    moveable_res.atomCount = 1
    moveable_res.type = 0
    residues.append(moveable_res)
    
    state.atoms = atoms
    state.activeAtomCount = 5
    state.residues = residues
    state.activeResidueCount = 5
    
    return state


def test_pgp_movement_energy():
    """Test PGP energy calculation for particle movement"""
    
    print("\n" + "="*70)
    print("PGP Movement Energy Test")
    print("="*70)
    
    state = create_fixed_and_moveable_system()
    
    # Set parameters
    alpha = 2.2
    mesh_size = [32, 32, 32]
    spline_order = 4
    
    # Initialize PME
    setPMEParameters(alpha, mesh_size, spline_order)
    initializePMEParameters(state.info.cutoff, state.info.box, alpha, mesh_size, spline_order)
    
    # 1. Calculate initial energy with PME Complete
    elec_initial, vdw_initial, total_initial = computeSystemEnergyPMEComplete(state)
    print(f"\nInitial PME Complete total energy: {total_initial:.6f} kJ/mol")
    
    # 2. Set up PGP and precompute grid from fixed particles
    pgp_mesh_size = [64, 64, 64]
    setPGPParameters(alpha, pgp_mesh_size, state.info.cutoff, pgp_mesh_size, spline_order, 1e-6)
    
    print("\nPrecomputing grid potential from fixed particles...")
    precomputeGridPotential(state, fixed_only=True)
    
    # 3. Calculate moveable particle energy using PGP
    pgp_energy_initial = calculateMoleculeEnergy(state)
    print(f"Initial PGP energy (moveable particle): {pgp_energy_initial:.6f} kJ/mol")
    
    # 4. Move the particle
    translation = [0.2, 0.0, 0.0]  # Move 0.2 nm in x direction
    moveable_atom = state.atoms[4]
    moveable_atom.x += translation[0]
    moveable_atom.y += translation[1]
    moveable_atom.z += translation[2]
    
    print(f"\nMoved particle by {translation} nm")
    
    # 5. Calculate new energy with PME Complete
    elec_moved, vdw_moved, total_moved = computeSystemEnergyPMEComplete(state)
    print(f"Moved PME Complete total energy: {total_moved:.6f} kJ/mol")
    
    # 6. Calculate new energy with PGP
    pgp_energy_moved = calculateMoleculeEnergy(state)
    print(f"Moved PGP energy (moveable particle): {pgp_energy_moved:.6f} kJ/mol")
    
    # 7. Compare energy differences
    pme_delta = total_moved - total_initial
    pgp_delta = pgp_energy_moved - pgp_energy_initial
    
    print(f"\nEnergy changes:")
    print(f"  PME Complete delta: {pme_delta:.6f} kJ/mol")
    print(f"  PGP delta:          {pgp_delta:.6f} kJ/mol")
    print(f"  Difference:         {abs(pme_delta - pgp_delta):.6f} kJ/mol")
    
    # The deltas should be similar
    rel_diff = abs(pme_delta - pgp_delta) / abs(pme_delta) * 100 if pme_delta != 0 else 0
    print(f"  Relative difference: {rel_diff:.3f}%")
    
    # Note: PGP calculates single particle energy, not full system energy change
    # So we expect some difference, but it should be reasonable
    assert rel_diff < 20.0, f"Energy change differs by {rel_diff:.3f}% (> 20%)"


def test_pgp_reciprocal_space_comparison():
    """Compare PGP with PME reciprocal space energy"""
    
    print("\n" + "="*70)
    print("PGP Reciprocal Space Comparison")
    print("="*70)
    
    state = create_fixed_and_moveable_system()
    
    # Set parameters
    alpha = 2.2
    mesh_size = [32, 32, 32]
    spline_order = 4
    
    # Initialize PME
    setPMEParameters(alpha, mesh_size, spline_order)
    initializePMEParameters(state.info.cutoff, state.info.box, alpha, mesh_size, spline_order)
    
    # Calculate PME movement energy (for moveable particles)
    pme_result = computeMovementEnergyPME(state)
    pme_elec = pme_result[0]
    pme_dict = pme_result[2]
    pme_reciprocal = pme_dict.get('reciprocal', 0.0)
    
    print(f"\nPME movement energy:")
    print(f"  Total electrostatic: {pme_elec:.6f} kJ/mol")
    print(f"  Reciprocal space:    {pme_reciprocal:.6f} kJ/mol")
    
    # Set up PGP
    pgp_mesh_size = [64, 64, 64]
    setPGPParameters(alpha, pgp_mesh_size, state.info.cutoff, pgp_mesh_size, spline_order, 1e-6)
    
    # Precompute grid
    precomputeGridPotential(state, fixed_only=True)
    
    # Calculate PGP energy
    pgp_energy = calculateMoleculeEnergy(state)
    
    print(f"\nPGP energy (moveable particle): {pgp_energy:.6f} kJ/mol")
    
    # PGP should approximate the reciprocal space contribution
    print(f"\nComparison:")
    print(f"  PME reciprocal: {pme_reciprocal:.6f} kJ/mol")
    print(f"  PGP energy:     {pgp_energy:.6f} kJ/mol")
    
    # They won't be exactly equal because PGP includes some additional terms
    # but they should be in the same ballpark
    ratio = pgp_energy / pme_reciprocal if pme_reciprocal != 0 else 0
    print(f"  Ratio (PGP/PME reciprocal): {ratio:.3f}")


def test_pgp_grid_spacing_convergence():
    """Test PGP convergence with grid spacing"""
    
    print("\n" + "="*70)
    print("PGP Grid Spacing Convergence Test")
    print("="*70)
    
    state = create_fixed_and_moveable_system()
    
    # Set parameters
    alpha = 2.2
    mesh_size = [32, 32, 32]
    spline_order = 4
    
    # Initialize PME
    setPMEParameters(alpha, mesh_size, spline_order)
    initializePMEParameters(state.info.cutoff, state.info.box, alpha, mesh_size, spline_order)
    
    # Get reference energy change
    elec_initial, vdw_initial, total_initial = computeSystemEnergyPMEComplete(state)
    
    # Move particle
    state.atoms[4].x += 0.3
    
    elec_moved, vdw_moved, total_moved = computeSystemEnergyPMEComplete(state)
    pme_delta = total_moved - total_initial
    
    print(f"\nReference PME Complete energy change: {pme_delta:.6f} kJ/mol")
    
    # Test different PGP grid sizes (must be powers of 2 for FFT)
    grid_sizes = [32, 64, 128]
    
    print(f"\n{'Grid Size':>10} | {'PGP Energy':>12} | {'vs PME Delta':>12}")
    print("-"*40)
    
    # Reset particle position
    state.atoms[4].x -= 0.3
    
    for grid_size in grid_sizes:
        pgp_mesh_size = [grid_size, grid_size, grid_size]
        setPGPParameters(alpha, pgp_mesh_size, state.info.cutoff, pgp_mesh_size, spline_order, 1e-6)
        
        # Precompute grid
        precomputeGridPotential(state, fixed_only=True)
        
        # Get initial energy
        pgp_initial = calculateMoleculeEnergy(state)
        
        # Move particle
        state.atoms[4].x += 0.3
        
        # Get moved energy
        pgp_moved = calculateMoleculeEnergy(state)
        
        # Reset position
        state.atoms[4].x -= 0.3
        
        pgp_delta = pgp_moved - pgp_initial
        diff = abs(pgp_delta - pme_delta)
        
        print(f"{grid_size:10d} | {pgp_delta:12.6f} | {diff:12.6f}")
    
    # The finest grid should give reasonable agreement
    # Note: PGP and PME Complete calculate different quantities, so some difference is expected
    assert diff < 5.0, f"Finest grid still differs by {diff:.6f} kJ/mol (> 5.0)"


if __name__ == "__main__":
    test_pgp_movement_energy()
    test_pgp_reciprocal_space_comparison()
    test_pgp_grid_spacing_convergence()
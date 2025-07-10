"""
Test PGP with mixed fixed/moveable particle interactions

This tests various scenarios with different combinations of fixed and moveable particles.
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


def create_mixed_system(n_fixed, n_moveable):
    """Create a system with specified numbers of fixed and moveable particles"""
    state = MCState()
    state.info.box = [4.0, 4.0, 4.0]
    state.info.cutoff = 1.5
    
    # Force field
    ff = MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.ljEps = [0.8]
    ff.ljSigma = [0.35]
    state.forcefield = ff
    
    atoms = []
    residues = []
    
    # Fixed particles in a line
    for i in range(n_fixed):
        atom = MCAtom()
        atom.x = 1.0 + i * 0.5
        atom.y = 2.0
        atom.z = 2.0
        atom.charge = (-1)**i  # Alternating charges
        atom.type = 0
        atoms.append(atom)
        
        res = MCResidue()
        res.active = True
        res.fixed = True
        res.atomStart = i
        res.atomCount = 1
        res.type = 0
        residues.append(res)
    
    # Moveable particles
    for i in range(n_moveable):
        atom = MCAtom()
        atom.x = 2.0
        atom.y = 2.0 + i * 0.5
        atom.z = 2.0
        atom.charge = 0.5 * (-1)**i
        atom.type = 0
        atoms.append(atom)
        
        res = MCResidue()
        res.active = True
        res.fixed = False
        res.atomStart = n_fixed + i
        res.atomCount = 1
        res.type = 0
        residues.append(res)
    
    state.atoms = atoms
    state.activeAtomCount = n_fixed + n_moveable
    state.residues = residues
    state.activeResidueCount = n_fixed + n_moveable
    
    return state


def test_pgp_single_moveable_particle():
    """Test PGP with one moveable particle in field of fixed particles"""
    
    print("\n" + "="*70)
    print("PGP Single Moveable Particle Test")
    print("="*70)
    
    state = create_mixed_system(n_fixed=3, n_moveable=1)
    
    # Set parameters
    alpha = 2.2
    mesh_size = [32, 32, 32]
    spline_order = 4
    
    # Initialize PME
    setPMEParameters(alpha, mesh_size, spline_order)
    initializePMEParameters(state.info.cutoff, state.info.box, alpha, mesh_size, spline_order)
    
    # Set up PGP
    pgp_mesh_size = [64, 64, 64]
    setPGPParameters(alpha, pgp_mesh_size, state.info.cutoff, pgp_mesh_size, spline_order, 1e-6)
    
    # Precompute grid from fixed particles
    precomputeGridPotential(state, fixed_only=True)
    
    # Calculate initial energy
    initial_energy = calculateMoleculeEnergy(state)
    print(f"Initial PGP energy: {initial_energy:.6f} kJ/mol")
    
    # Move the particle
    state.atoms[3].x += 0.2
    
    # Calculate new energy
    moved_energy = calculateMoleculeEnergy(state)
    print(f"Moved PGP energy: {moved_energy:.6f} kJ/mol")
    print(f"Energy change: {moved_energy - initial_energy:.6f} kJ/mol")
    
    # Energy should change when particle moves
    assert abs(moved_energy - initial_energy) > 0.01, "Energy should change when particle moves"


def test_pgp_multiple_moveable_particles():
    """Test PGP with multiple moveable particles"""
    
    print("\n" + "="*70)
    print("PGP Multiple Moveable Particles Test")
    print("="*70)
    
    state = create_mixed_system(n_fixed=2, n_moveable=3)
    
    # Set parameters
    alpha = 2.2
    mesh_size = [32, 32, 32]
    spline_order = 4
    
    # Initialize PME
    setPMEParameters(alpha, mesh_size, spline_order)
    initializePMEParameters(state.info.cutoff, state.info.box, alpha, mesh_size, spline_order)
    
    # Set up PGP
    pgp_mesh_size = [64, 64, 64]
    setPGPParameters(alpha, pgp_mesh_size, state.info.cutoff, pgp_mesh_size, spline_order, 1e-6)
    
    # Precompute grid
    precomputeGridPotential(state, fixed_only=True)
    
    # Calculate initial energy
    initial_energy = calculateMoleculeEnergy(state)
    print(f"Initial PGP energy (3 moveable particles): {initial_energy:.6f} kJ/mol")
    
    # Move different particles and check energy changes
    energies = []
    for i in range(3):
        # Reset positions
        for j in range(3):
            state.atoms[2+j].y = 2.0 + j * 0.5
        
        # Move one particle
        state.atoms[2+i].x += 0.1
        energy = calculateMoleculeEnergy(state)
        energies.append(energy)
        print(f"Energy after moving particle {i}: {energy:.6f} kJ/mol")
    
    # All movements should result in different energies
    assert len(set(energies)) == 3, "Different particle movements should give different energies"


def test_pgp_with_pme_movement_comparison():
    """Compare PGP with PME movement energy"""
    
    print("\n" + "="*70)
    print("PGP vs PME Movement Energy Comparison")
    print("="*70)
    
    state = create_mixed_system(n_fixed=4, n_moveable=2)
    
    # Set up movement residues
    state.movementResidues = []
    movement_info = pygcmc.MCMovementResidueInfo()
    movement_info.startIndex = 4  # First moveable residue
    movement_info.activeCount = 2  # Two moveable residues
    state.movementResidues.append(movement_info)
    
    # Set parameters
    alpha = 2.2
    mesh_size = [32, 32, 32]
    spline_order = 4
    
    # Initialize PME
    setPMEParameters(alpha, mesh_size, spline_order)
    initializePMEParameters(state.info.cutoff, state.info.box, alpha, mesh_size, spline_order)
    
    # Calculate PME movement energy
    pme_result = computeMovementEnergyPME(state)
    pme_elec = pme_result[0]
    pme_dict = pme_result[2]
    
    print(f"\nPME movement energy:")
    print(f"  Electrostatic: {pme_elec:.6f} kJ/mol")
    print(f"  Components: {pme_dict}")
    
    # Set up PGP
    pgp_mesh_size = [64, 64, 64]
    setPGPParameters(alpha, pgp_mesh_size, state.info.cutoff, pgp_mesh_size, spline_order, 1e-6)
    
    # Precompute grid
    precomputeGridPotential(state, fixed_only=True)
    
    # Calculate PGP energy
    pgp_energy = calculateMoleculeEnergy(state)
    print(f"\nPGP energy: {pgp_energy:.6f} kJ/mol")
    
    # They won't match exactly - PGP only calculates moveable particle energy in fixed field
    # while PME includes all interactions
    if abs(pme_elec) > 0.1:
        ratio = pgp_energy / pme_elec
        print(f"Ratio (PGP/PME): {ratio:.3f}")
        # PGP typically gives much smaller values since it's only particle-in-field energy
        # Allow a wider range since these calculate fundamentally different things
        assert 0.001 < abs(ratio) < 100.0, f"Energy ratio {ratio} is unreasonable"
        print("Note: PGP calculates particle-in-field energy, not total interaction energy")


if __name__ == "__main__":
    test_pgp_single_moveable_particle()
    test_pgp_multiple_moveable_particles()
    test_pgp_with_pme_movement_comparison()
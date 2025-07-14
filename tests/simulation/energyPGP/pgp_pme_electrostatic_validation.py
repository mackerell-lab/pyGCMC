"""
Test PGP vs PME Complete for pure electrostatic systems

For pure charge systems (no LJ), the energy change (dE) calculated by PGP 
and PME Complete should be identical, as both calculate the same electrostatic
interactions without the complications of LJ terms.
"""

import pytest
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pygcmc
from . import pgp_wrapper
from .pgp_wrapper import setPMEParameters, initializePMEParameters, setPGPParameters
from .pgp_wrapper import precomputeGridPotential, calculateMoleculeEnergy, computeMovementEnergyPME
from .pgp_wrapper import resetPGPState
from pygcmc import MCState, MCAtom, MCResidue
from pygcmc import MCForceField, MCMovementResidueInfo

def create_pure_electrostatic_system():
    """Create a system with only electrostatic interactions (no LJ)"""
    state = MCState()
    state.info.box = [5.0, 5.0, 5.0]
    state.info.cutoff = 2.0
    
    # Force field with zero LJ parameters
    ff = MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.ljEps = [0.0]  # Zero epsilon - no LJ interaction
    ff.ljSigma = [0.0]  # Zero sigma
    state.forcefield = ff
    
    atoms = []
    residues = []
    
    # Fixed charged particles in a square arrangement
    fixed_positions = [
        ([2.0, 2.0, 2.0], 1.0),   # Center positive
        ([1.0, 2.0, 2.0], -0.5),  # Left negative
        ([3.0, 2.0, 2.0], -0.5),  # Right negative
        ([2.0, 1.0, 2.0], -0.5),  # Bottom negative
        ([2.0, 3.0, 2.0], -0.5),  # Top negative
    ]
    
    for i, (pos, charge) in enumerate(fixed_positions):
        atom = MCAtom()
        atom.x, atom.y, atom.z = pos
        atom.charge = charge
        atom.type = 0
        atoms.append(atom)
        
        res = MCResidue()
        res.active = True
        res.fixed = True
        res.atomStart = i
        res.atomCount = 1
        res.type = 0
        residues.append(res)
    
    # Moveable charged particle
    moveable_atom = MCAtom()
    moveable_atom.x, moveable_atom.y, moveable_atom.z = 2.5, 2.5, 2.5
    moveable_atom.charge = 0.8
    moveable_atom.type = 0
    atoms.append(moveable_atom)
    
    moveable_res = MCResidue()
    moveable_res.active = True
    moveable_res.fixed = False
    moveable_res.atomStart = 5
    moveable_res.atomCount = 1
    moveable_res.type = 0
    residues.append(moveable_res)
    
    state.atoms = atoms
    state.activeAtomCount = 6
    state.residues = residues
    state.activeResidueCount = 6
    
    return state

def test_pgp_pme_pure_electrostatic_exact():
    """Test that PGP and PME movement energy give identical dE for pure electrostatic systems"""
    
    print("\n" + "="*70)
    print("PGP vs PME Movement Energy - Pure Electrostatic System Test")
    print("="*70)
    
    state = create_pure_electrostatic_system()
    
    # Set up movement residues for PME movement energy calculation
    state.movementResidues = []
    movement_info = pygcmc.MCMovementResidueInfo()
    movement_info.startIndex = 5  # The moveable residue
    movement_info.activeCount = 1
    state.movementResidues.append(movement_info)
    
    # Use identical parameters for both methods
    alpha = 2.2
    mesh_size = [64, 64, 64]  # Use same mesh size for both
    spline_order = 4
    
    # Initialize PME
    setPMEParameters(alpha, mesh_size, spline_order)
    initializePMEParameters(state.info.cutoff, state.info.box, alpha)
    
    # Calculate initial PME movement energy
    pme_move_init = computeMovementEnergyPME(state)
    pme_move_elec_init = pme_move_init[0]
    print(f"\nInitial PME movement energy: {pme_move_elec_init:.8f} kJ/mol")
    
    # Set up PGP with identical parameters
    setPGPParameters(alpha, mesh_size, state.info.cutoff, mesh_size, 4, 1e-6)
    
    # Precompute grid from fixed particles
    print("\nPrecomputing PGP grid from fixed particles...")
    precomputeGridPotential(state, fixed_only=True)
    
    # Calculate initial PGP energy
    pgp_initial = calculateMoleculeEnergy(state)
    print(f"Initial PGP energy: {pgp_initial:.8f} kJ/mol")
    
    # Move the particle
    move_distance = 0.3
    state.atoms[5].x += move_distance
    state.atoms[5].y += move_distance
    state.atoms[5].z += move_distance
    
    print(f"\nMoved particle by ({move_distance}, {move_distance}, {move_distance}) nm")
    
    # Calculate new PME movement energy
    pme_move_final = computeMovementEnergyPME(state)
    pme_move_elec_final = pme_move_final[0]
    print(f"\nMoved PME movement energy: {pme_move_elec_final:.8f} kJ/mol")
    
    # Calculate new PGP energy
    pgp_moved = calculateMoleculeEnergy(state)
    print(f"Moved PGP energy: {pgp_moved:.8f} kJ/mol")
    
    # Calculate energy changes
    pme_delta = pme_move_elec_final - pme_move_elec_init
    pgp_delta = pgp_moved - pgp_initial
    
    print(f"\nEnergy changes (dE):")
    print(f"  PME movement: {pme_delta:.8f} kJ/mol")
    print(f"  PGP:          {pgp_delta:.8f} kJ/mol")
    print(f"  Difference:   {abs(pme_delta - pgp_delta):.8f} kJ/mol")
    
    # For pure electrostatic systems, these should be nearly identical
    rel_error = abs(pme_delta - pgp_delta) / abs(pme_delta) * 100 if pme_delta != 0 else 0
    print(f"  Relative error: {rel_error:.4f}%")
    
    # Very strict tolerance when comparing with PME movement energy
    assert rel_error < 0.1, f"For pure electrostatic system, PGP and PME movement energy should give nearly identical dE. Got {rel_error:.4f}% difference"

def test_pgp_pme_different_movements():
    """Test multiple particle movements in pure electrostatic system"""
    
    print("\n" + "="*70)
    print("PGP vs PME Movement Energy - Multiple Movement Test")
    print("="*70)
    
    state = create_pure_electrostatic_system()
    
    # Set up movement residues
    state.movementResidues = []
    movement_info = pygcmc.MCMovementResidueInfo()
    movement_info.startIndex = 5
    movement_info.activeCount = 1
    state.movementResidues.append(movement_info)
    
    # Use identical parameters
    alpha = 2.2
    mesh_size = [64, 64, 64]
    spline_order = 4
    
    # Initialize both methods
    setPMEParameters(alpha, mesh_size, spline_order)
    initializePMEParameters(state.info.cutoff, state.info.box, alpha)
    setPGPParameters(alpha, mesh_size, state.info.cutoff, mesh_size, 4, 1e-6)
    
    # Precompute PGP grid
    precomputeGridPotential(state, fixed_only=True)
    
    # Test different movements
    movements = [
        (0.1, 0.0, 0.0),   # Move in x only
        (0.0, 0.1, 0.0),   # Move in y only
        (0.0, 0.0, 0.1),   # Move in z only
        (0.1, 0.1, 0.1),   # Move diagonally
        (-0.2, 0.0, 0.1),  # Mixed movement
    ]
    
    print(f"\n{'Movement':>15} | {'PME dE':>12} | {'PGP dE':>12} | {'Rel Error':>10}")
    print("-" * 55)
    
    for dx, dy, dz in movements:
        # Get initial energies
        pme_init = computeMovementEnergyPME(state)[0]
        pgp_init = calculateMoleculeEnergy(state)
        
        # Move particle
        state.atoms[5].x += dx
        state.atoms[5].y += dy
        state.atoms[5].z += dz
        
        # Get final energies
        pme_final = computeMovementEnergyPME(state)[0]
        pgp_final = calculateMoleculeEnergy(state)
        
        # Calculate dE
        pme_de = pme_final - pme_init
        pgp_de = pgp_final - pgp_init
        
        # Calculate relative error
        rel_err = abs(pme_de - pgp_de) / abs(pme_de) * 100 if pme_de != 0 else 0
        
        print(f"({dx:4.1f},{dy:4.1f},{dz:4.1f}) | {pme_de:12.6f} | {pgp_de:12.6f} | {rel_err:9.4f}%")
        
        # Reset position
        state.atoms[5].x -= dx
        state.atoms[5].y -= dy
        state.atoms[5].z -= dz
        
        # Should have very low error when comparing with PME movement energy
        assert rel_err < 0.1, f"Movement ({dx},{dy},{dz}) has {rel_err:.4f}% error (> 0.1%)"

def test_pgp_pme_convergence_with_mesh():
    """Test PGP-PME convergence as mesh size increases"""
    
    print("\n" + "="*70)
    print("PGP vs PME Movement Energy - Mesh Convergence Test")
    print("="*70)
    
    state = create_pure_electrostatic_system()
    
    # Set up movement residues
    state.movementResidues = []
    movement_info = MCMovementResidueInfo()
    movement_info.startIndex = 5
    movement_info.activeCount = 1
    state.movementResidues.append(movement_info)
    
    # Fixed parameters
    alpha = 2.2
    spline_order = 4
    
    # Test different mesh sizes
    mesh_sizes = [32, 64, 128]
    
    print(f"\n{'Mesh Size':>10} | {'PME dE':>12} | {'PGP dE':>12} | {'Rel Error':>10}")
    print("-" * 50)
    
    for mesh_size in mesh_sizes:
        mesh = [mesh_size, mesh_size, mesh_size]
        
        # Reset PGP state to avoid conflicts
        pgp_wrapper.resetPGPState()
        
        # Initialize both methods with same mesh
        setPMEParameters(alpha, mesh, spline_order)
        initializePMEParameters(state.info.cutoff, state.info.box, alpha)
        setPGPParameters(alpha, mesh, state.info.cutoff, mesh, 4, 1e-6)
        
        # Precompute PGP grid
        precomputeGridPotential(state, fixed_only=True)
        
        # Get initial energies
        pme_init = computeMovementEnergyPME(state)[0]
        pgp_init = calculateMoleculeEnergy(state)
        
        # Move particle
        state.atoms[5].x += 0.2
        state.atoms[5].y += 0.2
        
        # Get final energies
        pme_final = computeMovementEnergyPME(state)[0]
        pgp_final = calculateMoleculeEnergy(state)
        
        # Reset position
        state.atoms[5].x -= 0.2
        state.atoms[5].y -= 0.2
        
        # Calculate dE
        pme_de = pme_final - pme_init
        pgp_de = pgp_final - pgp_init
        
        # Calculate relative error
        rel_err = abs(pme_de - pgp_de) / abs(pme_de) * 100 if pme_de != 0 else 0
        
        print(f"{mesh_size:10d} | {pme_de:12.6f} | {pgp_de:12.6f} | {rel_err:9.4f}%")
        
        # Higher mesh sizes should have very low error
        if mesh_size >= 64:
            assert rel_err < 0.05, f"Mesh size {mesh_size} has {rel_err:.4f}% error (> 0.05%)"

if __name__ == "__main__":
    test_pgp_pme_pure_electrostatic_exact()
    test_pgp_pme_different_movements()
    test_pgp_pme_convergence_with_mesh()

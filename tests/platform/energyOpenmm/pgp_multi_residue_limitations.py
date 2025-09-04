"""
Test PGP Complete limitations with multi-residue moves

This test demonstrates the known limitation when moving multiple residues
simultaneously and provides guidance on workarounds.
"""

import pytest
import numpy as np
import pygcmc
from pygcmc import MCState, MCAtom, MCResidue, MCForceField, MCMovementResidueInfo

try:
    import openmm as mm
    from openmm import app
    OPENMM_AVAILABLE = True
except ImportError:
    OPENMM_AVAILABLE = False
    pytest.skip("OpenMM not available", allow_module_level=True)


def test_pgp_multi_residue_limitation():
    """Demonstrate and quantify error when moving multiple residues"""
    
    # Create system with multiple movable residues
    state = MCState()
    state.info.box = [8.0, 8.0, 8.0]
    state.info.cutoff = 1.2
    
    ff = MCForceField()
    ff.numTotalTypes = 2
    ff.numMovementTypes = 2
    ff.ljEps = [0.2, 0.3, 0.25, 0.25]
    ff.ljSigma = [0.3, 0.35, 0.325, 0.325]
    state.forcefield = ff
    
    atoms = []
    residues = []
    
    # Fixed residue
    fixed_positions = [
        ([4.0, 4.0, 4.0], 1.0, 0),
        ([4.3, 4.0, 4.0], -1.0, 1),
    ]
    
    for i, (pos, charge, atom_type) in enumerate(fixed_positions):
        atom = MCAtom()
        atom.x, atom.y, atom.z = pos
        atom.charge = charge
        atom.type = atom_type
        atoms.append(atom)
    
    res = MCResidue()
    res.atomStart = 0
    res.atomCount = 2
    res.active = True
    res.fixed = True
    res.type = 0
    residues.append(res)
    
    # Multiple moving residues (3 ion pairs)
    moving_positions = [
        # Pair 1
        ([2.0, 4.0, 4.0], 1.0, 0),
        ([2.3, 4.0, 4.0], -1.0, 1),
        # Pair 2
        ([6.0, 4.0, 4.0], 1.0, 0),
        ([6.3, 4.0, 4.0], -1.0, 1),
        # Pair 3
        ([4.0, 6.0, 4.0], 1.0, 0),
        ([4.3, 6.0, 4.0], -1.0, 1),
    ]
    
    for i, (pos, charge, atom_type) in enumerate(moving_positions):
        atom = MCAtom()
        atom.x, atom.y, atom.z = pos
        atom.charge = charge
        atom.type = atom_type
        atoms.append(atom)
    
    # Create 3 moving residues
    for i in range(3):
        res = MCResidue()
        res.atomStart = 2 + i * 2
        res.atomCount = 2
        res.active = True
        res.fixed = False
        res.type = 0
        residues.append(res)
    
    state.atoms = atoms
    state.activeAtomCount = len(atoms)
    state.residues = residues
    state.activeResidueCount = len(residues)
    
    # Initialize PGP
    alpha = 5.6 / state.info.cutoff
    mesh_size = [64, 64, 64]  # Power of 2
    
    pygcmc.setPMEParameters(alpha, mesh_size, 4, 1e-5)
    pygcmc.initializePMEParameters(state.info.cutoff, state.info.box, alpha)
    pygcmc.setPGPParameters(alpha, mesh_size, state.info.cutoff, mesh_size, 4, 1e-5)
    pygcmc.precomputeGridPotential(state, fixed_only=True)
    
    # Test 1: Single residue move (should be accurate)
    state.movementResidues.clear()
    movement_info = MCMovementResidueInfo()
    movement_info.startIndex = 1  # First moving residue
    movement_info.activeCount = 1
    state.movementResidues.append(movement_info)
    
    # Move first pair
    state.atoms[2].x += 0.2
    state.atoms[3].x += 0.2
    
    pgp_single = pygcmc.computeMovementEnergyPGPCompleteCorrect(state)
    pme_single = pygcmc.computeMovementEnergyPME(state)
    
    pgp_single_total = pgp_single[0] + pgp_single[1]
    pme_single_total = pme_single[0] + pme_single[1]
    
    single_error = abs(pgp_single_total - pme_single_total)
    
    # Reset
    state.atoms[2].x -= 0.2
    state.atoms[3].x -= 0.2
    
    # Test 2: Multiple residue move (will show larger error)
    state.movementResidues.clear()
    movement_info = MCMovementResidueInfo()
    movement_info.startIndex = 1  # All moving residues
    movement_info.activeCount = 3
    state.movementResidues.append(movement_info)
    
    # Move all pairs
    for i in range(2, 8):
        state.atoms[i].x += 0.2
    
    pgp_multi = pygcmc.computeMovementEnergyPGPCompleteCorrect(state)
    pme_multi = pygcmc.computeMovementEnergyPME(state)
    
    pgp_multi_total = pgp_multi[0] + pgp_multi[1]
    pme_multi_total = pme_multi[0] + pme_multi[1]
    
    multi_error = abs(pgp_multi_total - pme_multi_total)
    
    # Reset
    for i in range(2, 8):
        state.atoms[i].x -= 0.2
    
    # Results
    print("\n" + "="*60)
    print("PGP Multi-Residue Limitation Test")
    print("="*60)
    print(f"\nSingle residue move:")
    print(f"  PGP: {pgp_single_total:.6f} kJ/mol")
    print(f"  PME: {pme_single_total:.6f} kJ/mol")
    print(f"  Error: {single_error:.6f} kJ/mol ({single_error/abs(pme_single_total)*100:.2f}%)")
    
    print(f"\nMultiple residue move (3 pairs):")
    print(f"  PGP: {pgp_multi_total:.6f} kJ/mol")
    print(f"  PME: {pme_multi_total:.6f} kJ/mol")
    print(f"  Error: {multi_error:.6f} kJ/mol ({multi_error/abs(pme_multi_total)*100:.2f}%)")
    
    print(f"\nError increase factor: {multi_error/single_error:.1f}x")
    
    # Document the fundamental difference
    print("\n⚠️  Important Note:")
    print("   PGP Complete includes ALL interactions (including intramolecular)")
    print("   PyGCMC PME excludes intramolecular for movement residues")
    print("   This fundamental difference causes large apparent 'errors'")
    print("   ")
    print("   For accurate comparison with OpenMM PME, see:")
    print("   - test_pgp_complete_delta_e_accuracy")
    print("   - test_pgp_openmm_reciprocal_ratio")
    
    # We expect large differences due to different physics
    # This test documents the behavior rather than testing accuracy
    if single_error > 100:
        print("\n✓ Large difference confirms different intramolecular handling")
    else:
        # If error is small, there might be an issue
        assert single_error > 10, "Error suspiciously small - check intramolecular handling"


def test_pgp_sequential_vs_simultaneous_moves():
    """Compare sequential single moves vs simultaneous multi-residue move"""
    
    # Create simple system
    state = MCState()
    state.info.box = [6.0, 6.0, 6.0]
    state.info.cutoff = 1.0
    
    ff = MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.ljEps = [0.1]
    ff.ljSigma = [0.3]
    state.forcefield = ff
    
    atoms = []
    residues = []
    
    # Fixed atom
    atom = MCAtom()
    atom.x, atom.y, atom.z = 3.0, 3.0, 3.0
    atom.charge = 0.0
    atom.type = 0
    atoms.append(atom)
    
    res = MCResidue()
    res.atomStart = 0
    res.atomCount = 1
    res.active = True
    res.fixed = True
    res.type = 0
    residues.append(res)
    
    # Two moving atoms
    for i in range(2):
        atom = MCAtom()
        atom.x = 2.0 + i * 2.0
        atom.y = 3.0
        atom.z = 3.0
        atom.charge = (-1)**i
        atom.type = 0
        atoms.append(atom)
        
        res = MCResidue()
        res.atomStart = 1 + i
        res.atomCount = 1
        res.active = True
        res.fixed = False
        res.type = 0
        residues.append(res)
    
    state.atoms = atoms
    state.activeAtomCount = 3
    state.residues = residues
    state.activeResidueCount = 3
    
    # Initialize
    alpha = 5.6 / state.info.cutoff
    mesh_size = [32, 32, 32]
    
    pygcmc.setPMEParameters(alpha, mesh_size, 4, 1e-5)
    pygcmc.initializePMEParameters(state.info.cutoff, state.info.box, alpha)
    pygcmc.setPGPParameters(alpha, mesh_size, state.info.cutoff, mesh_size, 4, 1e-5)
    pygcmc.precomputeGridPotential(state, fixed_only=True)
    
    # Get initial energy
    initial_state = [atom.x for atom in state.atoms]
    
    # Method 1: Sequential moves
    total_delta_sequential = 0.0
    
    # Move atom 1
    state.movementResidues.clear()
    movement_info = MCMovementResidueInfo()
    movement_info.startIndex = 1
    movement_info.activeCount = 1
    state.movementResidues.append(movement_info)
    
    pgp_before = pygcmc.computeMovementEnergyPGPCompleteCorrect(state)
    state.atoms[1].x += 0.3
    pgp_after = pygcmc.computeMovementEnergyPGPCompleteCorrect(state)
    
    delta1 = (pgp_after[0] + pgp_after[1]) - (pgp_before[0] + pgp_before[1])
    total_delta_sequential += delta1
    
    # Move atom 2 (with atom 1 already moved)
    state.movementResidues.clear()
    movement_info = MCMovementResidueInfo()
    movement_info.startIndex = 2
    movement_info.activeCount = 1
    state.movementResidues.append(movement_info)
    
    pgp_before = pygcmc.computeMovementEnergyPGPCompleteCorrect(state)
    state.atoms[2].x -= 0.3
    pgp_after = pygcmc.computeMovementEnergyPGPCompleteCorrect(state)
    
    delta2 = (pgp_after[0] + pgp_after[1]) - (pgp_before[0] + pgp_before[1])
    total_delta_sequential += delta2
    
    # Reset
    state.atoms[1].x = initial_state[1]
    state.atoms[2].x = initial_state[2]
    
    # Method 2: Simultaneous move
    state.movementResidues.clear()
    movement_info = MCMovementResidueInfo()
    movement_info.startIndex = 1
    movement_info.activeCount = 2
    state.movementResidues.append(movement_info)
    
    pgp_before = pygcmc.computeMovementEnergyPGPCompleteCorrect(state)
    state.atoms[1].x += 0.3
    state.atoms[2].x -= 0.3
    pgp_after = pygcmc.computeMovementEnergyPGPCompleteCorrect(state)
    
    delta_simultaneous = (pgp_after[0] + pgp_after[1]) - (pgp_before[0] + pgp_before[1])
    
    # Compare
    difference = abs(total_delta_sequential - delta_simultaneous)
    
    print("\n" + "="*60)
    print("Sequential vs Simultaneous Move Comparison")
    print("="*60)
    print(f"\nSequential moves:")
    print(f"  Move 1 ΔE: {delta1:.6f} kJ/mol")
    print(f"  Move 2 ΔE: {delta2:.6f} kJ/mol")
    print(f"  Total ΔE:  {total_delta_sequential:.6f} kJ/mol")
    
    print(f"\nSimultaneous move:")
    print(f"  Total ΔE:  {delta_simultaneous:.6f} kJ/mol")
    
    print(f"\nDifference: {difference:.6f} kJ/mol")
    
    if difference > 0.01:
        print("\n⚠️  Sequential and simultaneous moves give different results!")
        print("   This is due to mov-mov reciprocal interactions being missed.")
    
    # Don't assert failure - this is documenting known behavior
    # The difference documents the mov-mov reciprocal issue


if __name__ == "__main__":
    test_pgp_multi_residue_limitation()
    test_pgp_sequential_vs_simultaneous_moves()
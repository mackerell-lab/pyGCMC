"""
Diagnose PGP LJ calculation issues

Since PME Complete matches OpenMM for LJ, but PGP shows large errors,
we need to understand what PGP is doing differently.
"""

import pytest
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pygcmc
from pygcmc import MCState, MCAtom, MCResidue, MCForceField, MCMovementResidueInfo
from pygcmc import initializePMEParameters, computeSystemEnergyPMEComplete
from pygcmc import setPMEParameters, setPGPParameters
from pygcmc import precomputeGridPotential, calculateMoleculeEnergy
from pygcmc import computeMovementEnergyPME
from pygcmc import computeSystemEnergyCutoff, getTotalEnergyComponents


def create_lj_only_system():
    """Create a system with only LJ interactions (no charges)"""
    state = MCState()
    state.info.box = [5.0, 5.0, 5.0]
    state.info.cutoff = 2.0
    
    # Force field with LJ but no charges
    ff = MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.ljEps = [1.0]  # epsilon in kJ/mol
    ff.ljSigma = [0.35]  # sigma in nm
    state.forcefield = ff
    
    atoms = []
    residues = []
    
    # Fixed particles
    fixed_positions = [
        [2.0, 2.5, 2.5],
        [3.0, 2.5, 2.5],
    ]
    
    for i, pos in enumerate(fixed_positions):
        atom = MCAtom()
        atom.x, atom.y, atom.z = pos
        atom.charge = 0.0  # No charge
        atom.type = 0
        atoms.append(atom)
        
        res = MCResidue()
        res.active = True
        res.fixed = True
        res.atomStart = i
        res.atomCount = 1
        res.type = 0
        residues.append(res)
    
    # Moveable particle
    moveable_atom = MCAtom()
    moveable_atom.x, moveable_atom.y, moveable_atom.z = 2.5, 2.5, 2.5
    moveable_atom.charge = 0.0  # No charge
    moveable_atom.type = 0
    atoms.append(moveable_atom)
    
    moveable_res = MCResidue()
    moveable_res.active = True
    moveable_res.fixed = False
    moveable_res.atomStart = 2
    moveable_res.atomCount = 1
    moveable_res.type = 0
    residues.append(moveable_res)
    
    state.atoms = atoms
    state.activeAtomCount = 3
    state.residues = residues
    state.activeResidueCount = 3
    
    # Set up movement residues
    movement_info = MCMovementResidueInfo()
    movement_info.startIndex = 2
    movement_info.activeCount = 1
    state.movementResidues = [movement_info]
    
    return state


def test_pgp_lj_only_system():
    """Test PGP with pure LJ system"""
    
    print("\n" + "="*70)
    print("PGP LJ-Only System Test")
    print("="*70)
    
    state = create_lj_only_system()
    
    # Calculate direct energy for reference
    computeSystemEnergyCutoff(state)
    elec_direct, vdw_direct = getTotalEnergyComponents(state)
    print(f"\nDirect calculation (cutoff):")
    print(f"  Electrostatic: {elec_direct:.8f} kJ/mol (should be 0)")
    print(f"  VdW:          {vdw_direct:.8f} kJ/mol")
    
    # Set up PME
    alpha = 2.2
    mesh_size = [64, 64, 64]
    spline_order = 4
    
    setPMEParameters(alpha, mesh_size, spline_order)
    initializePMEParameters(state.info.cutoff, state.info.box, alpha, mesh_size, spline_order)
    
    # Calculate PME Complete
    elec_pme, vdw_pme, total_pme = computeSystemEnergyPMEComplete(state)
    print(f"\nPME Complete:")
    print(f"  Electrostatic: {elec_pme:.8f} kJ/mol (should be 0)")
    print(f"  VdW:          {vdw_pme:.8f} kJ/mol")
    print(f"  Total:        {total_pme:.8f} kJ/mol")
    
    # Calculate PME movement energy
    pme_move = computeMovementEnergyPME(state)
    pme_move_elec = pme_move[0]
    pme_move_vdw = pme_move[1]
    print(f"\nPME Movement Energy:")
    print(f"  Electrostatic: {pme_move_elec:.8f} kJ/mol (should be 0)")
    print(f"  VdW:          {pme_move_vdw:.8f} kJ/mol")
    
    # Set up PGP
    setPGPParameters(alpha, mesh_size, state.info.cutoff, mesh_size, spline_order, 1e-6)
    
    print("\nPrecomputing PGP grid...")
    precomputeGridPotential(state, fixed_only=True)
    
    # Calculate PGP energy
    pgp_energy = calculateMoleculeEnergy(state)
    print(f"\nPGP energy: {pgp_energy:.8f} kJ/mol")
    
    # Compare
    print(f"\nComparison:")
    print(f"  Direct VdW:      {vdw_direct:.8f} kJ/mol")
    print(f"  PME Complete VdW: {vdw_pme:.8f} kJ/mol")
    print(f"  PME Movement VdW: {pme_move_vdw:.8f} kJ/mol")
    print(f"  PGP total:       {pgp_energy:.8f} kJ/mol")
    
    # PGP should match PME movement energy for LJ
    diff = abs(pgp_energy - pme_move_vdw)
    print(f"\nPGP vs PME Movement VdW difference: {diff:.8f} kJ/mol")
    
    # Move particle and check energy change
    print("\n" + "-"*50)
    print("Testing energy change on movement...")
    
    initial_pgp = pgp_energy
    initial_pme_move_vdw = pme_move_vdw
    
    # Move particle
    state.atoms[2].x += 0.1
    
    # Recalculate energies
    pgp_moved = calculateMoleculeEnergy(state)
    pme_move_moved = computeMovementEnergyPME(state)[1]
    
    pgp_delta = pgp_moved - initial_pgp
    pme_delta = pme_move_moved - initial_pme_move_vdw
    
    print(f"\nAfter moving particle by 0.1 nm:")
    print(f"  PGP delta:         {pgp_delta:.8f} kJ/mol")
    print(f"  PME Movement delta: {pme_delta:.8f} kJ/mol")
    print(f"  Difference:        {abs(pgp_delta - pme_delta):.8f} kJ/mol")
    
    rel_error = abs(pgp_delta - pme_delta) / abs(pme_delta) * 100 if pme_delta != 0 else 0
    print(f"  Relative error:    {rel_error:.4f}%")


def test_pgp_mixed_system():
    """Test PGP with both charges and LJ"""
    
    print("\n" + "="*70)
    print("PGP Mixed System Test (Charges + LJ)")
    print("="*70)
    
    state = MCState()
    state.info.box = [5.0, 5.0, 5.0]
    state.info.cutoff = 2.0
    
    # Force field with both
    ff = MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.ljEps = [1.0]
    ff.ljSigma = [0.35]
    state.forcefield = ff
    
    atoms = []
    residues = []
    
    # Fixed particles with charges
    fixed_data = [
        ([2.0, 2.5, 2.5], 1.0),
        ([3.0, 2.5, 2.5], -1.0),
    ]
    
    for i, (pos, charge) in enumerate(fixed_data):
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
    
    # Moveable particle with charge
    moveable_atom = MCAtom()
    moveable_atom.x, moveable_atom.y, moveable_atom.z = 2.5, 2.5, 2.5
    moveable_atom.charge = 0.5
    moveable_atom.type = 0
    atoms.append(moveable_atom)
    
    moveable_res = MCResidue()
    moveable_res.active = True
    moveable_res.fixed = False
    moveable_res.atomStart = 2
    moveable_res.atomCount = 1
    moveable_res.type = 0
    residues.append(moveable_res)
    
    state.atoms = atoms
    state.activeAtomCount = 3
    state.residues = residues
    state.activeResidueCount = 3
    
    movement_info = MCMovementResidueInfo()
    movement_info.startIndex = 2
    movement_info.activeCount = 1
    state.movementResidues = [movement_info]
    
    # Set up PME
    alpha = 2.2
    mesh_size = [64, 64, 64]
    spline_order = 4
    
    setPMEParameters(alpha, mesh_size, spline_order)
    initializePMEParameters(state.info.cutoff, state.info.box, alpha, mesh_size, spline_order)
    
    # Calculate PME movement energy
    pme_move = computeMovementEnergyPME(state)
    pme_move_elec = pme_move[0]
    pme_move_vdw = pme_move[1]
    pme_move_total = pme_move_elec + pme_move_vdw
    
    print(f"\nPME Movement Energy:")
    print(f"  Electrostatic: {pme_move_elec:.8f} kJ/mol")
    print(f"  VdW:          {pme_move_vdw:.8f} kJ/mol")
    print(f"  Total:        {pme_move_total:.8f} kJ/mol")
    
    # Set up PGP
    setPGPParameters(alpha, mesh_size, state.info.cutoff, mesh_size, spline_order, 1e-6)
    precomputeGridPotential(state, fixed_only=True)
    
    # Calculate PGP energy
    pgp_energy = calculateMoleculeEnergy(state)
    print(f"\nPGP energy: {pgp_energy:.8f} kJ/mol")
    
    # Break down the difference
    print(f"\nBreakdown:")
    print(f"  PGP total:           {pgp_energy:.8f} kJ/mol")
    print(f"  PME elec + vdw:      {pme_move_total:.8f} kJ/mol")
    print(f"  Difference:          {pgp_energy - pme_move_total:.8f} kJ/mol")
    
    # Test if PGP includes LJ at all
    print(f"\nDoes PGP include LJ?")
    print(f"  PGP - PME elec:      {pgp_energy - pme_move_elec:.8f} kJ/mol")
    print(f"  Expected (PME vdw):  {pme_move_vdw:.8f} kJ/mol")
    
    # Move and test again
    initial_pgp = pgp_energy
    initial_pme_total = pme_move_total
    
    state.atoms[2].x += 0.1
    
    pgp_moved = calculateMoleculeEnergy(state)
    pme_moved = computeMovementEnergyPME(state)
    pme_moved_total = pme_moved[0] + pme_moved[1]
    
    pgp_delta = pgp_moved - initial_pgp
    pme_delta = pme_moved_total - initial_pme_total
    
    print(f"\nAfter movement:")
    print(f"  PGP delta:      {pgp_delta:.8f} kJ/mol")
    print(f"  PME delta:      {pme_delta:.8f} kJ/mol")
    print(f"  Difference:     {abs(pgp_delta - pme_delta):.8f} kJ/mol")
    
    rel_error = abs(pgp_delta - pme_delta) / abs(pme_delta) * 100 if pme_delta != 0 else 0
    print(f"  Relative error: {rel_error:.4f}%")


def test_pgp_lj_distance_scan():
    """Scan PGP LJ energy at different distances"""
    
    print("\n" + "="*70)
    print("PGP LJ Distance Scan")
    print("="*70)
    
    state = MCState()
    state.info.box = [10.0, 10.0, 10.0]  # Large box
    state.info.cutoff = 4.0
    
    # Force field
    ff = MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.ljEps = [1.0]
    ff.ljSigma = [0.35]
    state.forcefield = ff
    
    atoms = []
    residues = []
    
    # Fixed particle at origin
    atom1 = MCAtom()
    atom1.x, atom1.y, atom1.z = 5.0, 5.0, 5.0
    atom1.charge = 0.0
    atom1.type = 0
    atoms.append(atom1)
    
    res1 = MCResidue()
    res1.active = True
    res1.fixed = True
    res1.atomStart = 0
    res1.atomCount = 1
    res1.type = 0
    residues.append(res1)
    
    # Moveable particle
    atom2 = MCAtom()
    atom2.x, atom2.y, atom2.z = 5.35, 5.0, 5.0  # Start at sigma
    atom2.charge = 0.0
    atom2.type = 0
    atoms.append(atom2)
    
    res2 = MCResidue()
    res2.active = True
    res2.fixed = False
    res2.atomStart = 1
    res2.atomCount = 1
    res2.type = 0
    residues.append(res2)
    
    state.atoms = atoms
    state.activeAtomCount = 2
    state.residues = residues
    state.activeResidueCount = 2
    
    movement_info = MCMovementResidueInfo()
    movement_info.startIndex = 1
    movement_info.activeCount = 1
    state.movementResidues = [movement_info]
    
    # Set up PME
    alpha = 2.2
    mesh_size = [64, 64, 64]
    spline_order = 4
    
    setPMEParameters(alpha, mesh_size, spline_order)
    initializePMEParameters(state.info.cutoff, state.info.box, alpha, mesh_size, spline_order)
    setPGPParameters(alpha, mesh_size, state.info.cutoff, mesh_size, spline_order, 1e-6)
    
    # Distances to test (in nm)
    distances = [0.35, 0.4, 0.5, 0.6, 0.8, 1.0, 1.5, 2.0, 2.5, 3.0]
    
    print(f"\n{'Distance':>10} | {'Direct LJ':>12} | {'PME Move LJ':>12} | {'PGP Total':>12} | {'PGP-PME':>12}")
    print("-" * 65)
    
    for d in distances:
        # Set position
        state.atoms[1].x = 5.0 + d
        
        # Calculate direct LJ
        computeSystemEnergyCutoff(state)
        _, vdw_direct = getTotalEnergyComponents(state)
        
        # Calculate PME movement LJ
        pme_move = computeMovementEnergyPME(state)
        pme_move_vdw = pme_move[1]
        
        # Precompute PGP grid and calculate
        precomputeGridPotential(state, fixed_only=True)
        pgp_energy = calculateMoleculeEnergy(state)
        
        diff = pgp_energy - pme_move_vdw
        
        print(f"{d:10.2f} | {vdw_direct:12.6f} | {pme_move_vdw:12.6f} | {pgp_energy:12.6f} | {diff:12.6f}")


if __name__ == "__main__":
    test_pgp_lj_only_system()
    test_pgp_mixed_system()
    test_pgp_lj_distance_scan()
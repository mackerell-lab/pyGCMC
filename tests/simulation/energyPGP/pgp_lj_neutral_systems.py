# tests/simulation/energyPGP/pgp_lj_neutral_systems.py
"""
Test PGP with neutral particles (no charges) to verify LJ calculations.

This test ensures that PGP correctly calculates Lennard-Jones interactions
when there are no electrostatic contributions in the system.
"""

import pytest
import math
import pygcmc
from pygcmc import MCAtom, MCResidue, MCState
from pygcmc import setPGPParameters, initializePMEParameters, precomputeGridPotential
from pygcmc import computeSystemEnergyPGP, computeMovementEnergyPGP
from pygcmc import computeSystemVdwEnergyCutoff
import numpy as np


def create_lj_only_system(box_size=5.0, cutoff=1.2):
    """Create a simple system with neutral atoms that have LJ parameters."""
    state = MCState()
    
    # System info
    state.info.box = [box_size, box_size, box_size]
    state.info.cutoff = cutoff
    state.info.use_switching = False  # Disable switching for simplicity
    
    # Force field with LJ parameters
    state.forcefield.numTotalTypes = 1
    state.forcefield.numMovementTypes = 1
    
    # LJ parameters for argon: epsilon = 0.996 kJ/mol, sigma = 0.3405 nm
    eps_ar = 0.996  # kJ/mol
    sigma_ar = 0.3405  # nm
    state.forcefield.ljEps = [eps_ar]  # For 1x1 matrix
    state.forcefield.ljSigma = [sigma_ar]
    
    # Create atoms with NO charges
    atoms = []
    
    # Fixed argon atoms in a small cluster
    fixed_positions = [
        [2.0, 2.5, 2.5],
        [2.8, 2.5, 2.5],  # 0.8 nm apart
        [2.4, 3.2, 2.5],  # Forms a triangle
    ]
    
    # Add fixed atoms
    for i, pos in enumerate(fixed_positions):
        atom = MCAtom()
        atom.x, atom.y, atom.z = pos
        atom.charge = 0.0  # NO CHARGE
        atom.type = 0  # AR type
        atoms.append(atom)
    
    # Add one moving atom
    moving_atom = MCAtom()
    moving_atom.x = 3.5
    moving_atom.y = 2.5
    moving_atom.z = 2.5
    moving_atom.charge = 0.0  # NO CHARGE
    moving_atom.type = 0  # AR type
    atoms.append(moving_atom)
    
    state.atoms = atoms
    state.activeAtomCount = len(atoms)
    
    # Create residues
    residues = []
    
    # Fixed residues (one atom each)
    for i in range(3):
        res = MCResidue()
        res.active = True
        res.fixed = True  # Fixed for PGP precomputation
        res.atomStart = i
        res.atomCount = 1
        res.type = 0
        residues.append(res)
    
    # Moving residue
    moving_res = MCResidue()
    moving_res.active = True
    moving_res.fixed = False  # This will move
    moving_res.atomStart = 3
    moving_res.atomCount = 1
    moving_res.type = 0
    residues.append(moving_res)
    
    state.residues = residues
    state.activeResidueCount = len(residues)
    
    return state


def calculate_lj_energy_manual(atoms, epsilon, sigma, cutoff):
    """Manually calculate LJ energy for verification."""
    total_energy = 0.0
    
    for i in range(len(atoms)):
        for j in range(i + 1, len(atoms)):
            dx = atoms[i].x - atoms[j].x
            dy = atoms[i].y - atoms[j].y
            dz = atoms[i].z - atoms[j].z
            r2 = dx*dx + dy*dy + dz*dz
            r = math.sqrt(r2)
            
            if r < cutoff:
                # Standard LJ calculation
                sigma_over_r = sigma / r
                sigma_over_r_6 = sigma_over_r ** 6
                sigma_over_r_12 = sigma_over_r_6 ** 2
                
                lj_energy = 4.0 * epsilon * (sigma_over_r_12 - sigma_over_r_6)
                total_energy += lj_energy
    
    return total_energy


def test_pgp_lj_only_system():
    """Test PGP with neutral atoms to verify LJ calculation."""
    
    print("\n=== Test: PGP with LJ-only (neutral) system ===")
    
    # Create system
    system = create_lj_only_system()
    eps_ar = system.forcefield.ljEps[0]
    sigma_ar = system.forcefield.ljSigma[0]
    
    print(f"Created system with {system.activeAtomCount} argon atoms (no charges)")
    print(f"LJ parameters: epsilon = {eps_ar} kJ/mol, sigma = {sigma_ar} nm")
    print(f"Box: {system.info.box}, Cutoff: {system.info.cutoff} nm")
    
    # Print atom positions
    print("\nAtom positions:")
    for i, atom in enumerate(system.atoms):
        residue = system.residues[i]
        status = "Fixed" if residue.fixed else "Moving"
        print(f"  Atom {i} ({status}): ({atom.x:.3f}, {atom.y:.3f}, {atom.z:.3f}), charge = {atom.charge}")
    
    # Calculate manual LJ energy for verification
    manual_lj = calculate_lj_energy_manual(
        system.atoms, eps_ar, sigma_ar, 
        system.info.cutoff
    )
    print(f"\nManual LJ calculation: {manual_lj:.6f} kJ/mol")
    
    # Method 1: Direct LJ calculation
    print("\n1. Direct LJ calculation:")
    computeSystemVdwEnergyCutoff(system)
    direct_lj_total = 0.0
    for i, res in enumerate(system.residues):
        print(f"   Residue {i}: LJ energy = {res.energy_vdw:.6f} kJ/mol")
        direct_lj_total += res.energy_vdw
    print(f"   Total direct LJ: {direct_lj_total:.6f} kJ/mol")
    
    # Method 2: PGP calculation
    print("\n2. PGP calculation:")
    
    # Initialize PGP parameters
    alpha = 0.3  # Doesn't matter much for neutral system
    mesh_size = [32, 32, 32]
    spline_order = 4
    tolerance = 1e-5
    
    setPGPParameters(
        alpha=alpha, 
        meshSize=mesh_size, 
        potential_cutoff=system.info.cutoff,
        potentialGridSize=mesh_size, 
        splineOrder=spline_order, 
        tolerance=tolerance
    )
    
    # Initialize PME parameters (required by PGP)
    initializePMEParameters(system.info.cutoff, system.info.box, alpha)
    
    # Precompute grid potential (should be zero for neutral atoms)
    print("   Precomputing grid potential for fixed atoms...")
    precomputeGridPotential(system, fixed_only=True)
    
    # Compute system energy using PGP
    computeSystemEnergyPGP(system)
    
    print(f"   PGP reciprocal energy: {system.ewald_energy.get('reciprocal', 0.0):.6f} kJ/mol (should be ~0)")
    print(f"   PGP real space energy: {system.ewald_energy.get('real_space', 0.0):.6f} kJ/mol (should be ~0)")
    print(f"   PGP self energy: {system.ewald_energy.get('self', 0.0):.6f} kJ/mol (should be 0)")
    
    # Get LJ from residues
    pgp_lj_total = 0.0
    for i, res in enumerate(system.residues):
        pgp_lj_total += res.energy_vdw
    print(f"   PGP LJ total: {pgp_lj_total:.6f} kJ/mol")
    print(f"   PGP total energy: {system.ewald_energy.get('total', 0.0):.6f} kJ/mol")
    
    # Verify electrostatic components are zero
    assert abs(system.ewald_energy.get('reciprocal', 0.0)) < 1e-10, \
        f"Reciprocal energy should be zero for neutral system, got {system.ewald_energy.get('reciprocal', 0.0)}"
    assert abs(system.ewald_energy.get('real_space', 0.0)) < 1e-10, \
        f"Real space energy should be zero for neutral system, got {system.ewald_energy.get('real_space', 0.0)}"
    assert abs(system.ewald_energy.get('self', 0.0)) < 1e-10, \
        f"Self energy should be zero for neutral system, got {system.ewald_energy.get('self', 0.0)}"
    
    # Verify LJ energy matches
    assert abs(pgp_lj_total - direct_lj_total) < 1e-6, \
        f"PGP LJ ({pgp_lj_total:.6f}) doesn't match direct LJ ({direct_lj_total:.6f})"
    
    # Test movement energy calculation
    print("\n3. Testing movement energy:")
    print("   Note: Movement energy functions return energy CHANGES from initial state")
    print("   The initial call returns 0 by design (reference state)")
    
    # Calculate energy before setting movement
    computeSystemEnergyPGP(system)
    energy_before = system.ewald_energy.get('total', 0.0)
    print(f"\n   Total system energy before move: {energy_before:.6f} kJ/mol")
    
    # Move the atom
    moving_atom = system.atoms[3]
    print(f"\n   Moving atom from ({moving_atom.x:.3f}, {moving_atom.y:.3f}, {moving_atom.z:.3f})")
    moving_atom.x += 0.2  # Move 0.2 nm away
    print(f"   to ({moving_atom.x:.3f}, {moving_atom.y:.3f}, {moving_atom.z:.3f})")
    
    # Calculate energy after movement
    computeSystemEnergyPGP(system)
    energy_after = system.ewald_energy.get('total', 0.0)
    print(f"\n   Total system energy after move: {energy_after:.6f} kJ/mol")
    
    # Calculate delta
    delta_total = energy_after - energy_before
    print(f"\n   Delta total energy: {delta_total:.6f} kJ/mol")
    
    # For LJ-only system, all energy change is from LJ
    print(f"   Since system is neutral, all energy change is from LJ interactions")
    
    # Verify delta is reasonable (should be positive as we moved away)
    assert delta_total > 0, f"LJ energy should increase when moving away, got delta = {delta_total}"
    
    print("\n✅ PGP correctly handles LJ-only systems!")



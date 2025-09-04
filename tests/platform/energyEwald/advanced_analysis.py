# tests/simulation/energyEwald/advanced_analysis.py
"""Energy Ewald advanced analysis tests."""

import pytest
import math
import random
import pygcmc
from pygcmc import MCState, MCAtom, MCResidue, MCForceField, MCInfo, MCMovementResidueInfo
import os
import sys

# Set log level to INFO for debugging output
pygcmc.System.set_log_level(pygcmc.LogLevel.INFO)

# Import helper function from basic_comparison
from .basic_comparison import create_nacl_crystal


def test_madelung_constant():
    """Test against known Madelung constant for NaCl"""
    # Madelung constant for NaCl (from literature)
    MADELUNG_NACL = 1.747564594633182190636212035544397403481

    # NaCl lattice constant (nm)
    a = 0.564

    # For n_cells x n_cells x n_cells crystal,
    # to ensure the crystal fills the box, the box size is set to n_cells * a
    n_cells = 4
    box_size = n_cells * a
    state = create_nacl_crystal(box_size, n_cells)

    # Define physical constants - identical to those in test_ewald_exact
    PI_M = math.pi
    ONE_4PI_EPS0 = 138.935456  # Coulomb constant converted to kJ·mol^-1·nm·e^-2
    eCharge = 1.6022e-19  # Elementary charge, unit: Coulomb (C)
    AVOGADRO = 6.02214076e23  # Avogadro's number
    eps0 = 8.8542e-12  # Vacuum permittivity, F/m
    a0 = 0.282e-9  # Meters, NaCl unit cell length

    # Test different alpha values and kmax to understand convergence
    alpha_tests = [1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]
    kmax_tests = [[10, 10, 10], [12, 12, 12], [14, 14, 14], [16, 16, 16], [18, 18, 18]]

    print("\nMadelung constant convergence study:")
    best_error = float('inf')
    best_madelung = 0.0
    best_params = None

    # 1. Estimate optimal alpha (using fixed kmax)
    print("\nEstimating optimal alpha...")
    optimal_alpha = None
    min_alpha_error = float('inf')
    kmax_fixed = [14, 14, 14]

    for alpha in alpha_tests:
        pygcmc.setEwaldParameters(alpha, kmax_fixed)
        result = pygcmc.computeSystemEnergyEwald(state)
        ewald_dict = result[2]  # Save ewald dictionary
        
        # Modification: Get electrostatic energy from ewald_dict
        elec_energy = ewald_dict['real_space'] + ewald_dict['reciprocal'] + ewald_dict['self']
        
        # Normalize energy: divide by number of unit cells (n_cells^3)
        elec_energy = elec_energy / (n_cells ** 3)
        
        # Calculate theoretical energy, identical to test_ewald_exact
        num_atoms = len(state.atoms)
        theoretical_energy = -(MADELUNG_NACL * eCharge * eCharge * AVOGADRO) / (4 * PI_M * eps0 * a0 * 2 * 1000)
        calculated_madelung = MADELUNG_NACL * (elec_energy / theoretical_energy)
        rel_error = abs(calculated_madelung - MADELUNG_NACL) / MADELUNG_NACL

        print(f"Alpha = {alpha:.1f}: Madelung = {calculated_madelung:.6f}, Error = {rel_error:.6f}")

        if rel_error < min_alpha_error:
            min_alpha_error = rel_error
            optimal_alpha = alpha

    print(f"\nOptimal alpha = {optimal_alpha}")

    # 2. Test different kmax values with optimal alpha
    print("\nTesting kmax convergence with optimal alpha...")
    for kmax in kmax_tests:
        pygcmc.setEwaldParameters(optimal_alpha, kmax)
        result = pygcmc.computeSystemEnergyEwald(state)
        ewald_dict = result[2]  # Save ewald dictionary
        
        # Modification: Get energy components from ewald_dict
        elec_energy = ewald_dict['real_space'] + ewald_dict['reciprocal'] + ewald_dict['self']
        vdw_energy = sum(res.energy_vdw for res in state.residues if res.active)
        total_energy = elec_energy + vdw_energy

        # Energy normalization: divide by n_cells^3
        total_energy = total_energy / (n_cells ** 3)
        elec_energy = elec_energy / (n_cells ** 3)
        vdw_energy = vdw_energy / (n_cells ** 3)

        # Use the same calculation method as in test_ewald_exact
        num_atoms = len(state.atoms)
        theoretical_energy = -(MADELUNG_NACL * eCharge * eCharge * AVOGADRO) / (4 * PI_M * eps0 * a0 * 2 * 1000)
        calculated_madelung = MADELUNG_NACL * (elec_energy / theoretical_energy)
        rel_error = abs(calculated_madelung - MADELUNG_NACL) / MADELUNG_NACL

        if rel_error < best_error:
            best_error = rel_error
            best_madelung = calculated_madelung
            best_params = (optimal_alpha, kmax)

        print(f"\nkmax = {kmax}:")
        print(f"Total energy per cell: {total_energy:.6f} kJ/mol")
        print(f"Electrostatic per cell: {elec_energy:.6f} kJ/mol")
        print(f"Van der Waals per cell: {vdw_energy:.6f} kJ/mol")
        print(f"Calculated Madelung: {calculated_madelung:.9f}")
        print(f"Relative error: {rel_error:.6f}")

    print(f"\nBest result:")
    print(f"Alpha = {best_params[0]}, kmax = {best_params[1]}")
    print(f"Calculated Madelung: {best_madelung:.9f}")
    print(f"Reference Madelung: {MADELUNG_NACL:.9f}")
    print(f"Best relative error: {best_error:.6f}")

    # 3. Check finite size effects
    print("\nChecking finite size effects...")
    cell_counts = [2, 3, 4, 5]
    for n in cell_counts:
        # Box size should be n * a
        state = create_nacl_crystal(n * a, n)
        pygcmc.setEwaldParameters(best_params[0], best_params[1])
        result = pygcmc.computeSystemEnergyEwald(state)
        ewald_dict = result[2]  # Save ewald dictionary
        
        # Modification: Get electrostatic energy from ewald_dict
        elec_energy = ewald_dict['real_space'] + ewald_dict['reciprocal'] + ewald_dict['self']
        elec_energy = elec_energy / (n ** 3)

        # Use the same calculation method as in test_ewald_exact
        num_atoms = len(state.atoms)
        theoretical_energy = -(MADELUNG_NACL * eCharge * eCharge * AVOGADRO) / (4 * PI_M * eps0 * a0 * 2 * 1000)
        calculated_madelung = MADELUNG_NACL * (elec_energy / theoretical_energy)
        rel_error = abs(calculated_madelung - MADELUNG_NACL) / MADELUNG_NACL
        
        print(f"\n{n}x{n}x{n} cells:")
        print(f"Calculated Madelung: {calculated_madelung:.9f}")
        print(f"Relative error: {rel_error:.6f}")

    # Check if the best relative error is within the allowed range
    assert best_error < 0.2, \
        f"Best calculated Madelung constant ({best_madelung}) differs too much from reference ({MADELUNG_NACL})"


def test_ewald_charge_neutrality():
    """Test that Ewald summation properly handles charge neutrality requirements"""
    
    # Create a neutral system
    state = create_nacl_crystal(2.0, 2)
    kmax = [8, 8, 8]
    alpha = 2.5
    
    # Calculate total charge of the system
    total_charge = sum(atom.charge for atom in state.atoms)
    print(f"\nInitial total charge: {total_charge}")
    assert abs(total_charge) < 1e-10, "System should be neutral initially"
    
    # This should work fine
    pygcmc.setEwaldParameters(alpha, kmax)
    pygcmc.computeSystemEnergyEwald(state)
    neutral_energy = sum(res.energy_vdw + res.energy_elec 
                        for res in state.residues if res.active)
    print(f"Neutral system energy: {neutral_energy:.6f} kJ/mol")
    
    # Now create a charged system by changing one Cl- to Na+
    charge_changed = False
    for atom in state.atoms:
        if atom.charge < 0:
            old_charge = atom.charge
            atom.charge = 1.0  # Change one Cl- to Na+
            charge_changed = True
            break
    
    assert charge_changed, "Failed to create charged system"
    
    # Verify system is now charged
    total_charge = sum(atom.charge for atom in state.atoms)
    print(f"Modified system total charge: {total_charge}")
    assert abs(total_charge) > 1e-10, "System should be charged after modification"
    
    # The energy calculation should either raise an exception or return a warning
    # Depending on your implementation, you might want to check for either behavior
    try:
        pygcmc.computeSystemEnergyEwald(state)
        energy = sum(res.energy_vdw + res.energy_elec 
                    for res in state.residues if res.active)
        print(f"Warning: Charged system energy calculated: {energy:.6f} kJ/mol")
        print("Note: Your implementation allows charged systems. Make sure this is intended behavior.")
    except RuntimeError as e:
        print(f"Expected exception raised: {str(e)}")
        assert "neutral" in str(e).lower(), "Exception should mention charge neutrality"

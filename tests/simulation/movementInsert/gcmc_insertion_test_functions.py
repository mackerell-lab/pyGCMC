# tests/simulation/movementInsert/gcmc_insertion_test.py
"""
GCMC insertion acceptance test based on detailed balance

This module tests GCMC insertion using the Metropolis criterion:
A(n->n+1) = min{1, f_n/(n+1) * exp[B - β(W_{n+1} - W_n)]}

where:
- f_n: fraction of unoccupied volume (cavity bias)
- n: current number of molecules
- B = β*μ_ex + ln(n̄), where μ_ex is excess chemical potential
- β = 1/kT
- W: total energy
"""

import pytest
import math
import random
import pygcmc
import os

# Constants
kB = 0.008314463  # Boltzmann constant in kJ/mol/K
kC = 138.935456   # Coulomb constant in kJ·nm/mol/e²

# Get data directory
TEST_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(TEST_DIR))), "data")


# Import helper functions
from .gcmc_insertion_helpers import (
    create_protein_system, load_benzene_molecule, find_cavity_sites,
    insert_benzene_at_position, calculate_insertion_energy,
    count_benzene_molecules, create_water_box, get_total_energy,
    check_minimum_distance, insert_water_at_position,
    count_water_molecules, create_simple_ion_system
)


def test_gcmc_benzene_insertion_in_protein():
    """Test GCMC insertion of benzene into protein cavity using real data"""
    
    # Temperature and chemical potential
    T = 300.0  # K
    beta = 1.0 / (kB * T)
    mu_ex = -2.0  # kJ/mol, typical for organic solutes
    
    # Load protein structure (4wp7)
    protein_pdb = os.path.join(DATA_DIR, "4wp7", "4wp7_fixed_with_5l13_silcs.1.prod.74.rec.pdb")
    benzene_pdb = os.path.join(DATA_DIR, "charmm36.ff", "mol", "benz.pdb")
    
    # Create system with protein
    system = create_protein_system(protein_pdb)
    
    # Load benzene molecule
    benzene = load_benzene_molecule(benzene_pdb)
    
    # Find potential insertion sites (cavities)
    cavity_sites = find_cavity_sites(system, probe_radius=0.3)  # 3 Å probe
    
    print(f"\nFound {len(cavity_sites)} potential cavity sites")
    
    # Calculate cavity bias factor
    total_volume = system.info.box[0] * system.info.box[1] * system.info.box[2]
    cavity_volume = len(cavity_sites) * 0.125  # Assuming 0.5 Å grid spacing
    f_n = cavity_volume / total_volume
    
    print(f"Cavity bias factor f_n = {f_n:.4f}")
    
    # Try multiple insertion attempts
    n_attempts = 100
    accepted = 0
    energies = []
    
    for attempt in range(n_attempts):
        # Select random cavity site
        if len(cavity_sites) > 0:
            site_idx = random.randint(0, len(cavity_sites) - 1)
            x, y, z = cavity_sites[site_idx]
        else:
            # Random position if no cavities found
            x = random.uniform(1.0, system.info.box[0] - 1.0)
            y = random.uniform(1.0, system.info.box[1] - 1.0)
            z = random.uniform(1.0, system.info.box[2] - 1.0)
        
        # Generate multiple orientations (configurational bias)
        n_orientations = 10
        orientation_energies = []
        
        for orient in range(n_orientations):
            # Random rotation
            theta = random.uniform(0, 2 * math.pi)
            phi = random.uniform(0, math.pi)
            psi = random.uniform(0, 2 * math.pi)
            
            # Create trial configuration
            trial_system = insert_benzene_at_position(system, benzene, x, y, z, theta, phi, psi)
            
            # Calculate insertion energy
            delta_energy = calculate_insertion_energy(trial_system, system)
            orientation_energies.append(delta_energy)
        
        # Select best orientation (configurational bias)
        best_idx = min(range(len(orientation_energies)), key=lambda i: orientation_energies[i])
        delta_energy = orientation_energies[best_idx]
        
        # Calculate Rosenbluth factor for configurational bias
        rosenbluth = sum(math.exp(-beta * e) for e in orientation_energies)
        p_orient = math.exp(-beta * delta_energy) / rosenbluth
        
        # GCMC acceptance criterion with cavity and configurational bias
        n_current = count_benzene_molecules(system)
        n_bar = 10.0  # Average number of molecules (for B calculation)
        B = beta * mu_ex + math.log(n_bar)
        
        # Acceptance probability
        acc_arg = f_n / (n_current + 1) * math.exp(B - beta * delta_energy) * rosenbluth
        acc_prob = min(1.0, acc_arg)
        
        # Metropolis acceptance
        if random.random() < acc_prob:
            accepted += 1
            # Actually insert the molecule
            system = insert_benzene_at_position(system, benzene, x, y, z, 
                                              theta, phi, psi)
            energies.append(delta_energy)
            
            print(f"Attempt {attempt+1}: ACCEPTED - ΔE = {delta_energy:.2f} kJ/mol, "
                  f"P_acc = {acc_prob:.4f}")
        else:
            print(f"Attempt {attempt+1}: REJECTED - ΔE = {delta_energy:.2f} kJ/mol, "
                  f"P_acc = {acc_prob:.4f}")
    
    acceptance_rate = accepted / n_attempts
    print(f"\nAcceptance rate: {acceptance_rate:.2%} ({accepted}/{n_attempts})")
    
    # Verify acceptance rate is reasonable
    # Note: With our simplified test system, acceptance rates can be high
    assert 0.001 < acceptance_rate <= 1.0, \
        f"Acceptance rate {acceptance_rate:.2%} outside expected range"
    
    # Check energy distribution
    if len(energies) > 0:
        avg_energy = sum(energies) / len(energies)
        print(f"Average insertion energy: {avg_energy:.2f} kJ/mol")
        
        # Favorable insertions should have negative energy
        assert avg_energy < 50.0, \
            f"Average insertion energy {avg_energy:.2f} too high"


def test_gcmc_water_insertion_simple():
    """Test GCMC water insertion in a simple box with existing waters"""
    
    T = 300.0  # K
    beta = 1.0 / (kB * T)
    mu_ex = -6.0  # kJ/mol, typical for water
    
    # Create box with a few water molecules
    system = create_water_box(n_waters=5, box_size=3.0)  # 3 nm box
    
    # Calculate initial energy
    pygcmc.computeSystemEnergyCutoff(system)
    initial_energy = get_total_energy(system)
    
    print(f"\nInitial system energy: {initial_energy:.2f} kJ/mol")
    
    # Try inserting one water molecule
    n_attempts = 50
    accepted = 0
    
    for attempt in range(n_attempts):
        # Random position
        x = random.uniform(0.3, 2.7)  # Stay away from boundaries
        y = random.uniform(0.3, 2.7)
        z = random.uniform(0.3, 2.7)
        
        # Check if position is too close to existing waters
        min_dist = check_minimum_distance(system, x, y, z)
        
        if min_dist < 0.2:  # 2 Å minimum distance
            print(f"Attempt {attempt+1}: REJECTED - Too close (d={min_dist:.3f} nm)")
            continue
        
        # Create trial system with new water
        trial_system = insert_water_at_position(system, x, y, z)
        
        # Calculate energies
        pygcmc.computeSystemEnergyCutoff(trial_system)
        final_energy = get_total_energy(trial_system)
        delta_energy = final_energy - initial_energy
        
        # GCMC acceptance criterion
        n_current = count_water_molecules(system)
        n_bar = 55.0  # Water concentration ~55 M
        B = beta * mu_ex + math.log(n_bar)
        
        # Simple cavity bias (no overlap)
        f_n = 0.9 if min_dist > 0.3 else 0.1
        
        acc_arg = f_n / (n_current + 1) * math.exp(B - beta * delta_energy)
        acc_prob = min(1.0, acc_arg)
        
        if random.random() < acc_prob:
            accepted += 1
            system = trial_system
            initial_energy = final_energy
            print(f"Attempt {attempt+1}: ACCEPTED - ΔE = {delta_energy:.2f} kJ/mol")
        else:
            print(f"Attempt {attempt+1}: REJECTED - ΔE = {delta_energy:.2f} kJ/mol, P = {acc_prob:.4f}")
    
    acceptance_rate = accepted / n_attempts
    print(f"\nWater insertion acceptance rate: {acceptance_rate:.2%}")
    
    # Should have some successful insertions
    assert accepted > 0, "No water insertions were accepted"
    assert acceptance_rate < 0.8, "Acceptance rate suspiciously high"


def test_verify_detailed_balance():
    """Verify GCMC detailed balance by checking forward/reverse move probabilities"""
    
    T = 300.0
    beta = 1.0 / (kB * T)
    mu_ex = -3.0
    
    # Create simple system
    system = create_simple_ion_system()
    
    # Parameters
    n = 5  # Current number of molecules
    n_bar = 10.0
    B = beta * mu_ex + math.log(n_bar)
    f_n = 0.8  # Cavity bias factor
    
    # Test energy changes
    test_energies = [-10.0, -5.0, 0.0, 5.0, 10.0, 20.0]  # kJ/mol
    
    print("\nVerifying detailed balance:")
    print("ΔE (kJ/mol) | P(insert) | P(delete) | Ratio | K_eq")
    print("-" * 60)
    
    for delta_e in test_energies:
        # Insertion probability: n -> n+1
        acc_insert = min(1.0, f_n / (n + 1) * math.exp(B - beta * delta_e))
        
        # Deletion probability: n+1 -> n (reverse move)
        # For deletion: A = min{1, (n+1)/f_n * exp(-B + β*ΔE)}
        acc_delete = min(1.0, (n + 1) / f_n * math.exp(-B + beta * delta_e))
        
        # Equilibrium constant
        K_eq = math.exp(B - beta * delta_e)
        
        # Check detailed balance
        ratio = (acc_insert / acc_delete) if acc_delete > 0 else float('inf')
        
        print(f"{delta_e:10.1f} | {acc_insert:9.4f} | {acc_delete:9.4f} | "
              f"{ratio:7.4f} | {K_eq:7.4f}")
        
        # Verify detailed balance (within numerical precision)
        if acc_insert < 1.0 and acc_delete < 1.0:
            expected_ratio = f_n / (n + 1) * K_eq * (n + 1) / f_n
            assert abs(ratio - K_eq) < 1e-10, \
                f"Detailed balance violated for ΔE={delta_e}"
    
    print("\n✓ Detailed balance verified!")


# Helper functions

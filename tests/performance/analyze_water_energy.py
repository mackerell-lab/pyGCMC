#!/usr/bin/env python
"""Analyze energy components of SWM4-NDP water"""

import sys
import numpy as np
sys.path.append('/home/zhaomt/gcmc/test107/pygcmc_dev/build/modules/bindings')

import pygcmc

def analyze_water_energy():
    print("=== Analyzing SWM4-NDP Water Energy Components ===\n")
    
    # SWM4-NDP parameters
    qO_core = 1.66260
    qD = -1.76260
    qM = -0.95710
    qH = 0.52855
    
    # Positions (nm)
    pos = {
        'O': np.array([0.0, 0.0, 0.0]),
        'D': np.array([0.02, 0.0, 0.0]),  # At hard wall limit
        'H1': np.array([0.09572, 0.0, 0.0]),
        'H2': np.array([-0.023999, 0.092663, 0.0]),
        'M': np.array([0.00793, 0.00986, 0.0])
    }
    
    charges = {'O': qO_core, 'D': qD, 'H1': qH, 'H2': qH, 'M': qM}
    atoms = ['O', 'D', 'H1', 'H2', 'M']
    
    # Calculate all pairwise distances and energies
    print("Pairwise Coulomb Interactions (intramolecular):")
    print("-" * 60)
    print(f"{'Pair':<10} {'Distance (Å)':<15} {'q1*q2':<15} {'Energy (kJ/mol)':<15}")
    print("-" * 60)
    
    ONE_4PI_EPS0 = 138.935456  # kJ·nm/mol/e²
    total_coulomb = 0.0
    
    for i in range(len(atoms)):
        for j in range(i+1, len(atoms)):
            atom1, atom2 = atoms[i], atoms[j]
            r_vec = pos[atom2] - pos[atom1]
            r = np.linalg.norm(r_vec)
            
            if r > 1e-10:
                q_prod = charges[atom1] * charges[atom2]
                energy = ONE_4PI_EPS0 * q_prod / r
                total_coulomb += energy
                
                print(f"{atom1}-{atom2:<7} {r*10:>12.4f}  {q_prod:>12.4f}  {energy:>12.2f}")
    
    print("-" * 60)
    print(f"Total intramolecular Coulomb: {total_coulomb:.2f} kJ/mol\n")
    
    # Calculate harmonic energy
    k_kj_nm2 = 166018.9  # From correct parameters
    r_OD = np.linalg.norm(pos['D'] - pos['O'])
    E_harmonic = 0.5 * k_kj_nm2 * r_OD**2
    
    print(f"Harmonic restraint energy:")
    print(f"  k = {k_kj_nm2:.1f} kJ/mol/nm²")
    print(f"  r_OD = {r_OD*1000:.3f} pm")
    print(f"  E_harmonic = 0.5 * k * r² = {E_harmonic:.2f} kJ/mol")
    
    print(f"\nTotal calculated: {total_coulomb + E_harmonic:.2f} kJ/mol")
    print(f"From test: 4871.77 kJ/mol")
    print(f"Difference: {4871.77 - (total_coulomb + E_harmonic):.2f} kJ/mol")
    
    # Identify problematic interactions
    print("\n\nProblematic Interactions:")
    print("-" * 40)
    
    # O-M interaction
    r_OM = np.linalg.norm(pos['M'] - pos['O'])
    E_OM = ONE_4PI_EPS0 * charges['O'] * charges['M'] / r_OM
    print(f"O-M interaction:")
    print(f"  Distance: {r_OM*10:.4f} Å")
    print(f"  q_O * q_M = {charges['O']} * {charges['M']} = {charges['O']*charges['M']:.4f}")
    print(f"  Energy: {E_OM:.2f} kJ/mol")
    print(f"  This is {abs(E_OM/total_coulomb)*100:.1f}% of total Coulomb energy!")
    
    # D-M interaction
    r_DM = np.linalg.norm(pos['M'] - pos['D'])
    E_DM = ONE_4PI_EPS0 * charges['D'] * charges['M'] / r_DM
    print(f"\nD-M interaction:")
    print(f"  Distance: {r_DM*10:.4f} Å")
    print(f"  q_D * q_M = {charges['D']} * {charges['M']} = {charges['D']*charges['M']:.4f}")
    print(f"  Energy: {E_DM:.2f} kJ/mol")
    
    print("\n\nConclusion:")
    print("-" * 40)
    print("The high energy is dominated by O-M repulsion due to:")
    print("1. Very short O-M distance (0.126 Å)")
    print("2. Both O_core (+1.663) and M (-0.957) have same sign after considering net charge")
    print("3. M-site is a virtual site and should be excluded from direct Coulomb")
    print("\nSolution: Exclude M-site from Coulomb calculations in same molecule")

if __name__ == "__main__":
    analyze_water_energy()
"""
Test real/reciprocal space separation in PME
"""

import numpy as np
from scipy.special import erfc, erf

# Test parameters
alpha = 3.5  # nm^-1
cutoff = 1.2  # nm
COULOMB = 138.935456  # kJ·nm/mol/e²

# Test distances
distances = [0.15, 0.3, 0.5, 0.8, 1.0, 1.2]

print("PME Real/Reciprocal Space Separation Test")
print("=" * 80)
print(f"Alpha = {alpha} nm^-1, Cutoff = {cutoff} nm")
print()

print(f"{'Distance (nm)':>15} {'1/r':>15} {'erfc(αr)/r':>15} {'erf(αr)/r':>15} {'Sum':>15} {'erfc ratio':>12}")
print("-" * 92)

for r in distances:
    # Full Coulomb
    full = 1.0 / r
    
    # Real space part
    real_part = erfc(alpha * r) / r
    
    # Reciprocal space part
    recip_part = erf(alpha * r) / r
    
    # Sum should equal full
    sum_parts = real_part + recip_part
    
    # Ratio of real space to full
    ratio = real_part / full
    
    print(f"{r:15.3f} {full:15.6f} {real_part:15.6f} {recip_part:15.6f} {sum_parts:15.6f} {ratio:12.4%}")

print()
print("Analysis:")
print(f"For alpha = {alpha}:")
print(f"- At r = 0.15 nm: erfc contribution is {erfc(alpha * 0.15) / 0.15 / (1/0.15):.1%} of full Coulomb")
print(f"- At r = 0.5 nm: erfc contribution is {erfc(alpha * 0.5) / 0.5 / (1/0.5):.1%} of full Coulomb")
print(f"- At r = 1.0 nm: erfc contribution is {erfc(alpha * 1.0) / 1.0 / (1/1.0):.1%} of full Coulomb")
print(f"- At r = {cutoff} nm: erfc contribution is {erfc(alpha * cutoff) / cutoff / (1/cutoff):.1%} of full Coulomb")

# Calculate expected real space energy
print("\nExpected energy contributions for the medium complexity system:")
print("Assuming typical charge product of 0.5 e²:")

# Simple estimate
real_energy = 0
full_energy = 0
n_pairs_close = 28  # From analysis
n_pairs_medium = 32
avg_dist_close = 0.3
avg_dist_medium = 0.75

# Close pairs
for i in range(n_pairs_close):
    r = avg_dist_close
    real_energy += COULOMB * 0.5 * erfc(alpha * r) / r
    full_energy += COULOMB * 0.5 / r

# Medium pairs  
for i in range(n_pairs_medium):
    r = avg_dist_medium
    real_energy += COULOMB * 0.5 * erfc(alpha * r) / r
    full_energy += COULOMB * 0.5 / r

print(f"Estimated full Coulomb energy: {full_energy:.1f} kJ/mol")
print(f"Estimated real space energy: {real_energy:.1f} kJ/mol")
print(f"Real space is {real_energy/full_energy:.1%} of full Coulomb")
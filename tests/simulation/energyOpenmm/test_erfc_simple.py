"""
Simple test to check erfc calculation implementation
"""

import math

# Test the erfc calculations
print("Testing erfc calculations for PME")
print("=" * 60)

alpha = 3.5
test_distances = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0]

print("r (nm) | α*r  | erfc(α*r) | erfc(α*r)/r | 1/r | Ratio")
print("-" * 60)

for r in test_distances:
    alpha_r = alpha * r
    erfc_val = math.erfc(alpha_r)
    erfc_over_r = erfc_val / r
    one_over_r = 1.0 / r
    ratio = one_over_r / erfc_over_r
    
    print(f"{r:6.1f} | {alpha_r:4.2f} | {erfc_val:9.6f} | {erfc_over_r:11.6f} | {one_over_r:5.1f} | {ratio:7.1f}x")

print()
print("Key insight:")
print("- The erfc term reduces the Coulomb interaction significantly")
print("- For r=0.5nm: 1/r is 75x larger than erfc(α*r)/r")
print("- If PyGCMC calculates 1/r instead of erfc(α*r)/r, it would be ~75x too large")
print()
print("Looking at test failures:")
print("- PME real space is 139-247x too large")
print("- This suggests PyGCMC might be calculating something closer to 1/r")
print("- Or the erfc implementation might be returning incorrect values")
"""
Analyze PME real space error patterns from test failures
"""

import math

# Data from test failures
test_data = [
    # Test name, Expected real space, Actual real space, Factor
    ("test_pgp_simple", -21.42, -2976.53, 139),
    ("test_pme_vs_ewald", 685.23, -169054.19, 247),  # sign error too
    ("test_pme_medium_complexity", None, -5642.33, None),  # dominating total
]

print("PME Real Space Error Analysis")
print("=" * 60)

# For α=3.5, analyze reduction factors at different distances
alpha = 3.5
distances = [0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0]

print(f"\nFor α={alpha}, erfc reduction factors:")
print("r (nm) | erfc(αr) | erfc(αr)/r | 1/r | Factor")
print("-" * 50)

for r in distances:
    erfc_val = math.erfc(alpha * r)
    erfc_over_r = erfc_val / r
    one_over_r = 1.0 / r
    factor = one_over_r / erfc_over_r
    print(f"{r:6.2f} | {erfc_val:8.5f} | {erfc_over_r:10.5f} | {one_over_r:5.1f} | {factor:6.0f}x")

print("\nObservations:")
print("1. Test failures show 139-247x larger real space energies")
print("2. This is between the reduction factors at r=0.5nm (75x) and r=0.6nm (336x)")
print("3. Suggests PyGCMC might be:")
print("   - Using wrong α value in erfcApprox")
print("   - Not applying erfc correctly")
print("   - Have a systematic error in the calculation")

# Check what alpha value would give 139x at r=0.5
print("\nReverse engineering:")
target_factor = 139
r = 0.5
# We want: (1/r) / (erfc(α*r)/r) = 139
# So: 1/erfc(α*r) = 139
# erfc(α*r) = 1/139 = 0.00719
# Need to find α such that erfc(α*0.5) = 0.00719

# Simple bisection method to replace scipy.optimize.fsolve
def bisect_find_alpha(r, target_factor, tol=1e-6):
    """Find alpha such that erfc(alpha * r) = 1/target_factor"""
    # Initial bounds
    a, b = 0.1, 10.0
    target = 1.0 / target_factor
    
    while b - a > tol:
        mid = (a + b) / 2
        val = math.erfc(mid * r)
        if val > target:
            a = mid
        else:
            b = mid
    return (a + b) / 2

alpha_139 = bisect_find_alpha(r, target_factor)
print(f"\nTo get 139x factor at r=0.5nm, α would need to be: {alpha_139:.2f}")
print(f"Verification: erfc({alpha_139:.2f} * 0.5) = {math.erfc(alpha_139 * 0.5):.5f}")

# Check for 247x
target_factor = 247
alpha_247 = bisect_find_alpha(r, target_factor)
print(f"\nTo get 247x factor at r=0.5nm, α would need to be: {alpha_247:.2f}")
print(f"Verification: erfc({alpha_247:.2f} * 0.5) = {math.erfc(alpha_247 * 0.5):.5f}")
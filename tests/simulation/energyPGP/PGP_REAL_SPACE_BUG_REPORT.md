# PGP Real-Space Calculation Bug Report

## Summary
The PGP (Precomputed Grid Potential) real-space energy calculation always returns 0.0, which is a critical bug that prevents PGP from correctly calculating short-range electrostatic interactions.

## Root Cause
The bug is caused by a design flaw in how the PMEParams class methods are implemented:

1. `PGPParams` inherits from `PMEParams`
2. When `pgp_params.erfcApprox(r)` is called in `PGPReal.cpp`, it uses the inherited method from PMEParams
3. The `erfcApprox()` method in PMEParams calls the global function `erfcApproximate(r)`
4. The global `erfcApproximate(r)` function uses the global `pme_params` instance instead of the calling instance
5. Since `pme_params.cutoff` might not be initialized when using PGP, the function returns 0.0

## Code Analysis

### Problem Location
File: `/home/zhaomt/gcmc/test107/pygcmc_dev/src/platform/cpu/energy/pme/PMECore.hpp`
```cpp
struct PMEParams {
    // ...
    double erfcApprox(double r) const {
        // Redirect to standalone function
        return erfcApproximate(r);  // BUG: Uses global pme_params
    }
};
```

File: `/home/zhaomt/gcmc/test107/pygcmc_dev/src/platform/cpu/energy/pme/PMEConfig.cpp`
```cpp
double erfcApproximate(double r) {
    if (r >= pme_params.cutoff) return 0.0;  // BUG: Uses pme_params instead of pgp_params
    // ... rest of implementation uses pme_params.erfcTable
}
```

### Impact
- PGP real-space energy is always 0.0
- Only affects inter-residue interactions (intra-residue loop issue is separate)
- PGP cannot correctly calculate short-range electrostatic interactions
- Tests comparing PGP with PME show significant discrepancies

### Test Results
1. Two opposite charges at 0.3 nm:
   - Expected real-space: -183.46 kJ/mol
   - PME real-space: -183.46 kJ/mol
   - PGP real-space: 0.0 kJ/mol ❌

2. Two opposite charges at 0.01 nm:
   - Expected real-space: -13580.04 kJ/mol
   - PME real-space: -13580.06 kJ/mol
   - PGP real-space: 0.0 kJ/mol ❌

## Proposed Solutions

### Solution 1: Override Methods in PGPParams
```cpp
struct PGPParams : public PMEParams {
    // Override the erfcApprox method
    double erfcApprox(double r) const {
        if (r >= cutoff) return 0.0;
        double x = r * erfcDXInv;
        int index = static_cast<int>(x);
        if (index >= static_cast<int>(erfcTable.size()) - 1) {
            return erfcTable.back();
        }
        double fraction = x - index;
        return erfcTable[index] + fraction * (erfcTable[index+1] - erfcTable[index]);
    }
};
```

### Solution 2: Make Global Functions Take Parameters
```cpp
double erfcApproximate(double r, const PMEParams& params) {
    if (r >= params.cutoff) return 0.0;
    // ... use params instead of pme_params
}
```

### Solution 3: Use Member Variables Directly
Remove the redirection to global functions and implement the methods directly in PMEParams using member variables.

## Additional Issue: Intra-Residue Interactions
The PGP real-space calculation loop only calculates inter-residue interactions (r2 = r1 + 1), missing intra-residue interactions. This is a separate issue that also needs to be addressed.

## Verification
Created multiple test cases that demonstrate:
1. PGP real-space returns 0.0 for all configurations
2. The issue persists regardless of initialization order
3. The issue is not related to cutoff distance
4. PME calculations work correctly with the same configurations
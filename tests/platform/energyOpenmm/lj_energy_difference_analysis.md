# LJ Energy Difference Analysis: PyGCMC vs OpenMM

## Summary

The 15% difference in LJ energy between PyGCMC and OpenMM is caused by **two distinct issues**:

1. **PyGCMC skips intramolecular LJ interactions**
2. **PyGCMC double-counts intermolecular LJ interactions**

## Issue 1: Missing Intramolecular Interactions

PyGCMC does not calculate LJ interactions between atoms within the same residue.

### Evidence from Code

In `src/platform/cpu/energy/common/EnergyDirectCore.cpp`, line 115:
```cpp
if (!residues[j].active || j == residue_idx) continue;
```

This explicitly skips when `j == residue_idx`, meaning atoms in the same residue don't interact.

### Test Results

Single molecule with 3 atoms:
- PyGCMC: 0.0 kJ/mol (no interactions calculated)
- OpenMM: -1.246 kJ/mol (includes all 3 pair interactions)
- Manual calculation: -1.246 kJ/mol

## Issue 2: Double-Counting Intermolecular Interactions

PyGCMC counts each intermolecular interaction twice - once for each residue involved.

### Evidence from Tests

Two-atom system (separate residues):
- PyGCMC: 2x the expected energy
- OpenMM: Correct single-counted energy
- Ratio consistently ~2.0 across different test cases

### How Double-Counting Occurs

When calculating energy for residue i:
1. PyGCMC loops through all atoms in residue i
2. For each atom, it calculates interactions with all atoms in other residues
3. The energy is added to residue i's total

When calculating energy for residue j:
1. The same interaction is calculated again (from j's perspective)
2. The energy is added to residue j's total

Total system energy = sum of all residue energies = 2x the actual energy

## Combined Effect

For a typical system with molecules containing multiple atoms:

1. **Intramolecular LJ missing**: Reduces energy (makes it less negative)
2. **Intermolecular LJ double-counted**: Increases energy magnitude (makes it more negative)

The net effect depends on the ratio of intramolecular to intermolecular interactions:
- Systems with many small molecules: Double-counting dominates → PyGCMC energy too negative
- Systems with large molecules: Missing intramolecular dominates → PyGCMC energy not negative enough

## OpenMM's Approach

OpenMM by default:
1. Calculates all non-bonded interactions (including intramolecular)
2. Counts each interaction exactly once
3. Uses exclusions/exceptions for bonded atoms (1-2, 1-3, scaled 1-4)

## Implications for GCMC

This affects:
- Absolute energies (incorrect by ~15% or more)
- Energy differences for insertion/deletion moves
- Acceptance probabilities
- Equilibrium distributions

## Recommendations

To fix PyGCMC:
1. Include intramolecular LJ calculations (remove the `j == residue_idx` skip)
2. Fix double-counting by either:
   - Calculating each pair only once (i < j approach)
   - Dividing the final energy by 2
   - Splitting pair energy equally between residues (current approach but ensure not double-counted)
#!/usr/bin/env python
"""Test to understand virtual site exclusions in SWM4-NDP"""

import numpy as np

def test_virtual_site_exclusion():
    print("=== Understanding Virtual Site Exclusions ===\n")
    
    print("OpenMM SWM4-NDP Implementation:")
    print("-" * 50)
    print("From TestDrudeSCFIntegrator.h line 79-81:")
    print("for (int j = 0; j < 5; j++)")
    print("    for (int k = 0; k < j; k++)")
    print("        nonbonded->addException(startIndex+j, startIndex+k, 0, 1, 0);")
    print("\nThis adds exclusions for ALL pairs within the molecule:")
    print("- O-D, O-H1, O-H2, O-M")
    print("- D-H1, D-H2, D-M") 
    print("- H1-H2, H1-M")
    print("- H2-M")
    print("\nTotal: 10 exclusions (all possible pairs)\n")
    
    print("Key Insight:")
    print("-" * 50)
    print("1. ALL intramolecular nonbonded interactions are excluded")
    print("2. This includes virtual site M interactions") 
    print("3. Only Drude harmonic restraint remains intramolecularly")
    print("4. This explains why single water should have low energy\n")
    
    print("For GCMC with rigid molecules:")
    print("-" * 50)
    print("- Intramolecular geometry is fixed")
    print("- Only intermolecular interactions matter")
    print("- Drude particles adjust via SCF")
    print("- Virtual sites participate in intermolecular Coulomb")
    print("- But NOT in intramolecular Coulomb (excluded)\n")
    
    # Calculate expected single water energy
    print("Expected single water energy:")
    print("-" * 30)
    k_kj_nm2 = 166018.9
    r_OD = 0.02  # nm (at hard wall)
    E_harmonic = 0.5 * k_kj_nm2 * r_OD**2
    print(f"Only harmonic restraint: {E_harmonic:.2f} kJ/mol")
    print("(No intramolecular Coulomb due to exclusions)")
    
    print("\n\nConclusion:")
    print("=" * 50)
    print("Our implementation needs to:")
    print("1. Exclude ALL intramolecular nonbonded interactions")
    print("2. Include virtual sites in these exclusions")
    print("3. Only calculate Drude harmonic restraint intramolecularly")
    print("4. Virtual sites DO participate in INTERmolecular interactions")

if __name__ == "__main__":
    test_virtual_site_exclusion()
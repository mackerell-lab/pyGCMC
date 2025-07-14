# tests/simulation/energyPGP/test_pgp_real_space_debug.py
"""
Debug test to understand why PGP real-space is zero.
"""

import pytest
import math
import pygcmc
from pygcmc import MCAtom, MCResidue, MCState
from pygcmc import setPGPParameters, initializePMEParameters, precomputeGridPotential
from pygcmc import computeSystemEnergyPGP, computeSystemEnergyPME


def test_pgp_real_space_debug():
    """Debug test with two atoms in DIFFERENT residues."""
    
    print("\n=== Debug: PGP Real-Space with Two Residues ===")
    
    # Create a simple system
    state = MCState()
    state.info.box = [5.0, 5.0, 5.0]
    state.info.cutoff = 1.2
    
    # Force field
    state.forcefield.numTotalTypes = 1
    state.forcefield.numMovementTypes = 1
    state.forcefield.ljEps = [0.0]
    state.forcefield.ljSigma = [0.3]
    
    # Create two atoms in DIFFERENT residues
    atoms = []
    
    atom1 = MCAtom()
    atom1.x = 2.5
    atom1.y = 2.5
    atom1.z = 2.5
    atom1.charge = 1.0
    atom1.type = 0
    atoms.append(atom1)
    
    atom2 = MCAtom()
    atom2.x = 3.0  # 0.5 nm away
    atom2.y = 2.5
    atom2.z = 2.5
    atom2.charge = -1.0
    atom2.type = 0
    atoms.append(atom2)
    
    state.atoms = atoms
    state.activeAtomCount = 2
    
    # Create TWO separate residues
    residues = []
    
    res1 = MCResidue()
    res1.active = True
    res1.fixed = False
    res1.atomStart = 0
    res1.atomCount = 1
    res1.type = 0
    residues.append(res1)
    
    res2 = MCResidue()
    res2.active = True
    res2.fixed = False
    res2.atomStart = 1
    res2.atomCount = 1
    res2.type = 0
    residues.append(res2)
    
    state.residues = residues
    state.activeResidueCount = 2
    
    print("System: Two atoms in DIFFERENT residues")
    print(f"Atom 0: pos=({atom1.x}, {atom1.y}, {atom1.z}), charge={atom1.charge}")
    print(f"Atom 1: pos=({atom2.x}, {atom2.y}, {atom2.z}), charge={atom2.charge}")
    print(f"Distance: 0.5 nm")
    
    # Initialize PGP with debug mode
    alpha = 2.0
    setPGPParameters(
        alpha=alpha,
        meshSize=[32, 32, 32],
        potential_cutoff=state.info.cutoff,
        potentialGridSize=[32, 32, 32],
        splineOrder=4,
        tolerance=1e-5
    )
    
    initializePMEParameters(state.info.cutoff, state.info.box, alpha)
    precomputeGridPotential(state, fixed_only=True)
    
    # Calculate with PGP
    print("\n--- PGP Calculation ---")
    computeSystemEnergyPGP(state)
    
    pgp_real = state.ewald_energy.get('real_space', 0.0)
    pgp_recip = state.ewald_energy.get('reciprocal', 0.0)
    pgp_self = state.ewald_energy.get('self', 0.0)
    pgp_total = state.ewald_energy.get('total', 0.0)
    
    print(f"PGP Real-space: {pgp_real:.6f} kJ/mol")
    print(f"PGP Reciprocal: {pgp_recip:.6f} kJ/mol")
    print(f"PGP Self: {pgp_self:.6f} kJ/mol")
    print(f"PGP Total: {pgp_total:.6f} kJ/mol")
    
    # Also test with PME for comparison
    print("\n--- PME Calculation ---")
    computeSystemEnergyPME(state)
    
    pme_real = state.ewald_energy.get('real_space', 0.0)
    pme_recip = state.ewald_energy.get('reciprocal', 0.0)
    pme_self = state.ewald_energy.get('self', 0.0)
    pme_total = state.ewald_energy.get('total', 0.0)
    
    print(f"PME Real-space: {pme_real:.6f} kJ/mol")
    print(f"PME Reciprocal: {pme_recip:.6f} kJ/mol")
    print(f"PME Self: {pme_self:.6f} kJ/mol")
    print(f"PME Total: {pme_total:.6f} kJ/mol")
    
    # Manual calculation
    r = 0.5
    kC = 138.935456
    erfc_val = math.erfc(alpha * r)
    expected_real = -erfc_val / r * kC
    print(f"\nExpected real-space: {expected_real:.6f} kJ/mol")
    
    # Check if PGP real-space matches PME
    if abs(pgp_real) < 0.1:
        print("\n❌ PGP real-space is near zero!")
    else:
        print("\n✅ PGP real-space is non-zero")
        
    if abs(pgp_real - pme_real) < 1.0:
        print("✅ PGP and PME real-space match")
    else:
        print("❌ PGP and PME real-space differ significantly")


def test_pgp_with_water_molecule():
    """Test with a water molecule (3 atoms in same residue)."""
    
    print("\n=== Test: PGP with Water Molecule ===")
    
    state = MCState()
    state.info.box = [5.0, 5.0, 5.0]
    state.info.cutoff = 1.2
    
    # Force field
    state.forcefield.numTotalTypes = 2
    state.forcefield.numMovementTypes = 2
    state.forcefield.ljEps = [0.6364, 0.0, 0.0, 0.0]  # O-O only
    state.forcefield.ljSigma = [0.3166, 0.0, 0.0, 0.0]
    
    # Create water molecule
    atoms = []
    
    # Oxygen
    o = MCAtom()
    o.x = 2.5
    o.y = 2.5
    o.z = 2.5
    o.charge = -0.834
    o.type = 0
    atoms.append(o)
    
    # Hydrogen 1
    h1 = MCAtom()
    h1.x = 2.5957
    h1.y = 2.5
    h1.z = 2.5
    h1.charge = 0.417
    h1.type = 1
    atoms.append(h1)
    
    # Hydrogen 2
    h2 = MCAtom()
    h2.x = 2.4243
    h2.y = 2.5587
    h2.z = 2.5
    h2.charge = 0.417
    h2.type = 1
    atoms.append(h2)
    
    state.atoms = atoms
    state.activeAtomCount = 3
    
    # Single residue for water
    residues = []
    res = MCResidue()
    res.active = True
    res.fixed = False
    res.atomStart = 0
    res.atomCount = 3
    res.type = 0
    residues.append(res)
    
    state.residues = residues
    state.activeResidueCount = 1
    
    print("Water molecule - all atoms in SAME residue")
    print("O-H1 distance: 0.0957 nm")
    print("O-H2 distance: 0.0957 nm")
    print("H1-H2 distance: ~0.1515 nm")
    
    # Initialize
    alpha = 2.0
    setPGPParameters(
        alpha=alpha,
        meshSize=[32, 32, 32],
        potential_cutoff=state.info.cutoff,
        potentialGridSize=[32, 32, 32],
        splineOrder=4,
        tolerance=1e-5
    )
    
    initializePMEParameters(state.info.cutoff, state.info.box, alpha)
    precomputeGridPotential(state, fixed_only=True)
    
    # Calculate
    computeSystemEnergyPGP(state)
    
    real_space = state.ewald_energy.get('real_space', 0.0)
    self_energy = state.ewald_energy.get('self', 0.0)
    
    print(f"\nPGP Real-space: {real_space:.6f} kJ/mol")
    print(f"PGP Self energy: {self_energy:.6f} kJ/mol")
    
    # For a single water molecule, real-space should be zero
    # because intra-residue interactions are not calculated in the loop
    if abs(real_space) < 0.1:
        print("\n⚠️  Real-space is near zero for single water molecule")
        print("This confirms that intra-residue interactions are NOT calculated in PGP real-space!")
    else:
        print("\n✅ Real-space is non-zero")


if __name__ == "__main__":
    test_pgp_real_space_debug()
    test_pgp_with_water_molecule()
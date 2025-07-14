# tests/simulation/energyPGP/test_pgp_real_minimal.py
"""
Minimal test to trace why PGP real-space is zero.
"""

import pytest
import math
import pygcmc
from . import pgp_wrapper
from .pgp_wrapper import initializePMEParameters, setPGPParameters, precomputeGridPotential
from .pgp_wrapper import computeSystemEnergyPGP

import os
from pygcmc import MCAtom, MCResidue, MCState

def test_pgp_real_space_minimal():
    """Minimal test with debug output enabled."""
    
    print("\n=== Minimal PGP Real-Space Test ===")
    
    # Enable debug mode
    os.environ['PYGCMC_DEBUG'] = '1'
    
    # Create minimal system
    state = MCState()
    state.info.box = [5.0, 5.0, 5.0]
    state.info.cutoff = 2.0  # Larger cutoff to ensure atoms are within range
    
    # Minimal force field
    state.forcefield.numTotalTypes = 1
    state.forcefield.numMovementTypes = 1
    state.forcefield.ljEps = [0.0]
    state.forcefield.ljSigma = [0.3]
    
    # Two atoms, different residues, close together
    atoms = []
    
    # Atom 1 - residue 0
    atom1 = MCAtom()
    atom1.x = 2.5
    atom1.y = 2.5
    atom1.z = 2.5
    atom1.charge = 1.0  # +1 charge
    atom1.type = 0
    atoms.append(atom1)
    
    # Atom 2 - residue 1
    atom2 = MCAtom()
    atom2.x = 2.8  # Only 0.3 nm away
    atom2.y = 2.5
    atom2.z = 2.5
    atom2.charge = -1.0  # -1 charge
    atom2.type = 0
    atoms.append(atom2)
    
    state.atoms = atoms
    state.activeAtomCount = 2
    
    # Two residues
    residues = []
    
    # Residue 0
    res0 = MCResidue()
    res0.active = True
    res0.fixed = False  # Movable
    res0.atomStart = 0
    res0.atomCount = 1
    res0.type = 0
    res0.energy_vdw = 0.0
    res0.energy_elec = 0.0
    residues.append(res0)
    
    # Residue 1
    res1 = MCResidue()
    res1.active = True
    res1.fixed = False  # Movable
    res1.atomStart = 1
    res1.atomCount = 1
    res1.type = 0
    res1.energy_vdw = 0.0
    res1.energy_elec = 0.0
    residues.append(res1)
    
    state.residues = residues
    state.activeResidueCount = 2
    
    # Print system info
    print(f"System: 2 atoms, 2 residues")
    print(f"Distance: 0.3 nm (cutoff = {state.info.cutoff} nm)")
    print(f"Charges: +1 and -1")
    print(f"Both residues are movable (fixed=False)")
    
    # Initialize PGP
    alpha = 2.0
    setPGPParameters(
        alpha=alpha,
        meshSize=[16, 16, 16],  # Smaller mesh for debugging
        potential_cutoff=state.info.cutoff,
        potentialGridSize=[16, 16, 16],
        splineOrder=4,
        tolerance=1e-5
    )
    
    # Also need to initialize PME parameters
    initializePMEParameters(state.info.cutoff, state.info.box, alpha)
    
    print("\nPrecomputing grid (should be empty since no fixed atoms)...")
    precomputeGridPotential(state)
    computeSystemEnergyPGP(state)
    real_space = state.ewald_energy.get("real_space", 0.0)
    reciprocal = state.ewald_energy.get('reciprocal')
    
    print(f"\nResults:")
    print(f"Real-space energy: {real_space:.6f} kJ/mol")
    print(f"Reciprocal energy: {reciprocal:.6f} kJ/mol")
    
    # Check residue energies
    print(f"\nResidue energies:")
    for i, res in enumerate(state.residues):
        print(f"  Residue {i}: elec={res.energy_elec:.6f}, vdw={res.energy_vdw:.6f}")
    
    # Manual calculation
    r = 0.3
    kC = 138.935456
    erfc_val = math.erfc(alpha * r)
    expected = -erfc_val / r * kC
    
    print(f"\nExpected real-space: {expected:.6f} kJ/mol")
    print(f"erfc({alpha}*{r}) = {erfc_val:.6f}")
    
    if abs(real_space) < 0.1:
        print("\n❌ PROBLEM: Real-space energy is zero!")
        print("This suggests the real-space calculation is not being performed correctly.")
    else:
        print("\n✅ Real-space energy is non-zero")

if __name__ == "__main__":
    test_pgp_real_space_minimal()

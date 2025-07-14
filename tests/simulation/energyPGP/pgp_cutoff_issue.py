# tests/simulation/energyPGP/test_pgp_cutoff_issue.py
"""
Test to verify if PGP real-space issue is due to cutoff check.
"""

import pytest
import math
import pygcmc
from . import pgp_wrapper
from .pgp_wrapper import initializePMEParameters, setPGPParameters, precomputeGridPotential, computeSystemEnergyPGP
from .pgp_wrapper import computeSystemEnergyPGP, computeSystemEnergyPME
from pygcmc import MCAtom, MCResidue, MCState

def test_pgp_cutoff_issue():
    """Test if PGP real-space issue is due to cutoff check in erfcApproximate."""
    
    print("\n=== Test: PGP Cutoff Issue ===")
    
    # Create system with atoms VERY close together
    state = MCState()
    state.info.box = [5.0, 5.0, 5.0]
    state.info.cutoff = 2.0  # Large cutoff
    
    # Force field
    state.forcefield.numTotalTypes = 1
    state.forcefield.numMovementTypes = 1
    state.forcefield.ljEps = [0.0]
    state.forcefield.ljSigma = [0.3]
    
    # Two atoms VERY close
    atoms = []
    
    atom1 = MCAtom()
    atom1.x = 2.5
    atom1.y = 2.5
    atom1.z = 2.5
    atom1.charge = 1.0
    atom1.type = 0
    atoms.append(atom1)
    
    atom2 = MCAtom()
    atom2.x = 2.51  # Only 0.01 nm away - well within cutoff
    atom2.y = 2.5
    atom2.z = 2.5
    atom2.charge = -1.0
    atom2.type = 0
    atoms.append(atom2)
    
    state.atoms = atoms
    state.activeAtomCount = 2
    
    # Two residues
    residues = []
    for i in range(2):
        res = MCResidue()
        res.active = True
        res.fixed = False
        res.atomStart = i
        res.atomCount = 1
        res.type = 0
        residues.append(res)
    
    state.residues = residues
    state.activeResidueCount = 2
    
    print(f"System: +1 and -1 charges at 0.01 nm distance")
    print(f"Cutoff: {state.info.cutoff} nm")
    
    # Initialize
    alpha = 2.0
    mesh_size = [32, 32, 32]
    
    # Initialize PME first
    initializePMEParameters(state.info.cutoff, state.info.box, alpha)
    
    # Then PGP
    setPGPParameters(alpha, mesh_size, state.info.cutoff, mesh_size, 4, 1e-6)
    
    precomputeGridPotential(state)
    
    # Calculate with PME
    print("\nPME calculation:")
    computeSystemEnergyPME(state)
    pme_real = state.ewald_energy.get('real_space', 0.0)
    print(f"PME real-space: {pme_real:.6f} kJ/mol")
    
    # Calculate with PGP
    print("\nPGP calculation:")
    computeSystemEnergyPGP(state)
    pgp_real = state.ewald_energy.get('real_space', 0.0)
    print(f"PGP real-space: {pgp_real:.6f} kJ/mol")
    
    # Manual calculation
    r = 0.01
    kC = 138.935456
    erfc_val = math.erfc(alpha * r)
    expected = -erfc_val / r * kC
    print(f"\nExpected real-space: {expected:.6f} kJ/mol")
    print(f"erfc({alpha}*{r}) = {erfc_val:.6f}")
    
    if abs(pgp_real) < 0.1:
        print("\n❌ PGP real-space is still zero even with very close atoms!")
        print("This confirms the issue is in the erfcApproximate function")
        print("which uses pme_params instead of pgp_params")
    else:
        print("\n✅ PGP real-space is non-zero")

def test_pgp_with_different_cutoffs():
    """Test PGP with different cutoff values."""
    
    print("\n=== Test: PGP with Different Cutoffs ===")
    
    # Create system
    def create_system(cutoff):
        state = MCState()
        state.info.box = [5.0, 5.0, 5.0]
        state.info.cutoff = cutoff
        
        state.forcefield.numTotalTypes = 1
        state.forcefield.numMovementTypes = 1
        state.forcefield.ljEps = [0.0]
        state.forcefield.ljSigma = [0.3]
        
        atoms = []
        
        atom1 = MCAtom()
        atom1.x = 2.5
        atom1.y = 2.5
        atom1.z = 2.5
        atom1.charge = 1.0
        atom1.type = 0
        atoms.append(atom1)
        
        atom2 = MCAtom()
        atom2.x = 2.8  # 0.3 nm away
        atom2.y = 2.5
        atom2.z = 2.5
        atom2.charge = -1.0
        atom2.type = 0
        atoms.append(atom2)
        
        state.atoms = atoms
        state.activeAtomCount = 2
        
        residues = []
        for i in range(2):
            res = MCResidue()
            res.active = True
            res.fixed = False
            res.atomStart = i
            res.atomCount = 1
            res.type = 0
            residues.append(res)
        
        state.residues = residues
        state.activeResidueCount = 2
        
        return state
    
    alpha = 2.0
    mesh_size = [32, 32, 32]
    distance = 0.3  # nm
    
    cutoffs = [0.5, 1.0, 1.5, 2.0]
    
    for cutoff in cutoffs:
        print(f"\n--- Cutoff = {cutoff} nm (distance = {distance} nm) ---")
        
        state = create_system(cutoff)
        
        # Initialize
        initializePMEParameters(cutoff, state.info.box, alpha)
        setPGPParameters(alpha, mesh_size, state.info.cutoff, mesh_size, 4, 1e-6)
        precomputeGridPotential(state)
        
        # Calculate
        computeSystemEnergyPGP(state)
        real_space = state.ewald_energy.get('real_space', 0.0)
        
        print(f"PGP real-space: {real_space:.6f} kJ/mol")
        
        if distance < cutoff and abs(real_space) < 0.1:
            print("❌ Zero despite distance < cutoff")
        elif distance >= cutoff and abs(real_space) < 0.1:
            print("✓ Zero as expected (distance >= cutoff)")
        else:
            print("✅ Non-zero")

if __name__ == "__main__":
    test_pgp_cutoff_issue()
    test_pgp_with_different_cutoffs()

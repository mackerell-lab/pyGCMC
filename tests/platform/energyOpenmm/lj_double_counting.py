"""
Test to confirm LJ double-counting issue in PyGCMC

This test creates simple systems to verify if PyGCMC is double-counting
LJ interactions.
"""

import pytest
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pygcmc
from pygcmc import MCState, MCAtom, MCResidue, MCForceField
from pygcmc import computeSystemEnergyCutoff

# Only import OpenMM if available
try:
    from openmm import *
    from openmm.app import *
    from openmm.unit import *
    OPENMM_AVAILABLE = True
except ImportError:
    OPENMM_AVAILABLE = False


def test_single_pair_lj():
    """Test LJ energy for a single pair of atoms"""
    
    print("\nSingle pair LJ test:")
    
    state = MCState()
    state.info.box = [10.0, 10.0, 10.0]  # Large box to avoid PBC
    state.info.cutoff = 5.0
    
    # Simple force field
    ff = MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.ljEps = [1.0]     # kJ/mol
    ff.ljSigma = [0.35]  # nm
    state.forcefield = ff
    
    # Two atoms at 0.5 nm distance
    atoms = []
    residues = []
    
    atom1 = MCAtom()
    atom1.x, atom1.y, atom1.z = 5.0, 5.0, 5.0
    atom1.charge = 0.0
    atom1.type = 0
    atoms.append(atom1)
    
    atom2 = MCAtom()
    atom2.x, atom2.y, atom2.z = 5.5, 5.0, 5.0  # 0.5 nm away
    atom2.charge = 0.0
    atom2.type = 0
    atoms.append(atom2)
    
    # Each atom in its own residue
    for i in range(2):
        res = MCResidue()
        res.active = True
        res.fixed = True
        res.atomStart = i
        res.atomCount = 1
        res.type = 0
        residues.append(res)
    
    state.atoms = atoms
    state.activeAtomCount = 2
    state.residues = residues
    state.activeResidueCount = 2
    
    # Calculate energy
    computeSystemEnergyCutoff(state)
    
    # Get individual residue energies
    res1_energy = state.residues[0].energy_vdw
    res2_energy = state.residues[1].energy_vdw
    total_energy = res1_energy + res2_energy
    
    # Calculate expected LJ energy
    r = 0.5
    sigma = 0.35
    epsilon = 1.0
    sr6 = (sigma/r)**6
    expected_single = 4.0 * epsilon * (sr6*sr6 - sr6)
    
    print(f"  Residue 1 energy: {res1_energy:.6f} kJ/mol")
    print(f"  Residue 2 energy: {res2_energy:.6f} kJ/mol")
    print(f"  Total energy: {total_energy:.6f} kJ/mol")
    print(f"  Expected (single count): {expected_single:.6f} kJ/mol")
    print(f"  Expected (double count): {2*expected_single:.6f} kJ/mol")
    
    # Check if energy is split between residues
    if abs(res1_energy - res2_energy) < 1e-6:
        print("\n  ✓ Energy is split equally between residues")
        print(f"  Each residue gets: {res1_energy:.6f} kJ/mol")
        if abs(total_energy - 2*expected_single) < 1e-3:
            print("  ⚠️  Total matches double-counted energy!")
        elif abs(total_energy - expected_single) < 1e-3:
            print("  ✓ Total matches single-counted energy")
    
    # Also check total system energy
    print(f"\n  Checking system totals:")
    total_vdw = 0.0
    for res in state.residues:
        total_vdw += res.energy_vdw
    print(f"  Sum of residue VDW: {total_vdw:.6f} kJ/mol")


def test_three_atom_system():
    """Test with 3 atoms to see interaction counting"""
    
    print("\n\nThree atom system test:")
    
    state = MCState()
    state.info.box = [10.0, 10.0, 10.0]
    state.info.cutoff = 5.0
    
    ff = MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.ljEps = [1.0]
    ff.ljSigma = [0.35]
    state.forcefield = ff
    
    # Three atoms in a line, 0.5 nm apart
    atoms = []
    residues = []
    
    positions = [
        [5.0, 5.0, 5.0],
        [5.5, 5.0, 5.0],
        [6.0, 5.0, 5.0]
    ]
    
    for i, pos in enumerate(positions):
        atom = MCAtom()
        atom.x, atom.y, atom.z = pos
        atom.charge = 0.0
        atom.type = 0
        atoms.append(atom)
        
        res = MCResidue()
        res.active = True
        res.fixed = True
        res.atomStart = i
        res.atomCount = 1
        res.type = 0
        residues.append(res)
    
    state.atoms = atoms
    state.activeAtomCount = 3
    state.residues = residues
    state.activeResidueCount = 3
    
    # Calculate energy
    computeSystemEnergyCutoff(state)
    
    # Get residue energies
    energies = [res.energy_vdw for res in state.residues]
    total = sum(energies)
    
    print(f"  Residue energies: {[f'{e:.6f}' for e in energies]}")
    print(f"  Total: {total:.6f} kJ/mol")
    
    # Calculate expected interactions:
    # 1-2: r=0.5
    # 1-3: r=1.0  
    # 2-3: r=0.5
    r_vals = [0.5, 1.0, 0.5]
    expected_total = 0.0
    for r in r_vals:
        sr6 = (0.35/r)**6
        expected_total += 4.0 * 1.0 * (sr6*sr6 - sr6)
    
    print(f"  Expected (single count): {expected_total:.6f} kJ/mol")
    print(f"  Expected (double count): {2*expected_total:.6f} kJ/mol")
    
    ratio = total / expected_total
    print(f"  Ratio (actual/expected): {ratio:.3f}")


def test_residue_assignment():
    """Test how LJ energy is assigned to residues"""
    
    print("\n\nResidue assignment test:")
    
    state = MCState()
    state.info.box = [10.0, 10.0, 10.0]
    state.info.cutoff = 5.0
    
    ff = MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.ljEps = [1.0]
    ff.ljSigma = [0.35]
    state.forcefield = ff
    
    # Create 2 residues with 2 atoms each
    atoms = []
    residues = []
    
    # Residue 1: atoms at (5,5,5) and (5.4,5,5)
    # Residue 2: atoms at (6,5,5) and (6.4,5,5)
    positions = [
        [5.0, 5.0, 5.0],    # Res 1, atom 1
        [5.4, 5.0, 5.0],    # Res 1, atom 2
        [6.0, 5.0, 5.0],    # Res 2, atom 1
        [6.4, 5.0, 5.0],    # Res 2, atom 2
    ]
    
    for pos in positions:
        atom = MCAtom()
        atom.x, atom.y, atom.z = pos
        atom.charge = 0.0
        atom.type = 0
        atoms.append(atom)
    
    # Two residues with 2 atoms each
    res1 = MCResidue()
    res1.active = True
    res1.fixed = True
    res1.atomStart = 0
    res1.atomCount = 2
    res1.type = 0
    residues.append(res1)
    
    res2 = MCResidue()
    res2.active = True
    res2.fixed = True
    res2.atomStart = 2
    res2.atomCount = 2
    res2.type = 0
    residues.append(res2)
    
    state.atoms = atoms
    state.activeAtomCount = 4
    state.residues = residues
    state.activeResidueCount = 2
    
    # Calculate energy
    computeSystemEnergyCutoff(state)
    
    print(f"  Residue 1 (2 atoms) energy: {state.residues[0].energy_vdw:.6f} kJ/mol")
    print(f"  Residue 2 (2 atoms) energy: {state.residues[1].energy_vdw:.6f} kJ/mol")
    print(f"  Total: {state.residues[0].energy_vdw + state.residues[1].energy_vdw:.6f} kJ/mol")
    
    # Calculate individual interactions
    print("\n  Individual distances:")
    for i in range(4):
        for j in range(i+1, 4):
            dx = atoms[i].x - atoms[j].x
            r = abs(dx)
            print(f"    Atom {i+1}-{j+1}: {r:.2f} nm")


if __name__ == "__main__":
    test_single_pair_lj()
    test_three_atom_system()
    test_residue_assignment()
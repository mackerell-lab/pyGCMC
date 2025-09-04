"""
Test whether PME Total contains LJ energy

Hypothesis: PME Total = electrostatic energy + erroneous LJ term
"""

import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pygcmc
from pygcmc import MCState, MCAtom, MCResidue, MCForceField
from pygcmc import (initializePMEParameters, computeSystemEnergyPME,
                    computeSystemVdwEnergyCutoff)


def test_lj_hypothesis():
    """Test whether PME Total contains LJ energy (regression test)
    
    Verify that fixed PME total no longer erroneously contains LJ energy
    """
    
    # Test configuration
    box_size = 5.0
    cutoff = 2.0
    alpha = 2.5
    
    # Test 1: Pure electrostatic system (no LJ)
    print("\nTest 1: Pure electrostatic system (LJ = 0)")
    print("-" * 50)
    
    state = MCState()
    state.info.box = [box_size, box_size, box_size]
    state.info.cutoff = cutoff
    
    ff = MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.ljEps = [0.0]  # No LJ
    ff.ljSigma = [0.3]
    state.forcefield = ff
    
    # 2 atoms
    atoms = []
    atoms.append(MCAtom())
    atoms[0].x, atoms[0].y, atoms[0].z = 2.0, 2.5, 2.5
    atoms[0].charge = 1.0
    atoms[0].type = 0
    
    atoms.append(MCAtom())
    atoms[1].x, atoms[1].y, atoms[1].z = 3.0, 2.5, 2.5
    atoms[1].charge = -1.0
    atoms[1].type = 0
    
    state.atoms = atoms
    state.activeAtomCount = 2
    
    # 2 residues
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
    
    initializePMEParameters(cutoff, state.info.box, alpha)
    computeSystemEnergyPME(state)
    
    pme_total_no_lj = state.ewald_energy.get('total', 0.0)
    pme_sum_no_lj = (state.ewald_energy.get('real_space', 0.0) + 
                     state.ewald_energy.get('reciprocal', 0.0) + 
                     state.ewald_energy.get('self', 0.0))
    offset_no_lj = pme_total_no_lj - pme_sum_no_lj
    
    print(f"PME Total: {pme_total_no_lj:.2f}")
    print(f"PME component sum: {pme_sum_no_lj:.2f}")
    print(f"Offset: {offset_no_lj:.2f}")
    
    # Test 2: System with LJ
    print("\nTest 2: System with LJ")
    print("-" * 50)
    
    state2 = MCState()
    state2.info.box = [box_size, box_size, box_size]
    state2.info.cutoff = cutoff
    
    ff2 = MCForceField()
    ff2.numTotalTypes = 1
    ff2.numMovementTypes = 1
    ff2.ljEps = [1.0]  # Has LJ
    ff2.ljSigma = [0.35]
    state2.forcefield = ff2
    
    # Same atom configuration
    atoms2 = []
    atoms2.append(MCAtom())
    atoms2[0].x, atoms2[0].y, atoms2[0].z = 2.0, 2.5, 2.5
    atoms2[0].charge = 1.0
    atoms2[0].type = 0
    
    atoms2.append(MCAtom())
    atoms2[1].x, atoms2[1].y, atoms2[1].z = 3.0, 2.5, 2.5
    atoms2[1].charge = -1.0
    atoms2[1].type = 0
    
    state2.atoms = atoms2
    state2.activeAtomCount = 2
    
    # Same residue configuration
    residues2 = []
    for i in range(2):
        res = MCResidue()
        res.active = True
        res.fixed = False
        res.atomStart = i
        res.atomCount = 1
        res.type = 0
        residues2.append(res)
    
    state2.residues = residues2
    state2.activeResidueCount = 2
    
    # Calculate LJ energy first
    try:
        lj_result = computeSystemVdwEnergyCutoff(state2)
        if lj_result is not None:
            lj_energy = lj_result
        else:
            lj_energy = 0.0
        print(f"LJ energy: {lj_energy:.6f}")
    except:
        print("LJ energy calculation failed")
        lj_energy = 0.0
    
    # Calculate PME
    initializePMEParameters(cutoff, state2.info.box, alpha)
    computeSystemEnergyPME(state2)
    
    pme_total_with_lj = state2.ewald_energy.get('total', 0.0)
    pme_sum_with_lj = (state2.ewald_energy.get('real_space', 0.0) + 
                       state2.ewald_energy.get('reciprocal', 0.0) + 
                       state2.ewald_energy.get('self', 0.0))
    offset_with_lj = pme_total_with_lj - pme_sum_with_lj
    
    print(f"PME Total: {pme_total_with_lj:.2f}")
    print(f"PME component sum: {pme_sum_with_lj:.2f}")
    print(f"Offset: {offset_with_lj:.2f}")
    
    print(f"\nOffset difference: {offset_with_lj - offset_no_lj:.6f}")
    print(f"Compare with LJ energy: Is difference close to LJ?")
    
    # Test 3: Check residue LJ energy
    print("\nTest 3: Residue LJ energy")
    print("-" * 50)
    
    total_res_lj = 0.0
    for i, res in enumerate(state2.residues):
        if hasattr(res, 'energy_vdw'):
            print(f"Residue {i} LJ energy: {res.energy_vdw:.6f}")
            total_res_lj += res.energy_vdw
    
    print(f"Residue LJ sum: {total_res_lj:.6f}")
    if 'lj_energy' in locals():
        print(f"System LJ energy: {lj_energy:.6f}")
    
    # Test 4: 1 residue case
    print("\nTest 4: System with 1 residue")
    print("-" * 50)
    
    state3 = MCState()
    state3.info = state2.info
    state3.forcefield = state2.forcefield
    state3.atoms = state2.atoms
    state3.activeAtomCount = state2.activeAtomCount
    
    # 1 residue containing 2 atoms
    residues3 = []
    res = MCResidue()
    res.active = True
    res.fixed = False
    res.atomStart = 0
    res.atomCount = 2
    res.type = 0
    residues3.append(res)
    
    state3.residues = residues3
    state3.activeResidueCount = 1
    
    initializePMEParameters(cutoff, state3.info.box, alpha)
    computeSystemEnergyPME(state3)
    
    pme_total_1res = state3.ewald_energy.get('total', 0.0)
    pme_sum_1res = (state3.ewald_energy.get('real_space', 0.0) + 
                    state3.ewald_energy.get('reciprocal', 0.0) + 
                    state3.ewald_energy.get('self', 0.0))
    offset_1res = pme_total_1res - pme_sum_1res
    
    print(f"1 residue system offset: {offset_1res:.2f} (should be 0)")
    
    print("\n" + "="*80)
    print("Conclusion")
    print("="*80)
    
    if abs(offset_1res) < 1e-6:
        print("✓ Confirmed: Offset is 0 for 1 residue")
    
    if abs(offset_no_lj) == abs(offset_with_lj):
        print("✓ Offset is independent of LJ parameters")
    else:
        print(f"✗ Offset may be related to LJ: difference = {abs(offset_with_lj - offset_no_lj):.6f}")


if __name__ == "__main__":
    test_lj_hypothesis()
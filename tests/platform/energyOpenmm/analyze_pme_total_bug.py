"""
In-depth analysis of PME Total field bug

Find the pattern and source of the offset
"""

import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pytest
import pygcmc
from pygcmc import MCState, MCAtom, MCResidue, MCForceField

def test_pme_total_configurations():
    """Test PME Total correctness under different configurations (regression test)
    
    Verify that the fixed PME total is calculated correctly under various atom/residue configurations
    Create new MCState instances for each test case to avoid state pollution
    """
    
    # Fixed parameters
    box_size = 5.0
    cutoff = 2.0
    alpha = 2.5
    mesh_size = [32, 32, 32]
    spline_order = 4
    
    # ------------------------------------------------------------
    # Initialize PME only once at the beginning of function (set α / mesh / spline order)
    # For different MCState instances, directly call computeSystemEnergyPME to reuse the same set of
    # global parameters, avoiding memory issues from repeatedly freeing/rebuilding grids
    # ------------------------------------------------------------
    try:
        pygcmc.setPMEParameters(alpha, mesh_size, spline_order)
        pygcmc.initializePMEParameters(
            cutoff,
            [box_size, box_size, box_size],
            alpha,
            mesh_size,
            spline_order
        )
    except Exception as e:
        # If initialization fails, it might be because it's already initialized
        print(f"PME initialization warning: {e}")
        pass
    
    # Test case 1: 2 atoms (electroneutral pair)
    print("Test case 1: 2 atoms (electroneutral pair)")
    
    # Create new state
    state1 = MCState()
    state1.info.box = [box_size, box_size, box_size]
    state1.info.cutoff = cutoff
    
    ff1 = MCForceField()
    ff1.numTotalTypes = 1
    ff1.numMovementTypes = 1
    ff1.ljEps = [0.0]
    ff1.ljSigma = [0.3]
    state1.forcefield = ff1
    
    # Create 2 atoms
    atom1 = MCAtom()
    atom1.x, atom1.y, atom1.z = 2.0, 2.0, 2.5
    atom1.charge = 1.0
    atom1.type = 0
    
    atom2 = MCAtom()
    atom2.x, atom2.y, atom2.z = 3.0, 2.0, 2.5
    atom2.charge = -1.0
    atom2.type = 0
    
    state1.atoms = [atom1, atom2]
    state1.activeAtomCount = 2
    
    # 2 residues
    res1_1 = MCResidue()
    res1_1.active = True
    res1_1.fixed = False
    res1_1.atomStart = 0
    res1_1.atomCount = 1
    res1_1.type = 0
    
    res1_2 = MCResidue()
    res1_2.active = True
    res1_2.fixed = False
    res1_2.atomStart = 1
    res1_2.atomCount = 1
    res1_2.type = 0
    
    state1.residues = [res1_1, res1_2]
    state1.activeResidueCount = 2
    
    # Calculate energy directly (PME already initialized at function start)
    pygcmc.computeSystemEnergyPME(state1)
    
    pme_total1 = state1.ewald_energy.get('total', 0.0)
    pme_sum1 = (state1.ewald_energy.get('real_space', 0.0) + 
                state1.ewald_energy.get('reciprocal', 0.0) + 
                state1.ewald_energy.get('self', 0.0))
    offset1 = pme_total1 - pme_sum1
    
    print(f"  Total: {pme_total1:.6f}, Sum: {pme_sum1:.6f}, Offset: {offset1:.6f}")
    assert abs(offset1) < 1e-6, f"2-atom system offset should be 0, actual: {offset1}"
    
    # Test case 2: 3 atoms
    print("\nTest case 2: 3 atoms")
    
    # Create new state
    state2 = MCState()
    state2.info.box = [box_size, box_size, box_size]
    state2.info.cutoff = cutoff
    
    ff2 = MCForceField()
    ff2.numTotalTypes = 1
    ff2.numMovementTypes = 1
    ff2.ljEps = [0.0]
    ff2.ljSigma = [0.3]
    state2.forcefield = ff2
    
    # Create 3 atoms
    positions2 = [[2.0, 2.0, 2.5], [3.0, 2.0, 2.5], [2.5, 3.0, 2.5]]
    charges2 = [1.0, -1.0, 0.5]
    
    atoms2 = []
    for i in range(3):
        atom = MCAtom()
        atom.x, atom.y, atom.z = positions2[i]
        atom.charge = charges2[i]
        atom.type = 0
        atoms2.append(atom)
    
    state2.atoms = atoms2
    state2.activeAtomCount = 3
    
    # 3 residues
    residues2 = []
    for i in range(3):
        res = MCResidue()
        res.active = True
        res.fixed = False
        res.atomStart = i
        res.atomCount = 1
        res.type = 0
        residues2.append(res)
    
    state2.residues = residues2
    state2.activeResidueCount = 3
    
    # Calculate energy directly, no need to reinitialize
    pygcmc.computeSystemEnergyPME(state2)
    
    pme_total2 = state2.ewald_energy.get('total', 0.0)
    pme_sum2 = (state2.ewald_energy.get('real_space', 0.0) + 
                state2.ewald_energy.get('reciprocal', 0.0) + 
                state2.ewald_energy.get('self', 0.0))
    offset2 = pme_total2 - pme_sum2
    
    print(f"  Total: {pme_total2:.6f}, Sum: {pme_sum2:.6f}, Offset: {offset2:.6f}")
    assert abs(offset2) < 1e-6, f"3-atom system offset should be 0, actual: {offset2}"
    
    # Test case 3: 4-atom system (1 residue with 4 atoms)
    print("\nTest case 3: 4-atom system (1 residue with 4 atoms)")
    
    # Create new state
    state3 = MCState()
    state3.info.box = [box_size, box_size, box_size]
    state3.info.cutoff = cutoff
    
    ff3 = MCForceField()
    ff3.numTotalTypes = 1
    ff3.numMovementTypes = 1
    ff3.ljEps = [0.0]
    ff3.ljSigma = [0.3]
    state3.forcefield = ff3
    
    # Create 4 atoms
    positions3 = [[2.0, 2.0, 2.5], [3.0, 2.0, 2.5], [3.0, 3.0, 2.5], [2.0, 3.0, 2.5]]
    charges3 = [1.0, -1.0, 1.0, -1.0]
    
    atoms3 = []
    for i in range(4):
        atom = MCAtom()
        atom.x, atom.y, atom.z = positions3[i]
        atom.charge = charges3[i]
        atom.type = 0
        atoms3.append(atom)
    
    state3.atoms = atoms3
    state3.activeAtomCount = 4
    
    # 1 residue with 4 atoms
    res3 = MCResidue()
    res3.active = True
    res3.fixed = False
    res3.atomStart = 0
    res3.atomCount = 4
    res3.type = 0
    
    state3.residues = [res3]
    state3.activeResidueCount = 1
    
    # Calculate energy directly, no need to reinitialize
    pygcmc.computeSystemEnergyPME(state3)
    
    pme_total3 = state3.ewald_energy.get('total', 0.0)
    pme_sum3 = (state3.ewald_energy.get('real_space', 0.0) + 
                state3.ewald_energy.get('reciprocal', 0.0) + 
                state3.ewald_energy.get('self', 0.0))
    offset3 = pme_total3 - pme_sum3
    
    print(f"  Total: {pme_total3:.6f}, Sum: {pme_sum3:.6f}, Offset: {offset3:.6f}")
    assert abs(offset3) < 1e-6, f"1 residue with 4 atoms offset should be 0, actual: {offset3}"
    
    # Verify all offsets are 0
    print(f"\nAll tests passed! Offsets: {offset1:.6f}, {offset2:.6f}, {offset3:.6f}")
    assert abs(offset1) < 1e-6 and abs(offset2) < 1e-6 and abs(offset3) < 1e-6, \
        "PME Total offset should be 0 under all configurations"


if __name__ == "__main__":
    test_pme_total_configurations()
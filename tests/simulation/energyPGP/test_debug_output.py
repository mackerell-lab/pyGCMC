# tests/simulation/energyPGP/test_debug_output.py
"""
Test with debug output to see what's happening in PME.
"""

import sys
sys.path.insert(0, '/home/zhaomt/gcmc/test107/pygcmc_dev/build')
import pygcmc
from pygcmc import MCState, MCAtom, MCResidue, MCForceField
from pygcmc import setPMEParameters, initializePMEParameters, computeSystemEnergyPME
import math


def test_debug_output():
    """Test to see debug output."""
    print("\n=== Running with Debug Output ===")
    
    # Create minimal system
    state = MCState()
    state.info.box = [10.0, 10.0, 10.0]
    state.info.cutoff = 1.2
    
    ff = MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.ljEps = [0.0]
    ff.ljSigma = [0.3]
    state.forcefield = ff
    
    atoms = []
    atom = MCAtom()
    atom.x = 5.0
    atom.y = 5.0
    atom.z = 5.0
    atom.charge = 1.0
    atom.type = 0
    atoms.append(atom)
    
    atom2 = MCAtom()
    atom2.x = 6.0
    atom2.y = 5.0
    atom2.z = 5.0
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
    
    # Initialize with explicit alpha
    alpha = 2.5
    mesh_size = [32, 32, 32]
    
    print("\nCalling setPMEParameters with alpha=2.5...")
    setPMEParameters(alpha, mesh_size, 4, 1e-6)
    
    print("\nCalling initializePMEParameters...")
    initializePMEParameters(state.info.cutoff, state.info.box, alpha)
    
    print("\nCalling computeSystemEnergyPME...")
    computeSystemEnergyPME(state)
    
    print(f"\nFinal energy: {state.ewald_energy.get('total', 0.0):.2f} kJ/mol")


if __name__ == "__main__":
    test_debug_output()
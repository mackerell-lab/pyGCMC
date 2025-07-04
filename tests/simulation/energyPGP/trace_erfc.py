# tests/simulation/energyPGP/trace_erfc.py
"""
Trace erfcApprox calls.
"""

import sys
sys.path.insert(0, '/home/zhaomt/gcmc/test107/pygcmc_dev/build')
import pygcmc
from pygcmc import MCState, MCAtom, MCResidue, MCForceField
from pygcmc import setPMEParameters, initializePMEParameters, computeSystemEnergyPME
import math


def trace_erfc():
    """Create a minimal test to trace erfcApprox."""
    print("\n=== Tracing erfcApprox ===")
    
    # Create state
    state = MCState()
    state.info.box = [10.0, 10.0, 10.0]
    state.info.cutoff = 1.2
    
    ff = MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.ljEps = [0.0]
    ff.ljSigma = [0.3]
    state.forcefield = ff
    
    # Just two atoms
    atoms = []
    atom1 = MCAtom()
    atom1.x = 5.0
    atom1.y = 5.0
    atom1.z = 5.0
    atom1.charge = 1.0
    atom1.type = 0
    atoms.append(atom1)
    
    atom2 = MCAtom()
    atom2.x = 6.0  # 1 nm apart
    atom2.y = 5.0
    atom2.z = 5.0
    atom2.charge = -1.0
    atom2.type = 0
    atoms.append(atom2)
    
    state.atoms = atoms
    state.activeAtomCount = 2
    
    # Just two residues
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
    
    # Initialize with debugging
    print("\nInitializing PME parameters...")
    setPMEParameters(2.5, [32, 32, 32], 4, 1e-6)
    initializePMEParameters(state.info.cutoff, state.info.box, 2.5)
    
    print("\nComputing energy (should trigger debug output)...")
    computeSystemEnergyPME(state)
    
    print(f"\nResidue 0 energy: {state.residues[0].energy_elec:.2f} kJ/mol")
    print(f"Residue 1 energy: {state.residues[1].energy_elec:.2f} kJ/mol")
    print(f"Total system energy: {state.systemEnergy:.2f} kJ/mol")


if __name__ == "__main__":
    trace_erfc()
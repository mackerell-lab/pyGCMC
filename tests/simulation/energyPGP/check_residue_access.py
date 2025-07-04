# tests/simulation/energyPGP/check_residue_access.py
"""
Check residue access pattern.
"""

import sys
sys.path.insert(0, '/home/zhaomt/gcmc/test107/pygcmc_dev/build')
import pygcmc
from pygcmc import MCState, MCAtom, MCResidue, MCForceField
from pygcmc import setPMEParameters, initializePMEParameters, computeSystemEnergyPGP


def check_residue_access():
    """Check how residues are accessed."""
    print("\n=== Checking Residue Access ===")
    
    # Create state
    state = MCState()
    state.info.box = [10.0, 10.0, 10.0]
    state.info.cutoff = 2.0
    
    ff = MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.ljEps = [1.0]
    ff.ljSigma = [0.34]
    state.forcefield = ff
    
    # Two atoms
    atoms = []
    atom1 = MCAtom()
    atom1.x = 5.0
    atom1.y = 5.0
    atom1.z = 5.0
    atom1.charge = 0.0
    atom1.type = 0
    atoms.append(atom1)
    
    atom2 = MCAtom()
    atom2.x = 5.068  # Very close
    atom2.y = 5.0
    atom2.z = 5.0
    atom2.charge = 0.0
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
    
    # Initialize
    setPMEParameters(2.5, [32, 32, 32], 4, 1e-6)
    initializePMEParameters(state.info.cutoff, state.info.box, 2.5)
    
    # Compute
    computeSystemEnergyPGP(state)
    
    # Access residues different ways
    print("\nDirect access:")
    print(f"  state.residues[0].energy_vdw = {state.residues[0].energy_vdw}")
    print(f"  state.residues[1].energy_vdw = {state.residues[1].energy_vdw}")
    
    print("\nList comprehension:")
    energies = [res.energy_vdw for res in state.residues if res.active]
    print(f"  energies = {energies}")
    
    print("\nSum:")
    lj_sum = sum(res.energy_vdw for res in state.residues if res.active)
    print(f"  sum = {lj_sum}")
    print(f"  sum/2 = {lj_sum/2}")
    
    print("\nManual iteration:")
    total = 0.0
    for i, res in enumerate(state.residues):
        if res.active:
            print(f"  residue {i}: energy_vdw = {res.energy_vdw}")
            total += res.energy_vdw
    print(f"  manual total = {total}")


if __name__ == "__main__":
    check_residue_access()
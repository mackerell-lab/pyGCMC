# tests/simulation/energyPGP/test_with_logs.py
"""
Test with logging enabled to see debug output.
"""

import sys
sys.path.insert(0, '/home/zhaomt/gcmc/test107/pygcmc_dev/build')
import pygcmc
from pygcmc import MCState, MCAtom, MCResidue, MCForceField
from pygcmc import setPMEParameters, initializePMEParameters, computeSystemEnergyPME
import math


def test_with_logs():
    """Run test to see log output."""
    print("\n=== Test with Logs ===")
    
    # Minimal system
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
    
    # Initialize
    alpha = 2.5
    mesh_size = [32, 32, 32]
    setPMEParameters(alpha, mesh_size, 4, 1e-6)
    initializePMEParameters(state.info.cutoff, state.info.box, alpha)
    
    # This should trigger the debug output
    computeSystemEnergyPME(state)
    
    # Print results
    print(f"\nResults:")
    print(f"Total energy: {state.ewald_energy.get('total', 0.0):.2f}")
    print(f"Residue[0] energy: {state.residues[0].energy_elec:.2f}")
    
    # If erfcApprox returns 2.0, then:
    # pair_energy = -1 * 2.0 / 1.0 = -2.0 (no COULOMB)
    # residue gets full energy: -2.0
    # after COULOMB: -2.0 * 138.935456 = -277.87
    # But we're seeing -19303, which is COULOMB²
    
    print(f"\nExpected if erfcApprox=2.0:")
    print(f"  Residue energy (1x COULOMB): {-2.0 * 138.935456:.2f}")
    print(f"  Residue energy (2x COULOMB): {-2.0 * 138.935456**2:.2f}")
    print(f"\nActual residue energy: {state.residues[0].energy_elec:.2f}")


if __name__ == "__main__":
    test_with_logs()
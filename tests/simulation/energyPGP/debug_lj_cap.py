# tests/simulation/energyPGP/debug_lj_cap.py
"""
Debug LJ energy capping issue.
"""

import sys
sys.path.insert(0, '/home/zhaomt/gcmc/test107/pygcmc_dev/build')
import pygcmc
from pygcmc import MCState, MCAtom, MCResidue, MCForceField
from pygcmc import setPMEParameters, initializePMEParameters, computeSystemEnergyPGP


def debug_lj_cap():
    """Debug why LJ energy is capped."""
    print("\n=== Debugging LJ Energy Cap ===")
    
    # Create a simple system with very short distance
    state = MCState()
    state.info.box = [10.0, 10.0, 10.0]
    state.info.cutoff = 2.0
    
    ff = MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.ljEps = [1.0]  # 1.0 kJ/mol
    ff.ljSigma = [0.34]  # nm
    state.forcefield = ff
    
    # Two atoms very close
    distance = 0.068  # 0.2 * sigma
    atoms = []
    
    atom1 = MCAtom()
    atom1.x = 5.0
    atom1.y = 5.0
    atom1.z = 5.0
    atom1.charge = 0.0
    atom1.type = 0
    atoms.append(atom1)
    
    atom2 = MCAtom()
    atom2.x = 5.0 + distance
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
    
    # Initialize and compute
    alpha = 2.5
    setPMEParameters(alpha, [32, 32, 32], 4, 1e-6)
    initializePMEParameters(state.info.cutoff, state.info.box, alpha)
    
    # Get state before calling C++
    print(f"\nBefore computeSystemEnergyPGP:")
    print(f"  Residue 0 VDW: {state.residues[0].energy_vdw}")
    print(f"  Residue 1 VDW: {state.residues[1].energy_vdw}")
    
    computeSystemEnergyPGP(state)
    
    print(f"\nAfter computeSystemEnergyPGP:")
    print(f"  Residue 0 VDW: {state.residues[0].energy_vdw}")
    print(f"  Residue 1 VDW: {state.residues[1].energy_vdw}")
    
    # Calculate expected
    sigma = 0.34
    epsilon = 1.0
    r_ratio = sigma / distance
    expected = 4 * epsilon * (r_ratio**12 - r_ratio**6)
    
    print(f"\nAnalytical calculation:")
    print(f"  Distance: {distance} nm")
    print(f"  Sigma/r: {r_ratio:.2f}")
    print(f"  (σ/r)^6: {r_ratio**6:.2e}")
    print(f"  (σ/r)^12: {r_ratio**12:.2e}")
    print(f"  Expected LJ: {expected:.2e} kJ/mol")
    
    total_vdw = state.residues[0].energy_vdw + state.residues[1].energy_vdw
    print(f"\nTotal VDW from residues: {total_vdw}")
    print(f"Per pair (÷2): {total_vdw/2}")


if __name__ == "__main__":
    debug_lj_cap()
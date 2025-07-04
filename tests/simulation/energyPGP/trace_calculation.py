# tests/simulation/energyPGP/trace_calculation.py
"""
Trace the calculation step by step.
"""

import sys
sys.path.insert(0, '/home/zhaomt/gcmc/test107/pygcmc_dev/build')
import pygcmc
from pygcmc import MCState, MCAtom, MCResidue, MCForceField
from pygcmc import setPMEParameters, initializePMEParameters, computeSystemEnergyPME
import math


def trace_calculation():
    """Trace PME calculation step by step."""
    print("\n=== Tracing PME Calculation ===")
    
    # Expected values
    r = 1.0  # nm
    q1, q2 = 1.0, -1.0
    alpha = 2.5
    
    print(f"\nExpected calculation:")
    print(f"  Distance: {r} nm")
    print(f"  Charges: q1={q1}, q2={q2}")
    print(f"  Alpha: {alpha}")
    
    expected_erfc = math.erfc(alpha * r)
    print(f"  erfc({alpha}*{r}) = {expected_erfc:.6f}")
    
    expected_pair = q1 * q2 * expected_erfc / r
    print(f"  Pair energy = {q1}*{q2}*{expected_erfc:.6f}/{r} = {expected_pair:.6f}")
    
    expected_with_coulomb = expected_pair * 138.935456
    print(f"  With COULOMB: {expected_with_coulomb:.6f} kJ/mol")
    
    expected_per_residue = expected_with_coulomb / 2.0
    print(f"  Per residue (half): {expected_per_residue:.6f} kJ/mol")
    
    # Now actual calculation
    print("\n\nActual calculation:")
    
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
    
    # Two atoms at 1 nm distance
    atoms = []
    atom1 = MCAtom()
    atom1.x = 5.0
    atom1.y = 5.0
    atom1.z = 5.0
    atom1.charge = q1
    atom1.type = 0
    atoms.append(atom1)
    
    atom2 = MCAtom()
    atom2.x = 5.0 + r
    atom2.y = 5.0
    atom2.z = 5.0
    atom2.charge = q2
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
    
    # Initialize PME
    setPMEParameters(alpha, [32, 32, 32], 4, 1e-6)
    initializePMEParameters(state.info.cutoff, state.info.box, alpha)
    
    # Compute energy
    computeSystemEnergyPME(state)
    
    # Get results
    res_energy = state.residues[0].energy_elec
    print(f"  Residue 0 energy: {res_energy:.6f} kJ/mol")
    print(f"  Residue 1 energy: {state.residues[1].energy_elec:.6f} kJ/mol")
    
    # If erfcApprox returns 1.0:
    print("\n\nIf erfcApprox returns 1.0:")
    wrong_pair = q1 * q2 * 1.0 / r
    print(f"  Pair energy = {q1}*{q2}*1.0/{r} = {wrong_pair:.6f}")
    print(f"  With COULOMB: {wrong_pair * 138.935456:.6f} kJ/mol")
    print(f"  Per residue: {wrong_pair * 138.935456 / 2.0:.6f} kJ/mol")
    
    # Check if COULOMB is applied twice
    print("\n\nIf COULOMB is applied twice:")
    double_coulomb = wrong_pair * 138.935456 * 138.935456 / 2.0
    print(f"  Per residue: {double_coulomb:.6f} kJ/mol")
    
    # This matches!
    if abs(res_energy - double_coulomb) < 1.0:
        print("\n✗ COULOMB is being applied twice AND erfcApprox returns 1.0!")


if __name__ == "__main__":
    trace_calculation()
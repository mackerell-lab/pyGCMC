# tests/simulation/energyPGP/debug_pme_detail.py
"""
Detailed debug of PME energy calculation.
"""

import pygcmc
from pygcmc import MCState, MCAtom, MCResidue, MCForceField
from pygcmc import setPMEParameters, initializePMEParameters, computeSystemEnergyPME
import math


def debug_pme_detail():
    """Debug PME energy calculation in detail."""
    print("\n=== Detailed PME Debug ===")
    
    # Create minimal system
    state = MCState()
    state.info.box = [10.0, 10.0, 10.0]
    state.info.cutoff = 1.2
    
    # Force field (no LJ)
    ff = MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.ljEps = [0.0]
    ff.ljSigma = [0.3]
    state.forcefield = ff
    
    # Two charged particles at exactly 1.0 nm
    atoms = []
    
    atom1 = MCAtom()
    atom1.x = 5.0
    atom1.y = 5.0
    atom1.z = 5.0
    atom1.charge = 1.0
    atom1.type = 0
    atoms.append(atom1)
    
    atom2 = MCAtom()
    atom2.x = 6.0  # Exactly 1.0 nm away
    atom2.y = 5.0
    atom2.z = 5.0
    atom2.charge = -1.0
    atom2.type = 0
    atoms.append(atom2)
    
    state.atoms = atoms
    state.activeAtomCount = 2
    
    # Create residues with zero initial energies
    residues = []
    for i in range(2):
        res = MCResidue()
        res.active = True
        res.fixed = False
        res.atomStart = i
        res.atomCount = 1
        res.type = 0
        res.energy_vdw = 0.0
        res.energy_elec = 0.0
        residues.append(res)
    
    state.residues = residues
    state.activeResidueCount = 2
    
    # Initialize PME
    alpha = 2.5
    mesh_size = [32, 32, 32]
    setPMEParameters(alpha, mesh_size, 4, 1e-6)
    initializePMEParameters(state.info.cutoff, state.info.box, alpha)
    
    print(f"System setup:")
    print(f"  Distance: 1.0 nm")
    print(f"  Charges: +1.0, -1.0")
    print(f"  Alpha: {alpha}")
    print(f"  Cutoff: {state.info.cutoff} nm")
    
    # Calculate expected values
    kC = 138.935456  # Coulomb constant
    r = 1.0  # nm
    erfc_val = math.erfc(alpha * r)
    
    print(f"\nExpected calculations:")
    print(f"  erfc(αr) = erfc({alpha} * {r}) = {erfc_val:.6f}")
    print(f"  Real-space pair energy (no COULOMB): q1*q2*erfc(αr)/r = {-1.0 * erfc_val / r:.6f}")
    print(f"  Real-space total (no COULOMB): {-1.0 * erfc_val / r:.6f}")
    print(f"  Real-space with COULOMB: {-kC * erfc_val / r:.6f} kJ/mol")
    
    # Compute energy
    computeSystemEnergyPME(state)
    
    print(f"\nActual PME results:")
    print(f"  real_space: {state.ewald_energy.get('real_space', 0.0):.6f}")
    print(f"  reciprocal: {state.ewald_energy.get('reciprocal', 0.0):.6f}")
    print(f"  self: {state.ewald_energy.get('self', 0.0):.6f}")
    print(f"  total: {state.ewald_energy.get('total', 0.0):.6f}")
    
    print(f"\nResidue energies:")
    for i, res in enumerate(state.residues):
        print(f"  Residue {i}: vdw={res.energy_vdw:.6f}, elec={res.energy_elec:.6f}")
    
    # Check if COULOMB was applied multiple times
    res0_elec = state.residues[0].energy_elec
    print(f"\nChecking for multiple COULOMB applications:")
    print(f"  Residue 0 elec energy: {res0_elec:.6f}")
    print(f"  Expected (half pair, no COULOMB): {-0.5 * erfc_val / r:.6f}")
    print(f"  Expected (half pair, 1x COULOMB): {-0.5 * kC * erfc_val / r:.6f}")
    print(f"  Expected (half pair, 2x COULOMB): {-0.5 * kC * kC * erfc_val / r:.6f}")
    
    # Check ratios
    if abs(res0_elec) > 1e-6:
        ratio1 = res0_elec / (-0.5 * erfc_val / r)
        ratio2 = res0_elec / (-0.5 * kC * erfc_val / r)
        print(f"\nRatios:")
        print(f"  Actual / (no COULOMB): {ratio1:.6f}")
        print(f"  Actual / (1x COULOMB): {ratio2:.6f}")
        
        if abs(ratio1 - kC * kC) < 0.1:
            print(f"\n⚠️ COULOMB appears to be applied TWICE!")
            print(f"  {kC}² = {kC * kC:.6f}")


if __name__ == "__main__":
    debug_pme_detail()
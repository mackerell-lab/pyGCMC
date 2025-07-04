# tests/simulation/energyPGP/fix_pme_test.py
"""
Test if removing COULOMB multiplication from residue energies fixes the issue.
"""

import pygcmc
from pygcmc import MCState, MCAtom, MCResidue, MCForceField
from pygcmc import setPMEParameters, initializePMEParameters, computeSystemEnergyPME
import math


def test_pme_fix():
    """Test PME with manual fix."""
    print("\n=== Testing PME Fix ===")
    
    # Create two-particle system
    state = MCState()
    state.info.box = [10.0, 10.0, 10.0]
    state.info.cutoff = 1.2
    
    ff = MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.ljEps = [0.0]
    ff.ljSigma = [0.3]
    state.forcefield = ff
    
    # Two particles
    atoms = []
    
    atom1 = MCAtom()
    atom1.x = 5.0
    atom1.y = 5.0
    atom1.z = 5.0
    atom1.charge = 1.0
    atom1.type = 0
    atoms.append(atom1)
    
    atom2 = MCAtom()
    atom2.x = 6.0  # 1 nm away
    atom2.y = 5.0
    atom2.z = 5.0
    atom2.charge = -1.0
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
    alpha = 2.5
    mesh_size = [32, 32, 32]
    setPMEParameters(alpha, mesh_size, 4, 1e-6)
    initializePMEParameters(state.info.cutoff, state.info.box, alpha)
    
    # Compute
    computeSystemEnergyPME(state)
    
    print("PME Results:")
    print(f"  Real-space: {state.ewald_energy.get('real_space', 0.0):.6f}")
    print(f"  Reciprocal: {state.ewald_energy.get('reciprocal', 0.0):.6f}")
    print(f"  Self: {state.ewald_energy.get('self', 0.0):.6f}")
    print(f"  Total: {state.ewald_energy.get('total', 0.0):.6f}")
    
    print(f"\nResidue energies:")
    print(f"  Residue[0] elec: {state.residues[0].energy_elec:.6f}")
    print(f"  Residue[1] elec: {state.residues[1].energy_elec:.6f}")
    
    # Manual fix - divide by COULOMB to get back original value
    kC = 138.935456
    res0_fixed = state.residues[0].energy_elec / kC
    res1_fixed = state.residues[1].energy_elec / kC
    
    print(f"\nDividing residue energies by COULOMB:")
    print(f"  Residue[0] / COULOMB: {res0_fixed:.6f}")
    print(f"  Residue[1] / COULOMB: {res1_fixed:.6f}")
    
    # Calculate expected
    r = 1.0
    erfc_val = math.erfc(alpha * r)
    expected_per_res_no_coulomb = -0.5 * erfc_val / r
    expected_per_res_with_coulomb = expected_per_res_no_coulomb * kC
    
    print(f"\nExpected values:")
    print(f"  Per residue (no COULOMB): {expected_per_res_no_coulomb:.6f}")
    print(f"  Per residue (with COULOMB): {expected_per_res_with_coulomb:.6f}")
    
    # Check if dividing by COULOMB gives the right value
    print(f"\nAnalysis:")
    print(f"  res0_fixed / expected_no_coulomb = {res0_fixed / expected_per_res_no_coulomb:.1f}")
    
    # Try dividing by COULOMB² 
    res0_fixed2 = state.residues[0].energy_elec / (kC * kC)
    print(f"  Residue[0] / COULOMB²: {res0_fixed2:.6f}")
    print(f"  res0_fixed2 / expected_no_coulomb = {res0_fixed2 / expected_per_res_no_coulomb:.6f}")
    
    # Calculate corrected total
    corrected_residue_total = 2 * res0_fixed  # Two residues with same energy
    corrected_total = corrected_residue_total + state.ewald_energy.get('reciprocal', 0.0) + state.ewald_energy.get('self', 0.0)
    
    print(f"\nCorrected total energy: {corrected_total:.6f}")
    
    # Expected total from direct Coulomb
    expected_total = -kC / r  # Direct Coulomb for opposite charges
    print(f"Expected total (direct Coulomb): {expected_total:.6f}")
    
    # The issue appears to be that residue.energy_elec is being multiplied by COULOMB²
    # This suggests COULOMB is applied once in PMEReal.cpp and once in PMEComposite.cpp
    
    print(f"\n⚠️ DIAGNOSIS: Residue energies have COULOMB² applied instead of COULOMB")
    print("   This means COULOMB is being applied twice somewhere in the code")


if __name__ == "__main__":
    test_pme_fix()
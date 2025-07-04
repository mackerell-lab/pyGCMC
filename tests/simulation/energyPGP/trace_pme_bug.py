# tests/simulation/energyPGP/trace_pme_bug.py
"""
Trace PME bug by checking intermediate values.
"""

import pygcmc
from pygcmc import MCState, MCAtom, MCResidue, MCForceField
from pygcmc import setPMEParameters, initializePMEParameters, computeSystemEnergyPME
import math


def trace_pme_bug():
    """Trace PME calculation step by step."""
    print("\n=== Tracing PME Bug ===")
    
    # Create minimal system - single particle first
    state = MCState()
    state.info.box = [10.0, 10.0, 10.0]
    state.info.cutoff = 1.2
    
    ff = MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.ljEps = [0.0]
    ff.ljSigma = [0.3]
    state.forcefield = ff
    
    # Single charged particle
    atoms = []
    atom = MCAtom()
    atom.x = 5.0
    atom.y = 5.0
    atom.z = 5.0
    atom.charge = 1.0
    atom.type = 0
    atoms.append(atom)
    state.atoms = atoms
    state.activeAtomCount = 1
    
    # Single residue
    residues = []
    res = MCResidue()
    res.active = True
    res.fixed = False
    res.atomStart = 0
    res.atomCount = 1
    res.type = 0
    residues.append(res)
    state.residues = residues
    state.activeResidueCount = 1
    
    # Initialize PME
    alpha = 2.5
    mesh_size = [32, 32, 32]
    setPMEParameters(alpha, mesh_size, 4, 1e-6)
    initializePMEParameters(state.info.cutoff, state.info.box, alpha)
    
    print("=== Test 1: Single particle ===")
    print(f"Charge: {atom.charge}")
    
    # Compute energy
    computeSystemEnergyPME(state)
    
    print(f"\nResults:")
    print(f"  Real-space: {state.ewald_energy.get('real_space', 0.0)}")
    print(f"  Reciprocal: {state.ewald_energy.get('reciprocal', 0.0)}")
    print(f"  Self: {state.ewald_energy.get('self', 0.0)}")
    print(f"  Total: {state.ewald_energy.get('total', 0.0)}")
    print(f"  Residue elec: {state.residues[0].energy_elec}")
    
    # For single particle, real-space should be 0 (no pairs)
    # Self energy should be -alpha/sqrt(pi) * q^2 * COULOMB
    kC = 138.935456
    expected_self = -alpha / math.sqrt(math.pi) * 1.0 * kC
    print(f"\nExpected self energy: {expected_self}")
    
    # Now test with particle pair
    print("\n\n=== Test 2: Two particles ===")
    
    # Reset state
    state = MCState()
    state.info.box = [10.0, 10.0, 10.0]
    state.info.cutoff = 1.2
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
    
    # Reinitialize PME (just to be sure)
    setPMEParameters(alpha, mesh_size, 4, 1e-6)
    initializePMEParameters(state.info.cutoff, state.info.box, alpha)
    
    print("Two particles at 1 nm distance")
    print("Charges: +1, -1")
    
    # Before computation
    print(f"\nBefore computation:")
    print(f"  Residue[0] elec: {state.residues[0].energy_elec}")
    print(f"  Residue[1] elec: {state.residues[1].energy_elec}")
    
    # Compute
    computeSystemEnergyPME(state)
    
    print(f"\nAfter computation:")
    print(f"  Real-space: {state.ewald_energy.get('real_space', 0.0)}")
    print(f"  Reciprocal: {state.ewald_energy.get('reciprocal', 0.0)}")
    print(f"  Self: {state.ewald_energy.get('self', 0.0)}")
    print(f"  Total: {state.ewald_energy.get('total', 0.0)}")
    print(f"  Residue[0] elec: {state.residues[0].energy_elec}")
    print(f"  Residue[1] elec: {state.residues[1].energy_elec}")
    
    # Calculate what values should be
    r = 1.0
    erfc_val = math.erfc(alpha * r)
    
    print(f"\nExpected values:")
    print(f"  erfc({alpha} * {r}) = {erfc_val}")
    print(f"  Real-space pair (no COULOMB): {-erfc_val / r}")
    print(f"  Real-space pair (with COULOMB): {-kC * erfc_val / r}")
    print(f"  Per residue (no COULOMB): {-0.5 * erfc_val / r}")
    print(f"  Per residue (with COULOMB): {-0.5 * kC * erfc_val / r}")
    
    # The issue appears to be that residue energies are getting a massive value
    # Let's check if it's related to some global state
    
    print(f"\nDiagnostics:")
    actual_res_energy = state.residues[0].energy_elec
    expected_with_coulomb = -0.5 * kC * erfc_val / r
    
    if abs(actual_res_energy) > 1000:
        print(f"⚠️ Residue energy is extremely large: {actual_res_energy}")
        print(f"   This suggests a bug in the PME real-space calculation")
        
        # Check if it's a simple multiplication issue
        factor = actual_res_energy / expected_with_coulomb
        print(f"   Factor: actual / expected = {factor}")
        
        # Check some common issues
        print(f"   Is it COULOMB²? {kC * kC} = {kC**2}")
        print(f"   Is it related to box volume? {10.0 * 10.0 * 10.0} = 1000")


if __name__ == "__main__":
    trace_pme_bug()
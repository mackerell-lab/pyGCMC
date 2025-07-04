# tests/simulation/energyPGP/debug_pme_minimal.py
"""
Minimal PME debug to isolate the issue.
"""

import pygcmc
from pygcmc import MCState, MCAtom, MCResidue, MCForceField
from pygcmc import setPMEParameters, initializePMEParameters
import math


def test_pme_minimal():
    """Minimal test to isolate PME issue."""
    print("\n=== Minimal PME Test ===")
    
    # Create state
    state = MCState()
    state.info.box = [10.0, 10.0, 10.0]
    state.info.cutoff = 1.2
    
    # Force field
    ff = MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.ljEps = [0.0]
    ff.ljSigma = [0.3]
    state.forcefield = ff
    
    # Single atom
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
    res.energy_vdw = 0.0
    res.energy_elec = 0.0
    residues.append(res)
    
    state.residues = residues
    state.activeResidueCount = 1
    
    # Initialize PME
    alpha = 2.5
    mesh_size = [32, 32, 32]
    setPMEParameters(alpha, mesh_size, 4, 1e-6)
    initializePMEParameters(state.info.cutoff, state.info.box, alpha)
    
    print(f"Single particle system:")
    print(f"  Charge: {atom.charge}")
    print(f"  Position: ({atom.x}, {atom.y}, {atom.z})")
    
    # Manually call computeRealSpacePME
    from pygcmc import computeRealSpacePME
    computeRealSpacePME(state, movement_only=False, store_in_residues=True)
    
    print(f"\nAfter computeRealSpacePME:")
    print(f"  state.ewald_energy['real_space']: {state.ewald_energy.get('real_space', 0.0)}")
    print(f"  residue.energy_elec: {state.residues[0].energy_elec}")
    
    # Now apply COULOMB manually
    COULOMB = 138.935456
    state.residues[0].energy_elec *= COULOMB
    
    print(f"\nAfter applying COULOMB:")
    print(f"  residue.energy_elec: {state.residues[0].energy_elec}")
    
    # Test with two particles
    print("\n\n=== Two Particle Test ===")
    
    # Add second atom
    atom2 = MCAtom()
    atom2.x = 6.0  # 1 nm away
    atom2.y = 5.0
    atom2.z = 5.0
    atom2.charge = -1.0
    atom2.type = 0
    atoms.append(atom2)
    
    state.atoms = atoms
    state.activeAtomCount = 2
    
    # Add second residue
    res2 = MCResidue()
    res2.active = True
    res2.fixed = False
    res2.atomStart = 1
    res2.atomCount = 1
    res2.type = 0
    res2.energy_vdw = 0.0
    res2.energy_elec = 0.0
    residues.append(res2)
    
    state.residues = residues
    state.activeResidueCount = 2
    
    # Reset energies
    state.residues[0].energy_elec = 0.0
    state.residues[1].energy_elec = 0.0
    
    # Call computeRealSpacePME again
    computeRealSpacePME(state, movement_only=False, store_in_residues=True)
    
    print(f"Two particles at 1 nm distance:")
    print(f"  Charges: +1, -1")
    print(f"\nAfter computeRealSpacePME:")
    print(f"  state.ewald_energy['real_space']: {state.ewald_energy.get('real_space', 0.0)}")
    print(f"  residue[0].energy_elec: {state.residues[0].energy_elec}")
    print(f"  residue[1].energy_elec: {state.residues[1].energy_elec}")
    
    # Calculate expected
    r = 1.0
    erfc_val = math.erfc(alpha * r)
    expected_pair = -1.0 * erfc_val / r  # q1*q2*erfc(αr)/r
    expected_per_res = expected_pair / 2.0
    
    print(f"\nExpected values:")
    print(f"  erfc(αr) = {erfc_val}")
    print(f"  Pair energy (no COULOMB): {expected_pair}")
    print(f"  Per residue (no COULOMB): {expected_per_res}")
    
    # Check ratio
    if abs(state.residues[0].energy_elec) > 1e-10:
        ratio = state.residues[0].energy_elec / expected_per_res
        print(f"\nRatio actual/expected: {ratio}")


if __name__ == "__main__":
    test_pme_minimal()
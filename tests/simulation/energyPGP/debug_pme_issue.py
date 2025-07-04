# tests/simulation/energyPGP/debug_pme_issue.py
"""
Debug PME total energy issue.
"""

import pygcmc
from pygcmc import MCState, MCAtom, MCResidue, MCForceField
from pygcmc import setPMEParameters, initializePMEParameters, computeSystemEnergyPME


def debug_pme_energy():
    """Debug PME energy calculation."""
    print("\n=== PME Energy Debug ===")
    
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
    
    # Two charged particles
    atoms = []
    
    atom1 = MCAtom()
    atom1.x = 5.0
    atom1.y = 5.0
    atom1.z = 5.0
    atom1.charge = 1.0
    atom1.type = 0
    atoms.append(atom1)
    
    atom2 = MCAtom()
    atom2.x = 6.0  # 1.0 nm away
    atom2.y = 5.0
    atom2.z = 5.0
    atom2.charge = -1.0
    atom2.type = 0
    atoms.append(atom2)
    
    state.atoms = atoms
    state.activeAtomCount = 2
    
    # Create residues
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
    
    # Compute energy
    computeSystemEnergyPME(state)
    
    # Debug output
    print(f"\nEwald energy dictionary:")
    print(f"  real_space: {state.ewald_energy.get('real_space', 0.0):.6f}")
    print(f"  reciprocal: {state.ewald_energy.get('reciprocal', 0.0):.6f}")
    print(f"  self: {state.ewald_energy.get('self', 0.0):.6f}")
    print(f"  total: {state.ewald_energy.get('total', 0.0):.6f}")
    
    print(f"\nResidue energies:")
    residue_total = 0.0
    for i, res in enumerate(state.residues):
        if res.active:
            print(f"  Residue {i}: vdw={res.energy_vdw:.6f}, elec={res.energy_elec:.6f}")
            residue_total += res.energy_vdw + res.energy_elec
    
    print(f"\nResidue total: {residue_total:.6f}")
    
    # Calculate expected total
    expected_total = residue_total + state.ewald_energy.get('reciprocal', 0.0) + state.ewald_energy.get('self', 0.0)
    print(f"Expected total (residue_total + recip + self): {expected_total:.6f}")
    
    # Check if real_space is included in residue energies
    if abs(residue_total - state.ewald_energy.get('real_space', 0.0)) < 1e-6:
        print("\n✓ Residue total equals real_space energy")
    else:
        print(f"\n⚠️ Residue total ({residue_total:.6f}) != real_space ({state.ewald_energy.get('real_space', 0.0):.6f})")


if __name__ == "__main__":
    debug_pme_energy()
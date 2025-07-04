# tests/simulation/energyPGP/pgp_lj_debug.py
"""
Debug LJ energy calculation in PGP.
"""

import pygcmc
from pygcmc import MCState, MCAtom, MCResidue, MCForceField
from pygcmc import setPGPParameters, initializePMEParameters, precomputeGridPotential
from pygcmc import computeSystemEnergyPGP


def test_lj_debug():
    """Debug LJ energy calculation."""
    print("\n=== LJ Energy Debug ===")
    
    # Create system with two LJ particles
    state = MCState()
    state.info.box = [5.0, 5.0, 5.0]
    state.info.cutoff = 1.2
    
    # Force field
    ff = MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.ljEps = [1.0]  # epsilon = 1.0 kJ/mol
    ff.ljSigma = [0.34]  # sigma = 0.34 nm
    state.forcefield = ff
    
    # Two atoms at r = 2^(1/6) * sigma (LJ minimum)
    r_min = 2**(1/6) * 0.34  # ~0.382 nm
    
    atoms = []
    
    atom1 = MCAtom()
    atom1.x = 2.5
    atom1.y = 2.5
    atom1.z = 2.5
    atom1.charge = 0.0
    atom1.type = 0
    atoms.append(atom1)
    
    atom2 = MCAtom()
    atom2.x = 2.5 + r_min
    atom2.y = 2.5
    atom2.z = 2.5
    atom2.charge = 0.0
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
        res.energy_vdw = 0.0  # Initialize
        residues.append(res)
    
    state.residues = residues
    state.activeResidueCount = 2
    
    print(f"Two LJ particles at distance r = {r_min:.4f} nm")
    print(f"Expected energy at minimum: -1.0 kJ/mol")
    
    # Initialize and calculate
    alpha = 2.5
    mesh_size = [32, 32, 32]
    initializePMEParameters(state.info.cutoff, state.info.box, alpha)
    setPGPParameters(alpha, mesh_size, state.info.cutoff, mesh_size, 4, 1e-6)
    precomputeGridPotential(state, fixed_only=True)
    computeSystemEnergyPGP(state)
    
    # Check individual residue energies
    print(f"\nResidue 0 LJ energy: {state.residues[0].energy_vdw:.6f} kJ/mol")
    print(f"Residue 1 LJ energy: {state.residues[1].energy_vdw:.6f} kJ/mol")
    
    total_lj = state.residues[0].energy_vdw + state.residues[1].energy_vdw
    print(f"Total LJ energy: {total_lj:.6f} kJ/mol")
    
    # The issue: each residue gets the full interaction energy
    # This causes double counting when we sum over residues
    print("\n⚠️  Each residue stores the full pair energy!")
    print("This leads to double counting when summing over residues.")
    
    # Correct way: take half of the sum, or just one residue's energy
    correct_lj = total_lj / 2.0
    print(f"\nCorrected LJ energy: {correct_lj:.6f} kJ/mol")
    
    if abs(correct_lj - (-1.0)) < 0.01:
        print("✅ Corrected value matches expected minimum energy")
    
    return total_lj, correct_lj


if __name__ == "__main__":
    test_lj_debug()
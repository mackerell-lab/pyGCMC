# tests/simulation/energyPGP/pgp_simple_test.py
"""
Simple test to verify PGP is working correctly after the fix.
"""

import pygcmc
from pygcmc import MCState, MCAtom, MCResidue, MCForceField
from pygcmc import setPGPParameters, setPMEParameters, initializePMEParameters
from pygcmc import precomputeGridPotential, computeSystemEnergyPGP, computeSystemEnergyPME


def test_pgp_simple():
    """Simple test with two charged particles."""
    print("\n=== Simple PGP Test ===")
    
    # Create system
    state = MCState()
    state.info.box = [5.0, 5.0, 5.0]
    state.info.cutoff = 1.2
    
    # Force field
    ff = MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.ljEps = [0.0]
    ff.ljSigma = [0.3]
    state.forcefield = ff
    
    # Two atoms
    atoms = []
    
    atom1 = MCAtom()
    atom1.x = 2.5
    atom1.y = 2.5
    atom1.z = 2.5
    atom1.charge = 1.0
    atom1.type = 0
    atoms.append(atom1)
    
    atom2 = MCAtom()
    atom2.x = 3.0  # 0.5 nm away
    atom2.y = 2.5
    atom2.z = 2.5
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
    
    # Initialize
    alpha = 2.5
    mesh_size = [32, 32, 32]
    
    # PGP calculation
    print("\nPGP calculation:")
    initializePMEParameters(state.info.cutoff, state.info.box, alpha)
    setPGPParameters(alpha, mesh_size, state.info.cutoff, mesh_size, 4, 1e-6)
    precomputeGridPotential(state, fixed_only=True)
    computeSystemEnergyPGP(state)
    
    pgp_real = state.ewald_energy.get('real_space', 0.0)
    pgp_recip = state.ewald_energy.get('reciprocal', 0.0)
    pgp_self = state.ewald_energy.get('self', 0.0)
    pgp_total = state.ewald_energy.get('total', 0.0)
    
    print(f"  Real-space: {pgp_real:.4f} kJ/mol")
    print(f"  Reciprocal: {pgp_recip:.4f} kJ/mol")
    print(f"  Self:       {pgp_self:.4f} kJ/mol")
    print(f"  Total:      {pgp_total:.4f} kJ/mol")
    
    # PME calculation for comparison
    print("\nPME calculation:")
    state2 = MCState()
    state2.info.box = [5.0, 5.0, 5.0]
    state2.info.cutoff = 1.2
    state2.forcefield = ff
    state2.atoms = [atom1, atom2]
    state2.activeAtomCount = 2
    state2.residues = residues
    state2.activeResidueCount = 2
    
    setPMEParameters(alpha, mesh_size, 4, 1e-6)
    initializePMEParameters(state2.info.cutoff, state2.info.box, alpha)
    computeSystemEnergyPME(state2)
    
    pme_real = state2.ewald_energy.get('real_space', 0.0)
    pme_recip = state2.ewald_energy.get('reciprocal', 0.0)
    pme_self = state2.ewald_energy.get('self', 0.0)
    pme_total = state2.ewald_energy.get('total', 0.0)
    
    print(f"  Real-space: {pme_real:.4f} kJ/mol")
    print(f"  Reciprocal: {pme_recip:.4f} kJ/mol")
    print(f"  Self:       {pme_self:.4f} kJ/mol")
    print(f"  Total:      {pme_total:.4f} kJ/mol")
    
    # Compare
    print("\nComparison:")
    print(f"Real-space diff: {abs(pgp_real - pme_real):.6f} kJ/mol")
    print(f"Reciprocal diff: {abs(pgp_recip - pme_recip):.6f} kJ/mol")
    print(f"Total diff:      {abs(pgp_total - pme_total):.6f} kJ/mol")
    
    # Check if real-space is correct after fix
    assert abs(pgp_real - pme_real) < 1.0, f"Real-space mismatch: PGP={pgp_real}, PME={pme_real}"
    print("\n✅ PGP real-space calculation is working correctly!")


if __name__ == "__main__":
    test_pgp_simple()
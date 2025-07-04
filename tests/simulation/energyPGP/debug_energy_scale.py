# tests/simulation/energyPGP/debug_energy_scale.py
"""
Debug the energy scale issue in PGP tests.
"""

import pygcmc
from pygcmc import MCState, MCAtom, MCResidue, MCForceField
from pygcmc import setPGPParameters, setPMEParameters, initializePMEParameters
from pygcmc import precomputeGridPotential, computeSystemEnergyPGP, computeSystemEnergyPME


def test_simple_two_particles():
    """Test with just two particles to understand energy scale."""
    print("\n=== Simple Two Particle Test ===")
    
    # Create minimal system
    state = MCState()
    state.info.box = [10.0, 10.0, 10.0]  # Larger box
    state.info.cutoff = 1.2
    
    # Force field
    ff = MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.ljEps = [0.0]  # No LJ
    ff.ljSigma = [0.3]
    state.forcefield = ff
    
    # Two particles far apart
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
    
    # Test parameters
    alpha = 2.5
    mesh_size = [32, 32, 32]
    
    print("\nSystem setup:")
    print(f"  Box: {state.info.box}")
    print(f"  Distance: 1.0 nm")
    print(f"  Charges: +1, -1")
    print(f"  Alpha: {alpha}")
    print(f"  Mesh: {mesh_size}")
    
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
    
    print(f"  Real-space: {pgp_real:.6f} kJ/mol")
    print(f"  Reciprocal: {pgp_recip:.6f} kJ/mol")
    print(f"  Self:       {pgp_self:.6f} kJ/mol")
    print(f"  Total:      {pgp_total:.6f} kJ/mol")
    
    # PME calculation
    print("\nPME calculation:")
    state2 = MCState()
    state2.info.box = [10.0, 10.0, 10.0]
    state2.info.cutoff = 1.2
    state2.forcefield = ff
    state2.atoms = atoms
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
    
    print(f"  Real-space: {pme_real:.6f} kJ/mol")
    print(f"  Reciprocal: {pme_recip:.6f} kJ/mol")
    print(f"  Self:       {pme_self:.6f} kJ/mol")
    print(f"  Total:      {pme_total:.6f} kJ/mol")
    
    # Expected values
    print("\nExpected (analytical):")
    import math
    kC = 138.935456  # Coulomb constant
    r = 1.0  # nm
    
    # Direct Coulomb
    direct = -kC / r
    print(f"  Direct Coulomb: {direct:.6f} kJ/mol")
    
    # Real-space with erfc
    erfc_val = math.erfc(alpha * r)
    real_expected = -erfc_val * kC / r
    print(f"  Real-space (erfc): {real_expected:.6f} kJ/mol")
    
    # Self energy
    self_expected = -2 * alpha * kC / math.sqrt(math.pi)
    print(f"  Self energy: {self_expected:.6f} kJ/mol")


def test_lj_cutoff():
    """Test LJ with cutoff to see if there's a maximum energy limit."""
    print("\n\n=== LJ Cutoff Test ===")
    
    state = MCState()
    state.info.box = [5.0, 5.0, 5.0]
    state.info.cutoff = 1.2
    
    ff = MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.ljEps = [1.0]
    ff.ljSigma = [0.34]
    state.forcefield = ff
    
    # Test very short distance
    distance = 0.068  # nm (0.2 * sigma)
    
    atoms = []
    atom1 = MCAtom()
    atom1.x = 2.5
    atom1.y = 2.5
    atom1.z = 2.5
    atom1.charge = 0.0
    atom1.type = 0
    atoms.append(atom1)
    
    atom2 = MCAtom()
    atom2.x = 2.5 + distance
    atom2.y = 2.5
    atom2.z = 2.5
    atom2.charge = 0.0
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
    
    print(f"Testing LJ at distance = {distance} nm")
    
    # Calculate
    alpha = 2.5
    mesh_size = [32, 32, 32]
    initializePMEParameters(state.info.cutoff, state.info.box, alpha)
    setPGPParameters(alpha, mesh_size, state.info.cutoff, mesh_size, 4, 1e-6)
    precomputeGridPotential(state, fixed_only=True)
    computeSystemEnergyPGP(state)
    
    # Get LJ energy
    lj_sum = sum(res.energy_vdw for res in state.residues if res.active)
    lj_energy = lj_sum / 2.0  # Correct for double counting
    
    # Expected
    r_ratio = 0.34 / distance
    expected = 4 * 1.0 * (r_ratio**12 - r_ratio**6)
    
    print(f"PGP LJ energy: {lj_energy:.2e} kJ/mol")
    print(f"Expected:      {expected:.2e} kJ/mol")
    print(f"Ratio:         {lj_energy/expected:.6f}")
    
    if abs(lj_energy - expected) > 0.01 * expected:
        print("⚠️  Large discrepancy detected!")
        print("Possible reasons:")
        print("  1. Energy cutoff at 10^6 kJ/mol")
        print("  2. Numerical precision issues")
        print("  3. Switching function applied")


if __name__ == "__main__":
    test_simple_two_particles()
    test_lj_cutoff()
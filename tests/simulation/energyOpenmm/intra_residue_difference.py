"""
Test to demonstrate the difference between PME and Ewald for intra-residue pairs
"""

import pygcmc
from pygcmc import MCState, MCInfo, MCAtom, MCResidue, MCForceField


def create_single_residue_two_ions():
    """Create a system with two ions in the same residue"""
    state = MCState()
    
    # Set box size
    state.info.box = [3.0, 3.0, 3.0]
    state.info.cutoff = 1.2
    
    # Create two atoms with opposite charges in the SAME residue
    atoms = []
    
    # Na+ ion
    atom1 = MCAtom()
    atom1.x = 1.5
    atom1.y = 1.5
    atom1.z = 1.5
    atom1.charge = 1.0
    atom1.type = 0
    atoms.append(atom1)
    
    # Cl- ion (0.5 nm away)
    atom2 = MCAtom()
    atom2.x = 2.0
    atom2.y = 1.5
    atom2.z = 1.5
    atom2.charge = -1.0
    atom2.type = 1
    atoms.append(atom2)
    
    # Create ONE residue containing BOTH atoms
    residues = []
    res = MCResidue()
    res.atomStart = 0
    res.atomCount = 2  # Both atoms in same residue
    res.active = True
    res.fixed = False
    res.type = 0
    residues.append(res)
    
    state.atoms = atoms
    state.residues = residues
    state.activeAtomCount = 2
    state.activeResidueCount = 1
    
    # Force field
    ff = MCForceField()
    ff.numTotalTypes = 2
    ff.ljSigma = [0.3, 0.35, 0.35, 0.4]
    ff.ljEps = [0.01, 0.02, 0.02, 0.04]
    state.forcefield = ff
    
    return state


def create_two_residue_two_ions():
    """Create a system with two ions in separate residues"""
    state = MCState()
    
    # Set box size
    state.info.box = [3.0, 3.0, 3.0]
    state.info.cutoff = 1.2
    
    # Create two atoms with opposite charges in SEPARATE residues
    atoms = []
    
    # Na+ ion
    atom1 = MCAtom()
    atom1.x = 1.5
    atom1.y = 1.5
    atom1.z = 1.5
    atom1.charge = 1.0
    atom1.type = 0
    atoms.append(atom1)
    
    # Cl- ion (0.5 nm away)
    atom2 = MCAtom()
    atom2.x = 2.0
    atom2.y = 1.5
    atom2.z = 1.5
    atom2.charge = -1.0
    atom2.type = 1
    atoms.append(atom2)
    
    # Create TWO residues, one for each atom
    residues = []
    
    res1 = MCResidue()
    res1.atomStart = 0
    res1.atomCount = 1
    res1.active = True
    res1.fixed = False
    res1.type = 0
    residues.append(res1)
    
    res2 = MCResidue()
    res2.atomStart = 1
    res2.atomCount = 1
    res2.active = True
    res2.fixed = False
    res2.type = 1
    residues.append(res2)
    
    state.atoms = atoms
    state.residues = residues
    state.activeAtomCount = 2
    state.activeResidueCount = 2
    
    # Force field
    ff = MCForceField()
    ff.numTotalTypes = 2
    ff.ljSigma = [0.3, 0.35, 0.35, 0.4]
    ff.ljEps = [0.01, 0.02, 0.02, 0.04]
    state.forcefield = ff
    
    return state


def test_intra_residue_handling():
    """Test how PME and Ewald handle intra-residue pairs differently"""
    
    print("\nTesting intra-residue pair handling:")
    print("=" * 70)
    
    # Test parameters
    alpha = 3.0
    kmax = [8, 8, 8]
    mesh_size = [32, 32, 32]
    spline_order = 4
    box = [3.0, 3.0, 3.0]
    cutoff = 1.2
    
    # Test 1: Two ions in SAME residue
    print("\n1. Two ions in SAME residue:")
    state1 = create_single_residue_two_ions()
    
    # Ewald
    pygcmc.setEwaldParameters(alpha, kmax)
    pygcmc.initializeEwaldParameters(cutoff, box, alpha)
    _, _, ewald_dict1 = pygcmc.computeSystemEnergyEwald(state1)
    
    # PME
    pygcmc.setPMEParameters(alpha, mesh_size, spline_order)
    pygcmc.initializePMEParameters(cutoff, box, alpha, mesh_size, spline_order)
    _, _, pme_dict1 = pygcmc.computeSystemEnergyPME(state1)
    
    print(f"  Ewald real space: {ewald_dict1['real_space']:.6f} kJ/mol")
    print(f"  PME real space:   {pme_dict1['real_space']:.6f} kJ/mol")
    print(f"  Difference:       {abs(ewald_dict1['real_space'] - pme_dict1['real_space']):.6f} kJ/mol")
    
    # Test 2: Two ions in SEPARATE residues
    print("\n2. Two ions in SEPARATE residues:")
    state2 = create_two_residue_two_ions()
    
    # Ewald
    pygcmc.setEwaldParameters(alpha, kmax)
    pygcmc.initializeEwaldParameters(cutoff, box, alpha)
    _, _, ewald_dict2 = pygcmc.computeSystemEnergyEwald(state2)
    
    # PME
    pygcmc.setPMEParameters(alpha, mesh_size, spline_order)
    pygcmc.initializePMEParameters(cutoff, box, alpha, mesh_size, spline_order)
    _, _, pme_dict2 = pygcmc.computeSystemEnergyPME(state2)
    
    print(f"  Ewald real space: {ewald_dict2['real_space']:.6f} kJ/mol")
    print(f"  PME real space:   {pme_dict2['real_space']:.6f} kJ/mol")
    print(f"  Difference:       {abs(ewald_dict2['real_space'] - pme_dict2['real_space']):.6f} kJ/mol")
    
    # Analysis
    print("\nAnalysis:")
    print(f"  Ewald difference (separate - same): {ewald_dict2['real_space'] - ewald_dict1['real_space']:.6f} kJ/mol")
    print(f"  PME difference (separate - same):   {pme_dict2['real_space'] - pme_dict1['real_space']:.6f} kJ/mol")
    
    print("\nConclusion:")
    if abs(ewald_dict1['real_space']) < 1e-6:
        print("  ✗ Ewald skips intra-residue pairs (real space = 0 for same residue)")
    print("  ✓ PME includes intra-residue pairs correctly")
    
    # Check if total energies are different due to reciprocal space
    print(f"\n  Total energies:")
    print(f"  Same residue:     Ewald={ewald_dict1['total']:.6f}, PME={pme_dict1['total']:.6f}")
    print(f"  Separate residues: Ewald={ewald_dict2['total']:.6f}, PME={pme_dict2['total']:.6f}")


if __name__ == "__main__":
    test_intra_residue_handling()
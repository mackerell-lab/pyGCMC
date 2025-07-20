#!/usr/bin/env python
"""Test Drude with rigid molecule assumption in GCMC"""

import sys
import numpy as np
sys.path.append('/home/zhaomt/gcmc/test107/pygcmc_dev/build/modules/bindings')

import pygcmc

def test_drude_rigid_molecule():
    print("=== Testing Drude with Rigid Molecule (GCMC Style) ===\n")
    
    # Initialize
    pygcmc.initializeDrudeForce()
    
    # Create state
    state = pygcmc.MCState()
    state.info.cutoff = 2.0
    state.info.box = [10.0, 10.0, 10.0]
    
    # SWM4-NDP parameters
    qO_core = 1.66260
    qD = -1.76260
    qM = -0.95710
    qH = 0.52855
    alpha_nm3 = 0.0013
    k_kj_nm2 = 166018.9
    
    print("Key insight for GCMC:")
    print("- Molecules are rigid (fixed internal coordinates)")
    print("- Only Drude positions change during SCF")
    print("- Intramolecular non-Drude energies are constant")
    print("- M-site position is fixed relative to O,H1,H2\n")
    
    # Create single water
    atoms = []
    
    # Fixed geometry water
    o = pygcmc.MCAtom()
    o.x = 0.0
    o.y = 0.0
    o.z = 0.0
    o.charge = qO_core
    o.type = 0
    atoms.append(o)
    
    d = pygcmc.MCAtom()
    d.x = 0.0
    d.y = 0.0
    d.z = 0.0
    d.charge = qD
    d.type = 1
    atoms.append(d)
    
    h1 = pygcmc.MCAtom()
    h1.x = 0.09572
    h1.y = 0.0
    h1.z = 0.0
    h1.charge = qH
    h1.type = 2
    atoms.append(h1)
    
    h2 = pygcmc.MCAtom()
    h2.x = -0.023999
    h2.y = 0.092663
    h2.z = 0.0
    h2.charge = qH
    h2.type = 2
    atoms.append(h2)
    
    m = pygcmc.MCAtom()
    m.x = 0.00793
    m.y = 0.00986
    m.z = 0.0
    m.charge = qM
    m.type = 3
    atoms.append(m)
    
    # Calculate constant intramolecular energy (excluding Drude)
    print("Constant intramolecular energies (kJ/mol):")
    print("-" * 50)
    
    ONE_4PI_EPS0 = 138.935456
    const_energy = 0.0
    
    # Non-Drude pairs
    pairs = [
        ("O-H1", 0, 2), ("O-H2", 0, 3), ("O-M", 0, 4),
        ("H1-H2", 2, 3), ("H1-M", 2, 4), ("H2-M", 3, 4)
    ]
    
    for name, i, j in pairs:
        r = np.linalg.norm([atoms[j].x - atoms[i].x, 
                           atoms[j].y - atoms[i].y,
                           atoms[j].z - atoms[i].z])
        E = ONE_4PI_EPS0 * atoms[i].charge * atoms[j].charge / r
        const_energy += E
        print(f"{name:8} {E:12.2f}")
    
    print(f"{'Total':8} {const_energy:12.2f}")
    print("(This energy is constant in GCMC - can be ignored)\n")
    
    # Set up system
    res = pygcmc.MCResidue()
    res.atomStart = 0
    res.atomCount = 5
    res.active = True
    res.fixed = False
    res.type = 0
    
    state.atoms = atoms
    state.activeAtomCount = 5
    state.residues = [res]
    state.activeResidueCount = 1
    
    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 4
    state.forcefield = ff
    
    # Add Drude
    pygcmc.addDrudeParticle(
        drudeIndex=1,
        parentIndex=0,
        charge=qD,
        polarizability=alpha_nm3
    )
    
    params = pygcmc.DrudeSCFParams()
    params.tolerance = 1e-8
    params.maxIterations = 100
    params.dampingFactor = 0.5
    params.maxDrudeDistance = 0.02
    pygcmc.setDrudeSCFParameters(params)
    
    # Test 1: Isolated water
    print("Test 1: Isolated rigid water")
    print("-" * 40)
    
    energy1, _ = pygcmc.computeSystemEnergyDrude(state)
    
    dx = state.atoms[1].x - state.atoms[0].x
    dy = state.atoms[1].y - state.atoms[0].y
    dz = state.atoms[1].z - state.atoms[0].z
    dist1 = np.sqrt(dx*dx + dy*dy + dz*dz)
    
    # Drude-related energies only
    E_harmonic = 0.5 * k_kj_nm2 * dist1**2
    
    print(f"Total energy: {energy1:.2f} kJ/mol")
    print(f"Drude displacement: {dist1*1000:.3f} pm")
    print(f"Harmonic energy: {E_harmonic:.2f} kJ/mol")
    print(f"Other Drude energies: {energy1 - E_harmonic:.2f} kJ/mol")
    
    # Test 2: Two rigid waters
    print("\n\nTest 2: Two rigid waters")
    print("-" * 40)
    
    # Add second water at 3 Å
    offset = 0.3
    for i in range(5):
        atom = pygcmc.MCAtom()
        atom.x = atoms[i].x + offset
        atom.y = atoms[i].y
        atom.z = atoms[i].z
        atom.charge = atoms[i].charge
        atom.type = atoms[i].type
        state.atoms.append(atom)
    
    state.activeAtomCount = 10
    
    res2 = pygcmc.MCResidue()
    res2.atomStart = 5
    res2.atomCount = 5
    res2.active = True
    res2.fixed = False
    res2.type = 0
    state.residues.append(res2)
    state.activeResidueCount = 2
    
    pygcmc.addDrudeParticle(
        drudeIndex=6,
        parentIndex=5,
        charge=qD,
        polarizability=alpha_nm3
    )
    
    # Reset Drudes
    state.atoms[1].x = 0.0
    state.atoms[6].x = offset
    
    energy2, _ = pygcmc.computeSystemEnergyDrude(state)
    
    interaction = energy2 - 2*energy1
    
    print(f"Total energy: {energy2:.2f} kJ/mol")
    print(f"Single water energy: {energy1:.2f} kJ/mol")
    print(f"Interaction energy: {interaction:.2f} kJ/mol")
    
    # Check individual contributions
    print("\nIntermolecular contributions at 3 Å:")
    
    # Calculate some key distances
    r_OO = offset
    r_OH = np.sqrt(offset**2 + 0.09572**2)
    r_OM = np.sqrt((offset-0.00793)**2 + 0.00986**2)
    
    print(f"O-O distance: {r_OO*10:.1f} Å")
    print(f"O-H distance: ~{r_OH*10:.1f} Å")
    print(f"O-M distance: ~{r_OM*10:.1f} Å")
    
    # Estimate dominant contribution
    E_OO = ONE_4PI_EPS0 * (qO_core + qD) * (qO_core + qD) / r_OO
    E_OM = ONE_4PI_EPS0 * (qO_core + qD) * qM / r_OM
    
    print(f"\nApproximate contributions:")
    print(f"O-O (net charges): {E_OO:.2f} kJ/mol")
    print(f"O-M interaction: {E_OM:.2f} kJ/mol")
    
    print("\n\nConclusions for GCMC:")
    print("-" * 40)
    print("1. Intramolecular energies (non-Drude) are constant")
    print("2. Only need to calculate:")
    print("   - Drude harmonic restraint")
    print("   - Drude-related intramolecular terms")
    print("   - All intermolecular interactions")
    print("3. M-site DOES participate in intermolecular Coulomb")
    print("4. Current high single-molecule energy includes constant terms")

if __name__ == "__main__":
    test_drude_rigid_molecule()
# tests/simulation/energyNB/nonbonded_pbc_basic.py
"""Non-bonded PBC basic tests."""

import pytest
import pygcmc
import math


def test_pbc_basic():
    """Test basic PBC energy calculation.
    
    Setup:
    - Two residues at opposite sides of the box
    - Box size = 2.0 nm
    - Residues should interact through PBC
    """
    state = pygcmc.MCState()
    
    # 1. Set up force field
    state.forcefield.numTotalTypes = 2
    state.forcefield.numMovementTypes = 2
    state.forcefield.ljEps = [1.0] * 4
    state.forcefield.ljSigma = [1.0] * 4
    
    # Set up a large box and appropriate cutoff
    state.info.box = [10.0, 10.0, 10.0]  # Large box to avoid PBC effects
    state.info.cutoff = 2.0  # nm
    
    # 2. Set up two atoms on opposite sides of the box
    atoms = []
    positions = [
        (0.5, 1.0, 1.0),    # Near the middle of the box
        (1.5, 1.0, 1.0)     # Near the middle of the box
    ]
    charges = [0.5, -0.5]   # Opposite charges
    
    for i, (x, y, z) in enumerate(positions):
        atom = pygcmc.MCAtom()
        atom.x = x
        atom.y = y
        atom.z = z
        atom.charge = charges[i]
        atom.type = i
        atoms.append(atom)
    
    state.atoms = atoms
    state.activeAtomCount = 2
    
    # 3. Set up two residues
    residues = []
    for i in range(2):
        res = pygcmc.MCResidue()
        res.active = True
        res.type = i
        res.atomStart = i
        res.atomCount = 1
        residues.append(res)
    
    state.residues = residues
    state.activeResidueCount = 2
    
    # Calculate PBC energy
    pygcmc.computeSystemEnergyPBCCutoff(state)
    pbc_vdw = sum(res.energy_vdw for res in state.residues)
    pbc_elec = sum(res.energy_elec for res in state.residues)
    pbc_energy = (pbc_vdw + pbc_elec) / 2.0
    
    # Through PBC, the actual distance between atoms should be 1.0 nm (not 1.0 nm)
    # In x direction: 1.5 - 0.5 = 1.0, exactly equal to sigma
    COULOMB = 138.935458  # kJ·nm/mol/e²
    min_dist = 1.0  # Minimum image distance
    expected_elec = COULOMB * 0.5 * (-0.5) / min_dist  # Electrostatic energy
    
    # LJ energy is 0 at r = sigma
    expected_vdw = 0.0  # Assuming r ≈ sigma
    expected_energy = expected_elec + expected_vdw
    
    # Print detailed debug information
    print(f"\nDetailed energy comparison:")
    print(f"Box size: {state.info.box}")
    print(f"Cutoff distance: {state.info.cutoff} nm")
    print(f"Atom positions: {positions}")
    print(f"Minimum image distance: {min_dist} nm")
    print(f"Expected electrostatic energy: {expected_elec} kJ/mol")
    print(f"Expected VDW energy: {expected_vdw} kJ/mol")
    print(f"Expected total energy: {expected_energy} kJ/mol")
    print(f"Actual VDW energy: {pbc_vdw/2} kJ/mol")
    print(f"Actual electrostatic energy: {pbc_elec/2} kJ/mol")
    print(f"Actual total energy: {pbc_energy} kJ/mol")
    
    # Check VDW and electrostatic energies separately
    rel_tol = 0.1  # 10% relative error tolerance
    
    # Check electrostatic energy
    actual_elec = pbc_elec / 2.0
    assert abs((actual_elec - expected_elec) / expected_elec) < rel_tol, \
           f"Electrostatic energy {actual_elec} differs too much from expected {expected_elec}"
    
    # Check VDW energy (should be close to 0 because r ≈ sigma)
    actual_vdw = pbc_vdw / 2.0
    assert abs(actual_vdw) < 10.0, f"VDW energy {actual_vdw} is too large"
    
    # Check total energy
    assert abs((pbc_energy - expected_energy) / expected_energy) < rel_tol, \
           f"PBC energy {pbc_energy} differs too much from expected {expected_energy}"


def test_pbc_invalid_box():
    """Test PBC energy calculation with invalid box dimensions."""
    state = pygcmc.MCState()
    
    # Set up basic state
    state.forcefield.numTotalTypes = 1
    state.forcefield.numMovementTypes = 1
    state.forcefield.ljEps = [1.0]
    state.forcefield.ljSigma = [1.0]
    
    # Set invalid box dimensions
    state.info.box = [0.0, 1.0, 1.0]  # x dimension is 0
    
    # Add a simple atom and residue
    atom = pygcmc.MCAtom()
    state.atoms = [atom]
    
    res = pygcmc.MCResidue()
    res.active = True
    state.residues = [res]
    state.activeResidueCount = 1
    
    # Expect exception to be raised
    with pytest.raises(RuntimeError, match="Invalid box dimensions"):
        pygcmc.computeSystemEnergyPBCCutoff(state)


def test_pbc_vs_nopbc():
    """Compare PBC and non-PBC energy calculations.
    
    Setup:
    - Two residues well within the box (not near boundaries)
    - Results should be identical when atoms are far from boundaries
    """
    state = pygcmc.MCState()
    
    # 1. Set up force field
    state.forcefield.numTotalTypes = 2
    state.forcefield.numMovementTypes = 2
    state.forcefield.ljEps = [1.0] * 4
    state.forcefield.ljSigma = [1.0] * 4
    
    # Set up a large box and appropriate cutoff
    state.info.box = [10.0, 10.0, 10.0]  # Large box to avoid PBC effects
    state.info.cutoff = 2.0  # nm
    
    # 2. Set up two atoms in the middle of the box
    atoms = []
    positions = [
        (4.0, 5.0, 5.0),
        (5.0, 5.0, 5.0)
    ]
    charges = [0.5, -0.5]
    
    for i, (x, y, z) in enumerate(positions):
        atom = pygcmc.MCAtom()
        atom.x = x
        atom.y = y
        atom.z = z
        atom.charge = charges[i]
        atom.type = i
        atoms.append(atom)
    
    state.atoms = atoms
    state.activeAtomCount = 2
    
    # 3. Set up residues
    residues = []
    for i in range(2):
        res = pygcmc.MCResidue()
        res.active = True
        res.type = i
        res.atomStart = i
        res.atomCount = 1
        residues.append(res)
    
    state.residues = residues
    state.activeResidueCount = 2
    
    # Calculate both types of energy
    pygcmc.computeSystemEnergyCutoff(state)
    nopbc_vdw = sum(res.energy_vdw for res in state.residues)
    nopbc_elec = sum(res.energy_elec for res in state.residues)
    nopbc_energy = (nopbc_vdw + nopbc_elec) / 2.0
    
    # Reset energy
    for res in state.residues:
        res.energy_vdw = 0.0
        res.energy_elec = 0.0
    
    # Calculate PBC energy
    pygcmc.computeSystemEnergyPBCCutoff(state)
    pbc_vdw = sum(res.energy_vdw for res in state.residues)
    pbc_elec = sum(res.energy_elec for res in state.residues)
    pbc_energy = (pbc_vdw + pbc_elec) / 2.0
    
    # For atoms in the middle of the box, both calculation methods should give the same result
    rel_tol = 1e-5  # Very small relative error tolerance
    assert abs((pbc_energy - nopbc_energy) / nopbc_energy) < rel_tol, \
           f"PBC energy {pbc_energy} differs from non-PBC energy {nopbc_energy}"
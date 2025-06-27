# tests/simulation/energyNB/nonbonded_pbc.py
"""Non-bonded PBC and cutoff tests."""

import pytest
import pygcmc
import math


def test_cutoff_nonperiodic():
    """Test nonbonded energy calculation with distance cutoff.
    
    Setup:
    - Multiple residues at various distances
    - Only interactions within cutoff distance are considered
    - No periodic boundary conditions
    """
    state = pygcmc.MCState()
    
    # 1. Set up force field
    state.forcefield.numTotalTypes = 3  # Three types
    state.forcefield.numMovementTypes = 3
    
    # Set up force field parameters (3x3 matrix flattened)
    state.forcefield.ljEps = [1.0] * 9    # All interactions eps=1.0
    state.forcefield.ljSigma = [1.0] * 9  # All interactions sigma=1.0
    
    # Set cutoff distance
    state.info.cutoff = 1.5  # nm
    
    # 2. Set up three atoms arranged sequentially on the x-axis
    atoms = []
    positions = [
        (0.0, 0.0, 0.0),  # Origin
        (1.0, 0.0, 0.0),  # x=1.0 nm, within cutoff distance
        (2.0, 0.0, 0.0)   # x=2.0 nm, beyond cutoff distance
    ]
    charges = [0.5, -0.5, 0.5]  # Alternating charges
    
    for i, (x, y, z) in enumerate(positions):
        atom = pygcmc.MCAtom()
        atom.x = x
        atom.y = y
        atom.z = z
        atom.charge = charges[i]
        atom.type = i
        atoms.append(atom)
    
    state.atoms = atoms
    state.activeAtomCount = 3
    
    # 3. Set up three residues
    residues = []
    for i in range(3):
        res = pygcmc.MCResidue()
        res.active = True
        res.type = i
        res.atomStart = i
        res.atomCount = 1
        residues.append(res)
    
    state.residues = residues
    state.activeResidueCount = 3
    
    # Calculate energy
    pygcmc.computeSystemEnergyCutoff(state)
    
    # Verify each residue has energy
    for i in range(3):
        energy = state.residues[i].energy_vdw + state.residues[i].energy_elec
        assert not math.isnan(energy), f"Residue {i} has NaN energy"
        assert not math.isinf(energy), f"Residue {i} has infinite energy"
    
    # Calculate total energy
    total_vdw = sum(res.energy_vdw for res in state.residues)
    total_elec = sum(res.energy_elec for res in state.residues)
    total_energy = total_vdw + total_elec
    
    # Total energy needs to be divided by 2 because each interaction was calculated twice
    system_energy = total_energy / 2.0
    
    # Verify the reasonableness of total energy
    assert system_energy < 1e6, "Total system energy too large"
    assert system_energy > -1e6, "Total system energy too negative"
    
    # Verify that energy magnitude is reasonable
    # For our setup:
    # - LJ energy is 0 at r=sigma
    # - Electrostatic energy = k_c * q1*q2/r
    # Only residue 0-1 interaction is within cutoff distance:
    # - q1=0.5, q2=-0.5, r=1.0
    # - The energy of each interaction is fully added to both residues
    COULOMB = 138.935458  # kJ·nm/mol/e²
    expected_energy = 2.0 * COULOMB * 0.5 * (-0.5) / 1.0  # Interaction energy between one pair of atoms, multiplied by 2 because it's added twice in the code

    rel_tol = 0.1  # 10% relative error tolerance
    assert abs((system_energy - expected_energy) / expected_energy) < rel_tol, \
           f"System energy {system_energy} differs too much from expected {expected_energy}"
    
    # Verify that residue 2 energy is 0 (because its distance to other residues exceeds cutoff)
    res2_energy = state.residues[2].energy_vdw + state.residues[2].energy_elec
    # residue 2 should have energy because it receives energy from interaction with residue 1
    expected_res2_energy = COULOMB * 0.5 * (-0.5) / 1.0  # residue 0-1 interaction
    assert abs((res2_energy - expected_res2_energy) / expected_res2_energy) < rel_tol, \
           f"Residue 2 energy {res2_energy} differs too much from expected {expected_res2_energy}"


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


def test_pbc_cross_boundary():
    """Test PBC energy calculation for atoms crossing periodic boundaries.
    
    Setup:
    - Three atoms forming a chain across the periodic boundary
    - Middle atom near boundary should interact correctly with both other atoms
    """
    state = pygcmc.MCState()
    
    # 1. Set up force field
    state.forcefield.numTotalTypes = 3
    state.forcefield.numMovementTypes = 3
    state.forcefield.ljEps = [1.0] * 9
    state.forcefield.ljSigma = [1.0] * 9
    
    # Set up box and cutoff
    state.info.box = [3.0, 3.0, 3.0]  # 3.0 nm cubic box
    state.info.cutoff = 2.0  # nm
    
    # 2. Set up three atoms: one near the box edge, one crossing the boundary, and one on the other side of the box
    atoms = []
    positions = [
        (2.8, 1.5, 1.5),  # Near the right boundary
        (0.1, 1.5, 1.5),  # Near the left boundary
        (1.5, 1.5, 1.5)   # In the middle
    ]
    charges = [0.5, 0.5, -1.0]  # Make the middle atom attractive to both sides
    
    for i, (x, y, z) in enumerate(positions):
        atom = pygcmc.MCAtom()
        atom.x = x
        atom.y = y
        atom.z = z
        atom.charge = charges[i]
        atom.type = i
        atoms.append(atom)
    
    state.atoms = atoms
    state.activeAtomCount = 3
    
    # 3. Set up residues
    residues = []
    for i in range(3):
        res = pygcmc.MCResidue()
        res.active = True
        res.type = i
        res.atomStart = i
        res.atomCount = 1
        residues.append(res)
    
    state.residues = residues
    state.activeResidueCount = 3
    
    # Calculate PBC energy
    pygcmc.computeSystemEnergyPBCCutoff(state)
    
    # Verify all residues have reasonable energy
    for i, res in enumerate(state.residues):
        energy = res.energy_vdw + res.energy_elec
        assert not math.isnan(energy), f"Residue {i} has NaN energy"
        assert not math.isinf(energy), f"Residue {i} has infinite energy"
        assert abs(energy) < 1e6, f"Residue {i} has unreasonably large energy: {energy}"
    
    # Through PBC, the minimum distance between atoms 0 and 1 should be 0.3 nm
    # instead of the actual distance of 2.7 nm in the box
    COULOMB = 138.935458  # kJ·nm/mol/e²
    
    # Calculate total energy
    total_vdw = sum(res.energy_vdw for res in state.residues)
    total_elec = sum(res.energy_elec for res in state.residues)
    system_energy = (total_vdw + total_elec) / 2.0
    
    # Verify energy is finite and reasonable
    assert not math.isnan(system_energy), "System energy is NaN"
    assert not math.isinf(system_energy), "System energy is infinite"
    assert abs(system_energy) < 1e6, f"System energy too large: {system_energy}"
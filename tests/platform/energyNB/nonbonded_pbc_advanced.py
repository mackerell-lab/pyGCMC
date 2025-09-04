# tests/simulation/energyNB/nonbonded_pbc_advanced.py
"""Non-bonded PBC advanced tests."""

import pytest
import pygcmc
import math


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
"""
PGP Internal Consistency Tests

These tests verify that PGP gives consistent results for its own algorithm,
without comparing to PME values.
"""

import pytest
import pygcmc
from . import pgp_wrapper
from .pgp_wrapper import initializePMEParameters, setPGPParameters, precomputeGridPotential
from .pgp_wrapper import computeSystemEnergyPGP
from pygcmc import MCState, MCAtom, MCResidue
from pygcmc import MCForceField

def test_pgp_energy_symmetry():
    """Test that PGP gives symmetric results for symmetric configurations."""
    # Create symmetric system
    state = MCState()
    state.info.box = [6.0, 6.0, 6.0]
    state.info.cutoff = 2.0
    
    ff = MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.ljEps = [0.0]
    ff.ljSigma = [0.3]
    state.forcefield = ff
    
    # Place 4 charges in symmetric square
    positions = [
        [2.0, 2.0, 3.0],
        [4.0, 2.0, 3.0],
        [4.0, 4.0, 3.0],
        [2.0, 4.0, 3.0]
    ]
    
    atoms = []
    for i, pos in enumerate(positions):
        atom = MCAtom()
        atom.x, atom.y, atom.z = pos
        atom.charge = 1.0 if i % 2 == 0 else -1.0
        atom.type = 0
        atoms.append(atom)
    
    state.atoms = atoms
    state.activeAtomCount = 4
    
    # Create residues
    residues = []
    for i in range(4):
        res = MCResidue()
        res.active = True
        res.fixed = False
        res.atomStart = i
        res.atomCount = 1
        res.type = 0
        residues.append(res)
    
    state.residues = residues
    state.activeResidueCount = 4
    
    # Calculate energy
    alpha = 2.5
    mesh_size = [32, 32, 32]
    initializePMEParameters(state.info.cutoff, state.info.box, alpha)
    setPGPParameters(alpha, mesh_size, state.info.cutoff, mesh_size, 4, 1e-6)
    precomputeGridPotential(state)
    computeSystemEnergyPGP(state)
    
    energy1 = state.ewald_energy.get('total', 0.0)
    
    # Rotate system 90 degrees (should give same energy)
    for i, atom in enumerate(atoms):
        x, y = positions[i][0], positions[i][1]
        atom.x = 3.0 - (y - 3.0)  # Rotate around center
        atom.y = 3.0 + (x - 3.0)
    
    computeSystemEnergyPGP(state)
    energy2 = state.ewald_energy.get('total', 0.0)
    
    print(f"Original energy: {energy1:.6f}")
    print(f"Rotated energy: {energy2:.6f}")
    print(f"Difference: {abs(energy1 - energy2):.6e}")
    
    # Should be symmetric
    assert abs(energy1 - energy2) < 0.1, "PGP should give symmetric results"

def test_pgp_energy_scaling():
    """Test that PGP energy scales correctly with charge."""
    base_charge = 1.0
    multipliers = [0.5, 1.0, 2.0]
    energies = []
    
    for mult in multipliers:
        state = MCState()
        state.info.box = [5.0, 5.0, 5.0]
        state.info.cutoff = 2.0
        
        ff = MCForceField()
        ff.numTotalTypes = 1
        ff.numMovementTypes = 1
        ff.ljEps = [0.0]
        ff.ljSigma = [0.3]
        state.forcefield = ff
        
        atoms = []
        atom1 = MCAtom()
        atom1.x = 2.5
        atom1.y = 2.5
        atom1.z = 2.5
        atom1.charge = base_charge * mult
        atom1.type = 0
        atoms.append(atom1)
        
        atom2 = MCAtom()
        atom2.x = 3.5
        atom2.y = 2.5
        atom2.z = 2.5
        atom2.charge = -base_charge * mult
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
        
        alpha = 2.5
        mesh_size = [32, 32, 32]
        initializePMEParameters(state.info.cutoff, state.info.box, alpha)
        setPGPParameters(alpha, mesh_size, state.info.cutoff, mesh_size, 4, 1e-6)
        precomputeGridPotential(state)
        computeSystemEnergyPGP(state)
        
        energy = state.ewald_energy.get('total', 0.0)
        energies.append(energy)
        print(f"Charge multiplier {mult}: Energy = {energy:.6f}")
    
    # Energy should scale with charge squared
    ratio1 = energies[0] / energies[1]  # 0.5^2 = 0.25
    ratio2 = energies[2] / energies[1]  # 2.0^2 = 4.0
    
    print(f"Energy ratio (0.5x/1x): {ratio1:.3f}, expected ~0.25")
    print(f"Energy ratio (2x/1x): {ratio2:.3f}, expected ~4.0")
    
    # Allow for PGP approximation errors
    assert abs(ratio1 - 0.25) < 0.1, "Energy should scale with charge squared"
    assert abs(ratio2 - 4.0) < 0.5, "Energy should scale with charge squared"

def test_pgp_grid_independence():
    """Test that PGP gives consistent results with different grid offsets."""
    # This tests that the interpolation is working correctly
    
    energies = []
    offsets = [0.0, 0.1, 0.2]  # Small position offsets
    
    for offset in offsets:
        state = MCState()
        state.info.box = [5.0, 5.0, 5.0]
        state.info.cutoff = 2.0
        
        ff = MCForceField()
        ff.numTotalTypes = 1
        ff.numMovementTypes = 1
        ff.ljEps = [0.0]
        ff.ljSigma = [0.3]
        state.forcefield = ff
        
        atoms = []
        atom1 = MCAtom()
        atom1.x = 2.5 + offset
        atom1.y = 2.5 + offset
        atom1.z = 2.5
        atom1.charge = 1.0
        atom1.type = 0
        atoms.append(atom1)
        
        atom2 = MCAtom()
        atom2.x = 3.5 + offset
        atom2.y = 2.5 + offset
        atom2.z = 2.5
        atom2.charge = -1.0
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
        
        alpha = 2.5
        mesh_size = [32, 32, 32]
        initializePMEParameters(state.info.cutoff, state.info.box, alpha)
        setPGPParameters(alpha, mesh_size, state.info.cutoff, mesh_size, 4, 1e-6)
        precomputeGridPotential(state)
        computeSystemEnergyPGP(state)
        
        energy = state.ewald_energy.get('total', 0.0)
        energies.append(energy)
        print(f"Offset {offset}: Energy = {energy:.6f}")
    
    # Energies should be very similar (interpolation error only)
    max_diff = max(energies) - min(energies)
    avg_energy = sum(energies) / len(energies)
    rel_diff = max_diff / abs(avg_energy) if avg_energy != 0 else max_diff
    
    print(f"Max energy difference: {max_diff:.6f}")
    print(f"Relative difference: {rel_diff:.6e}")
    
    assert rel_diff < 0.01, "PGP should give consistent results with position shifts"

if __name__ == "__main__":
    test_pgp_energy_symmetry()
    test_pgp_energy_scaling() 
    test_pgp_grid_independence()

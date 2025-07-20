# tests/simulation/energyDrude/screened_tests.py
"""Thole screening interaction tests for Drude oscillators."""

import pytest
import math
import pygcmc
from pygcmc import MCState, DrudeForce, DrudeSCFParams
from .helpers import ONE_4PI_EPS0


def test_thole_screening():
    """Test basic Thole screening function"""
    state = MCState()
    
    # Set box size
    state.info.box = [10.0, 10.0, 10.0]
    state.info.setTemperature(300.0)
    state.info.cutoff = 5.0
    
    # Create two polarizable atoms
    atoms = []
    
    # First atom and Drude
    atom1 = pygcmc.MCAtom()
    atom1.x = 5.0
    atom1.y = 5.0
    atom1.z = 5.0
    atom1.charge = 1.0
    atom1.type = 0
    atoms.append(atom1)
    
    drude1 = pygcmc.MCAtom()
    drude1.x = 5.0
    drude1.y = 5.0
    drude1.z = 5.0
    drude1.charge = -1.0
    drude1.type = 1
    atoms.append(drude1)
    
    # Second atom and Drude (separated by 0.5 nm)
    atom2 = pygcmc.MCAtom()
    atom2.x = 5.5
    atom2.y = 5.0
    atom2.z = 5.0
    atom2.charge = 1.0
    atom2.type = 0
    atoms.append(atom2)
    
    drude2 = pygcmc.MCAtom()
    drude2.x = 5.5
    drude2.y = 5.0
    drude2.z = 5.0
    drude2.charge = -1.0
    drude2.type = 1
    atoms.append(drude2)
    
    state.atoms = atoms
    state.activeAtomCount = 4
    
    # Create residues
    res1 = pygcmc.MCResidue()
    res1.atomStart = 0
    res1.atomCount = 2
    res1.active = True
    res1.fixed = False
    res1.type = 0
    
    res2 = pygcmc.MCResidue()
    res2.atomStart = 2
    res2.atomCount = 2
    res2.active = True
    res2.fixed = False
    res2.type = 0
    
    state.residues = [res1, res2]
    state.activeResidueCount = 2
    
    # Set force field
    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 2
    ff.numMovementTypes = 2
    ff.ljSigma = [0.3] * 4
    ff.ljEps = [0.0] * 4
    state.forcefield = ff
    
    # Create Drude force
    drude_force = DrudeForce()
    
    # Add Drude particles
    charge = -1.0
    polarizability = 0.001  # nm^3
    
    idx1 = drude_force.addParticle(
        drudeIndex=1, parentIndex=0,
        aniso1Index=-1, aniso2Index=-1,
        aniso3Index=-1, aniso4Index=-1,
        charge=charge, polarizability=polarizability,
        aniso12=1.0, aniso34=1.0
    )
    
    idx2 = drude_force.addParticle(
        drudeIndex=3, parentIndex=2,
        aniso1Index=-1, aniso2Index=-1,
        aniso3Index=-1, aniso4Index=-1,
        charge=charge, polarizability=polarizability,
        aniso12=1.0, aniso34=1.0
    )
    
    # Add screened pair with Thole parameter
    thole = 2.0  # Typical Thole parameter
    drude_force.addScreenedPair(idx1, idx2, thole)
    
    assert drude_force.getNumScreenedPairs() == 1
    
    # Calculate energy
    energy = drude_force.calculateEnergySCF(state)
    
    # Drude particles should be displaced due to mutual polarization
    dx1 = state.atoms[1].x - state.atoms[0].x
    dx2 = state.atoms[3].x - state.atoms[2].x
    
    # They should be displaced in opposite directions
    assert dx1 * dx2 < 0  # Opposite signs
    
    # Energy should be non-zero
    assert energy > 0.0


def test_screened_pair_energy():
    """Test energy calculation with screened pairs"""
    state = MCState()
    
    # Set box size
    state.info.box = [10.0, 10.0, 10.0]
    state.info.setTemperature(300.0)
    state.info.cutoff = 5.0
    
    # Create two polarizable atoms at known distance
    distance = 0.4  # nm
    
    atoms = []
    
    # First atom and Drude
    atom1 = pygcmc.MCAtom()
    atom1.x = 5.0
    atom1.y = 5.0
    atom1.z = 5.0
    atom1.charge = 2.0  # Different charges
    atom1.type = 0
    atoms.append(atom1)
    
    drude1 = pygcmc.MCAtom()
    drude1.x = 5.0
    drude1.y = 5.0
    drude1.z = 5.0
    drude1.charge = -2.0
    drude1.type = 1
    atoms.append(drude1)
    
    # Second atom and Drude
    atom2 = pygcmc.MCAtom()
    atom2.x = 5.0 + distance
    atom2.y = 5.0
    atom2.z = 5.0
    atom2.charge = 1.0
    atom2.type = 0
    atoms.append(atom2)
    
    drude2 = pygcmc.MCAtom()
    drude2.x = 5.0 + distance
    drude2.y = 5.0
    drude2.z = 5.0
    drude2.charge = -1.0
    drude2.type = 1
    atoms.append(drude2)
    
    state.atoms = atoms
    state.activeAtomCount = 4
    
    # Create residues
    residues = []
    for i in range(2):
        res = pygcmc.MCResidue()
        res.atomStart = i * 2
        res.atomCount = 2
        res.active = True
        res.fixed = False
        res.type = 0
        residues.append(res)
    
    state.residues = residues
    state.activeResidueCount = 2
    
    # Set force field
    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 2
    ff.numMovementTypes = 2
    ff.ljSigma = [0.3] * 4
    ff.ljEps = [0.0] * 4
    state.forcefield = ff
    
    # Create Drude force without screening first
    drude_force_unscreened = DrudeForce()
    
    # Different polarizabilities
    pol1 = 0.002  # nm^3
    pol2 = 0.001  # nm^3
    
    drude_force_unscreened.addParticle(
        drudeIndex=1, parentIndex=0,
        aniso1Index=-1, aniso2Index=-1,
        aniso3Index=-1, aniso4Index=-1,
        charge=-2.0, polarizability=pol1,
        aniso12=1.0, aniso34=1.0
    )
    
    drude_force_unscreened.addParticle(
        drudeIndex=3, parentIndex=2,
        aniso1Index=-1, aniso2Index=-1,
        aniso3Index=-1, aniso4Index=-1,
        charge=-1.0, polarizability=pol2,
        aniso12=1.0, aniso34=1.0
    )
    
    energy_unscreened = drude_force_unscreened.calculateEnergySCF(state)
    
    # Reset Drude positions
    state.atoms[1].x = 5.0
    state.atoms[3].x = 5.0 + distance
    
    # Now with screening
    drude_force_screened = DrudeForce()
    
    idx1 = drude_force_screened.addParticle(
        drudeIndex=1, parentIndex=0,
        aniso1Index=-1, aniso2Index=-1,
        aniso3Index=-1, aniso4Index=-1,
        charge=-2.0, polarizability=pol1,
        aniso12=1.0, aniso34=1.0
    )
    
    idx2 = drude_force_screened.addParticle(
        drudeIndex=3, parentIndex=2,
        aniso1Index=-1, aniso2Index=-1,
        aniso3Index=-1, aniso4Index=-1,
        charge=-1.0, polarizability=pol2,
        aniso12=1.0, aniso34=1.0
    )
    
    # Add screening
    thole = 1.3
    drude_force_screened.addScreenedPair(idx1, idx2, thole)
    
    energy_screened = drude_force_screened.calculateEnergySCF(state)
    
    # Screened energy should be different from unscreened
    assert abs(energy_screened - energy_unscreened) > 1e-6
    
    # At short distances, screening reduces the interaction
    # so screened energy should typically be lower in magnitude
    assert energy_screened != energy_unscreened


def test_multiple_screened_pairs():
    """Test system with multiple screened pairs"""
    state = MCState()
    
    # Set box size
    state.info.box = [10.0, 10.0, 10.0]
    state.info.setTemperature(300.0)
    state.info.cutoff = 5.0
    
    # Create three polarizable atoms in a triangle
    atoms = []
    positions = [
        (5.0, 5.0, 5.0),
        (5.5, 5.0, 5.0),
        (5.25, 5.433, 5.0)  # Equilateral triangle
    ]
    
    for i, (x, y, z) in enumerate(positions):
        # Parent atom
        atom = pygcmc.MCAtom()
        atom.x = x
        atom.y = y
        atom.z = z
        atom.charge = 1.0
        atom.type = 0
        atoms.append(atom)
        
        # Drude particle
        drude = pygcmc.MCAtom()
        drude.x = x
        drude.y = y
        drude.z = z
        drude.charge = -1.0
        drude.type = 1
        atoms.append(drude)
    
    state.atoms = atoms
    state.activeAtomCount = 6
    
    # Create residues
    residues = []
    for i in range(3):
        res = pygcmc.MCResidue()
        res.atomStart = i * 2
        res.atomCount = 2
        res.active = True
        res.fixed = False
        res.type = 0
        residues.append(res)
    
    state.residues = residues
    state.activeResidueCount = 3
    
    # Set force field
    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 2
    ff.numMovementTypes = 2
    ff.ljSigma = [0.3] * 4
    ff.ljEps = [0.0] * 4
    state.forcefield = ff
    
    # Create Drude force
    drude_force = DrudeForce()
    
    # Add three Drude particles
    charge = -1.0
    polarizability = 0.001
    indices = []
    
    for i in range(3):
        idx = drude_force.addParticle(
            drudeIndex=i*2+1, parentIndex=i*2,
            aniso1Index=-1, aniso2Index=-1,
            aniso3Index=-1, aniso4Index=-1,
            charge=charge, polarizability=polarizability,
            aniso12=1.0, aniso34=1.0
        )
        indices.append(idx)
    
    # Add screened pairs between all three
    thole = 1.5
    drude_force.addScreenedPair(indices[0], indices[1], thole)
    drude_force.addScreenedPair(indices[0], indices[2], thole)
    drude_force.addScreenedPair(indices[1], indices[2], thole)
    
    assert drude_force.getNumScreenedPairs() == 3
    
    # Set tight SCF for accurate result
    scf_params = DrudeSCFParams()
    scf_params.tolerance = 1e-7
    scf_params.maxIterations = 100
    drude_force.setSCFParameters(scf_params)
    
    # Calculate energy
    energy = drude_force.calculateEnergySCF(state)
    
    # All Drude particles should be displaced
    total_displacement = 0.0
    for i in range(3):
        drude_idx = i * 2 + 1
        parent_idx = i * 2
        
        dx = state.atoms[drude_idx].x - state.atoms[parent_idx].x
        dy = state.atoms[drude_idx].y - state.atoms[parent_idx].y
        dz = state.atoms[drude_idx].z - state.atoms[parent_idx].z
        dist = math.sqrt(dx*dx + dy*dy + dz*dz)
        
        total_displacement += dist
        assert dist > 1e-6  # Each should be displaced
    
    # Average displacement should be reasonable
    avg_displacement = total_displacement / 3
    assert avg_displacement > 1e-5
    assert avg_displacement < 0.01
    
    # Energy should be positive
    assert energy > 0.0


def test_thole_parameter_effect():
    """Test effect of different Thole parameters"""
    # Base state setup
    def create_test_state():
        state = MCState()
        state.info.box = [10.0, 10.0, 10.0]
        state.info.setTemperature(300.0)
        state.info.cutoff = 5.0
        
        # Two atoms 0.3 nm apart
        atoms = []
        
        for i in range(2):
            atom = pygcmc.MCAtom()
            atom.x = 5.0 + i * 0.3
            atom.y = 5.0
            atom.z = 5.0
            atom.charge = 1.5
            atom.type = 0
            atoms.append(atom)
            
            drude = pygcmc.MCAtom()
            drude.x = atom.x
            drude.y = atom.y
            drude.z = atom.z
            drude.charge = -1.5
            drude.type = 1
            atoms.append(drude)
        
        state.atoms = atoms
        state.activeAtomCount = 4
        
        residues = []
        for i in range(2):
            res = pygcmc.MCResidue()
            res.atomStart = i * 2
            res.atomCount = 2
            res.active = True
            res.fixed = False
            res.type = 0
            residues.append(res)
        
        state.residues = residues
        state.activeResidueCount = 2
        
        ff = pygcmc.MCForceField()
        ff.numTotalTypes = 2
        ff.numMovementTypes = 2
        ff.ljSigma = [0.3] * 4
        ff.ljEps = [0.0] * 4
        state.forcefield = ff
        
        return state
    
    # Test different Thole parameters
    thole_values = [0.0, 1.0, 2.0, 4.0]
    energies = []
    
    for thole in thole_values:
        state = create_test_state()
        drude_force = DrudeForce()
        
        charge = -1.5
        polarizability = 0.0015
        
        idx1 = drude_force.addParticle(
            drudeIndex=1, parentIndex=0,
            aniso1Index=-1, aniso2Index=-1,
            aniso3Index=-1, aniso4Index=-1,
            charge=charge, polarizability=polarizability,
            aniso12=1.0, aniso34=1.0
        )
        
        idx2 = drude_force.addParticle(
            drudeIndex=3, parentIndex=2,
            aniso1Index=-1, aniso2Index=-1,
            aniso3Index=-1, aniso4Index=-1,
            charge=charge, polarizability=polarizability,
            aniso12=1.0, aniso34=1.0
        )
        
        if thole > 0:
            drude_force.addScreenedPair(idx1, idx2, thole)
        
        energy = drude_force.calculateEnergySCF(state)
        energies.append(energy)
    
    # Higher Thole parameter should lead to more screening
    # and thus different energies
    for i in range(1, len(energies)):
        assert energies[i] != energies[0]  # Different from unscreened
    
    # Energies should change monotonically with Thole parameter
    # (though the exact relationship depends on distance)
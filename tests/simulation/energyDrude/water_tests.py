# tests/simulation/energyDrude/water_tests.py
"""SWM4-NDP water model tests for Drude oscillators."""

import pytest
import math
import pygcmc
from pygcmc import MCState, DrudeForce, DrudeSCFParams
from .helpers import (
    create_swm4_water,
    create_water_box,
    ONE_4PI_EPS0
)


def test_swm4_water_single():
    """Test single SWM4-NDP water molecule"""
    state = MCState()
    
    # Set box size
    state.info.box = [10.0, 10.0, 10.0]
    state.info.setTemperature(300.0)
    state.info.cutoff = 5.0
    
    # Create single water molecule
    water_atoms, water_res = create_swm4_water(5.0, 5.0, 5.0)
    state.atoms = water_atoms
    state.activeAtomCount = len(water_atoms)
    
    water_res.atomStart = 0
    state.residues = [water_res]
    state.activeResidueCount = 1
    
    # Set force field
    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 4  # O, D, H, M
    ff.numMovementTypes = 4
    ff.ljSigma = [0.318395] * 16
    ff.ljEps = [0.0] * 16
    ff.ljEps[0] = 0.21094 * 4.184  # O-O interaction
    state.forcefield = ff
    
    # Create Drude force
    drude_force = DrudeForce()
    
    # Add Drude particle for oxygen
    # SWM4-NDP parameters
    charge = -1.71636
    # Force constant from SWM4-NDP: 100000 kcal/mol/Å²
    k_kcal = 100000.0  # kcal/mol/Å²
    k_kj = k_kcal * 4.184 * 100.0  # Convert to kJ/mol/nm²
    polarizability = charge * charge / (ONE_4PI_EPS0 * k_kj)
    
    drude_idx = drude_force.addParticle(
        drudeIndex=1,  # Drude particle
        parentIndex=0,  # Oxygen
        aniso1Index=-1, aniso2Index=-1,
        aniso3Index=-1, aniso4Index=-1,
        charge=charge,
        polarizability=polarizability,
        aniso12=1.0, aniso34=1.0
    )
    
    # Calculate energy
    energy = drude_force.calculateEnergySCF(state)
    
    # For isolated water, Drude should be very close to oxygen
    dx = state.atoms[1].x - state.atoms[0].x
    dy = state.atoms[1].y - state.atoms[0].y
    dz = state.atoms[1].z - state.atoms[0].z
    dist = math.sqrt(dx*dx + dy*dy + dz*dz)
    
    # Should be very small after SCF
    assert dist < 1e-4


def test_swm4_water_dimer():
    """Test SWM4-NDP water dimer interaction"""
    state = MCState()
    
    # Set box size
    state.info.box = [10.0, 10.0, 10.0]
    state.info.setTemperature(300.0)
    state.info.cutoff = 5.0
    
    # Create two water molecules
    water1_atoms, water1_res = create_swm4_water(5.0, 5.0, 5.0)
    water2_atoms, water2_res = create_swm4_water(5.3, 5.0, 5.0)  # 3 Å apart
    
    # Combine atoms
    atoms = water1_atoms + water2_atoms
    water2_res.atomStart = len(water1_atoms)
    
    state.atoms = atoms
    state.activeAtomCount = len(atoms)
    state.residues = [water1_res, water2_res]
    state.activeResidueCount = 2
    
    # Set force field
    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 4
    ff.numMovementTypes = 4
    ff.ljSigma = [0.318395] * 16
    ff.ljEps = [0.0] * 16
    ff.ljEps[0] = 0.21094 * 4.184
    state.forcefield = ff
    
    # Create Drude force
    drude_force = DrudeForce()
    
    # Add Drude particles for both waters
    charge = -1.71636
    k_kj = 100000.0 * 4.184 * 100.0
    polarizability = charge * charge / (ONE_4PI_EPS0 * k_kj)
    
    # First water
    drude_force.addParticle(
        drudeIndex=1, parentIndex=0,
        aniso1Index=-1, aniso2Index=-1,
        aniso3Index=-1, aniso4Index=-1,
        charge=charge, polarizability=polarizability,
        aniso12=1.0, aniso34=1.0
    )
    
    # Second water
    drude_force.addParticle(
        drudeIndex=6, parentIndex=5,
        aniso1Index=-1, aniso2Index=-1,
        aniso3Index=-1, aniso4Index=-1,
        charge=charge, polarizability=polarizability,
        aniso12=1.0, aniso34=1.0
    )
    
    # Calculate energy
    energy = drude_force.calculateEnergySCF(state)
    
    # Both Drude particles should be displaced due to interaction
    # First water Drude
    dx1 = state.atoms[1].x - state.atoms[0].x
    dy1 = state.atoms[1].y - state.atoms[0].y
    dz1 = state.atoms[1].z - state.atoms[0].z
    dist1 = math.sqrt(dx1*dx1 + dy1*dy1 + dz1*dz1)
    
    # Second water Drude
    dx2 = state.atoms[6].x - state.atoms[5].x
    dy2 = state.atoms[6].y - state.atoms[5].y
    dz2 = state.atoms[6].z - state.atoms[5].z
    dist2 = math.sqrt(dx2*dx2 + dy2*dy2 + dz2*dz2)
    
    # Both should have some displacement due to intermolecular interactions
    assert dist1 > 1e-5
    assert dist2 > 1e-5
    
    # Energy should be non-zero
    assert energy > 0.0


def test_swm4_water_box():
    """Test small box of SWM4-NDP water molecules"""
    # Create box with 8 waters (2x2x2)
    state = create_water_box(8, box_size=2.0)
    
    # Create Drude force
    drude_force = DrudeForce()
    
    # Add Drude particles for all waters
    charge = -1.71636
    k_kj = 100000.0 * 4.184 * 100.0
    polarizability = charge * charge / (ONE_4PI_EPS0 * k_kj)
    
    for i in range(8):
        drude_idx = i * 5 + 1  # Drude is second atom in each water
        parent_idx = i * 5      # Oxygen is first atom
        
        drude_force.addParticle(
            drudeIndex=drude_idx,
            parentIndex=parent_idx,
            aniso1Index=-1, aniso2Index=-1,
            aniso3Index=-1, aniso4Index=-1,
            charge=charge,
            polarizability=polarizability,
            aniso12=1.0, aniso34=1.0
        )
    
    # Set tight SCF parameters
    scf_params = DrudeSCFParams()
    scf_params.tolerance = 1e-6
    scf_params.maxIterations = 100
    drude_force.setSCFParameters(scf_params)
    
    # Calculate energy
    energy = drude_force.calculateEnergySCF(state)
    
    # Check that all Drude particles found equilibrium
    total_displacement = 0.0
    for i in range(8):
        drude_idx = i * 5 + 1
        parent_idx = i * 5
        
        dx = state.atoms[drude_idx].x - state.atoms[parent_idx].x
        dy = state.atoms[drude_idx].y - state.atoms[parent_idx].y
        dz = state.atoms[drude_idx].z - state.atoms[parent_idx].z
        dist = math.sqrt(dx*dx + dy*dy + dz*dz)
        total_displacement += dist
    
    avg_displacement = total_displacement / 8
    
    # Average displacement should be small but non-zero
    assert avg_displacement > 1e-5
    assert avg_displacement < 0.01  # Less than 0.1 Å
    
    # Energy should be substantial due to many interactions
    assert energy > 0.0
    
    # Energy per water should be reasonable
    energy_per_water = energy / 8
    assert energy_per_water < 1000.0  # Reasonable upper bound


def test_swm4_water_polarization():
    """Test water polarization in external field"""
    state = MCState()
    
    # Set box size
    state.info.box = [10.0, 10.0, 10.0]
    state.info.setTemperature(300.0)
    state.info.cutoff = 5.0
    
    # Create single water molecule
    water_atoms, water_res = create_swm4_water(5.0, 5.0, 5.0)
    
    # Add external charges to create field
    # Positive charge on left, negative on right
    pos_charge = pygcmc.MCAtom()
    pos_charge.x = 3.0
    pos_charge.y = 5.0
    pos_charge.z = 5.0
    pos_charge.charge = 10.0
    pos_charge.type = 4
    
    neg_charge = pygcmc.MCAtom()
    neg_charge.x = 7.0
    neg_charge.y = 5.0
    neg_charge.z = 5.0
    neg_charge.charge = -10.0
    neg_charge.type = 4
    
    # Combine all atoms
    atoms = water_atoms + [pos_charge, neg_charge]
    water_res.atomStart = 0
    
    state.atoms = atoms
    state.activeAtomCount = len(atoms)
    state.residues = [water_res]
    state.activeResidueCount = 1
    
    # Set force field
    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 5  # O, D, H, M, external
    ff.numMovementTypes = 5
    ff.ljSigma = [0.3] * 25
    ff.ljEps = [0.0] * 25
    state.forcefield = ff
    
    # Create Drude force
    drude_force = DrudeForce()
    
    # Add Drude particle
    charge = -1.71636
    k_kj = 100000.0 * 4.184 * 100.0
    polarizability = charge * charge / (ONE_4PI_EPS0 * k_kj)
    
    drude_force.addParticle(
        drudeIndex=1, parentIndex=0,
        aniso1Index=-1, aniso2Index=-1,
        aniso3Index=-1, aniso4Index=-1,
        charge=charge, polarizability=polarizability,
        aniso12=1.0, aniso34=1.0
    )
    
    # Calculate energy
    energy = drude_force.calculateEnergySCF(state)
    
    # Drude should be displaced in field direction (negative x)
    dx = state.atoms[1].x - state.atoms[0].x
    
    # Should be displaced toward positive charge (negative x direction)
    assert dx < -0.001  # At least 0.001 nm displacement
    
    # Energy should be lowered by polarization
    assert energy > 0.0
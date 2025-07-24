"""
Tests comparing PyGCMC Drude implementation with OpenMM
Based on SWM4-NDP water model parameters
"""

import pytest
import numpy as np
import pygcmc
import math


# OpenMM SWM4-NDP parameters (from DRUDE_WATER_TESTS_OPENMM.md)
OPENMM_PARAMS = {
    'charges': {
        'O': 1.71636,
        'D': -1.71636,
        'H': 0.55733,
        'M': -1.11466
    },
    'polarizability': 0.978e-3,  # nm^3
    'thole': 1.3,
    'k_spring': 418400.0,  # kJ/mol/nm^2 (converted from 1000 kcal/mol/Å^2)
    'box_size': 1.86206,  # nm
    'cutoff': 0.9,  # nm
}


def create_swm4_water(origin=(0, 0, 0)):
    """Create SWM4-NDP water molecule at given origin"""
    x0, y0, z0 = origin
    
    atoms = []
    
    # Oxygen
    oxygen = pygcmc.MCAtom()
    oxygen.x, oxygen.y, oxygen.z = x0, y0, z0
    oxygen.charge = OPENMM_PARAMS['charges']['O']
    oxygen.type = 0
    atoms.append(oxygen)
    
    # Drude (initially at oxygen position)
    drude = pygcmc.MCAtom()
    drude.x, drude.y, drude.z = x0, y0, z0
    drude.charge = OPENMM_PARAMS['charges']['D']
    drude.type = 1
    atoms.append(drude)
    
    # Hydrogen 1 - along x-axis
    h1 = pygcmc.MCAtom()
    h1.x = x0 + 0.09572
    h1.y = y0
    h1.z = z0
    h1.charge = OPENMM_PARAMS['charges']['H']
    h1.type = 2
    atoms.append(h1)
    
    # Hydrogen 2 - 104.52 degree angle
    angle = 104.52 * 3.14159265359 / 180.0
    h2 = pygcmc.MCAtom()
    h2.x = x0 + 0.09572 * math.cos(angle)
    h2.y = y0 + 0.09572 * math.sin(angle) 
    h2.z = z0
    h2.charge = OPENMM_PARAMS['charges']['H']
    h2.type = 2
    atoms.append(h2)
    
    # Virtual site M - on bisector
    # Position is 0.024034 nm from O along bisector
    bisector_angle = angle / 2.0
    m_site = pygcmc.MCAtom()
    m_site.x = x0 + 0.024034 * math.cos(bisector_angle)
    m_site.y = y0 + 0.024034 * math.sin(bisector_angle)
    m_site.z = z0
    m_site.charge = OPENMM_PARAMS['charges']['M']
    m_site.type = 3
    atoms.append(m_site)
    
    return atoms


def create_water_system(num_waters=1, box_size=None):
    """Create a system with multiple water molecules"""
    if box_size is None:
        box_size = OPENMM_PARAMS['box_size']
    
    state = pygcmc.MCState()
    state.info.box = [box_size, box_size, box_size]
    state.info.cutoff = OPENMM_PARAMS['cutoff']
    
    # Force field
    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 2
    ff.numMovementTypes = 2
    ff.ljEps = [0.0, 0.0]
    ff.ljSigma = [0.1, 0.1]
    state.forcefield = ff
    
    all_atoms = []
    residues = []
    
    # Create water molecules in a grid
    spacing = box_size / (num_waters ** (1/3))
    positions = []
    for i in range(num_waters):
        x = (i % num_waters) * spacing
        y = ((i // num_waters) % num_waters) * spacing  
        z = (i // (num_waters * num_waters)) * spacing
        positions.append((x, y, z))
    
    for i, pos in enumerate(positions[:num_waters]):
        atoms = create_swm4_water(pos)
        all_atoms.extend(atoms)
        
        # Create residue for this water
        res = pygcmc.MCResidue()
        res.active = True
        res.fixed = False
        res.atomStart = i * 4
        res.atomCount = 4
        res.type = 0
        residues.append(res)
    
    state.atoms = all_atoms
    state.activeAtomCount = len(all_atoms)
    state.residues = residues
    state.activeResidueCount = len(residues)
    
    return state


def test_single_water_energy():
    """Test single water molecule energy (should be ~0 with only harmonic term)"""
    state = pygcmc.MCState()
    
    # Create single water
    state.atoms = create_swm4_water()
    state.activeAtomCount = 5
    state.info.box = [2.0, 2.0, 2.0]
    state.info.cutoff = OPENMM_PARAMS['cutoff']
    
    # Setup residue to mark all atoms as same molecule
    res = pygcmc.MCResidue()
    res.atomStart = 0
    res.atomCount = 5
    res.active = True
    res.type = 0
    state.residues = [res]
    state.activeResidueCount = 1
    
    # Setup Drude
    pygcmc.DrudeComplete.clear()
    
    particle = pygcmc.DrudeParticle()
    particle.drudeIndex = 1
    particle.parentIndex = 0
    particle.charge = OPENMM_PARAMS['charges']['D']
    particle.polarizability = OPENMM_PARAMS['polarizability']
    particle.computeSpringConstants()
    
    pygcmc.DrudeComplete.addParticle(particle)
    
    # Set SCF parameters matching OpenMM
    params = pygcmc.DrudeSCFParams()
    params.tolerance = 1e-5  # 0.01 kJ/mol/nm as in OpenMM
    params.maxIterations = 100
    params.maxDrudeDistance = 0.02
    params.enableHardWall = True  # Explicitly enable hard wall
    pygcmc.DrudeComplete.setParameters(params)
    
    # Calculate energy
    energy = pygcmc.DrudeComplete.calculateEnergy(state)
    
    # With Drude at parent position, energy should be very small
    # (only from intramolecular Coulomb if not properly excluded)
    assert abs(energy) < 1.0, f"Single water energy too large: {energy} kJ/mol"
    
    pygcmc.DrudeComplete.clear()


def test_water_dimer_interaction():
    """Test water dimer interaction energy"""
    state = pygcmc.MCState()
    
    # Create two waters separated by ~3 Å
    water1 = create_swm4_water((0, 0, 0))
    water2 = create_swm4_water((0.3, 0, 0))
    
    state.atoms = water1 + water2
    state.activeAtomCount = 10
    state.info.box = [2.0, 2.0, 2.0]
    state.info.cutoff = OPENMM_PARAMS['cutoff']
    
    # Setup residues for proper intramolecular exclusions
    res1 = pygcmc.MCResidue()
    res1.atomStart = 0
    res1.atomCount = 5
    res1.active = True
    res1.type = 0
    
    res2 = pygcmc.MCResidue()
    res2.atomStart = 5
    res2.atomCount = 5
    res2.active = True
    res2.type = 0
    
    state.residues = [res1, res2]
    state.activeResidueCount = 2
    
    # Setup Drude particles
    pygcmc.DrudeComplete.clear()
    
    # Water 1 Drude
    p1 = pygcmc.DrudeParticle()
    p1.drudeIndex = 1
    p1.parentIndex = 0
    p1.charge = OPENMM_PARAMS['charges']['D']
    p1.polarizability = OPENMM_PARAMS['polarizability']
    p1.computeSpringConstants()
    pygcmc.DrudeComplete.addParticle(p1)
    
    # Water 2 Drude
    p2 = pygcmc.DrudeParticle()
    p2.drudeIndex = 6  # Second water's Drude
    p2.parentIndex = 5  # Second water's oxygen
    p2.charge = OPENMM_PARAMS['charges']['D']
    p2.polarizability = OPENMM_PARAMS['polarizability']
    p2.computeSpringConstants()
    pygcmc.DrudeComplete.addParticle(p2)
    
    # Add Thole screening between waters
    pair = pygcmc.ScreenedPair()
    pair.dipole1 = 0
    pair.dipole2 = 1
    pair.thole = OPENMM_PARAMS['thole']
    pygcmc.DrudeComplete.addScreenedPair(pair)
    
    # Set SCF parameters
    params = pygcmc.DrudeSCFParams()
    params.tolerance = 1e-5
    params.maxIterations = 100
    params.maxDrudeDistance = 0.02
    params.enableHardWall = True  # Explicitly enable hard wall
    pygcmc.DrudeComplete.setParameters(params)
    
    # Calculate total energy
    total_energy = pygcmc.DrudeComplete.calculateEnergy(state)
    
    # Calculate individual water energies
    # Water 1 alone
    state1 = pygcmc.MCState()
    state1.atoms = water1
    state1.activeAtomCount = 5
    state1.info.box = state.info.box
    
    # Setup residue for water 1
    res1_only = pygcmc.MCResidue()
    res1_only.atomStart = 0
    res1_only.atomCount = 5
    res1_only.active = True
    res1_only.type = 0
    state1.residues = [res1_only]
    state1.activeResidueCount = 1
    
    pygcmc.DrudeComplete.clear()
    pygcmc.DrudeComplete.addParticle(p1)
    pygcmc.DrudeComplete.setParameters(params)
    energy1 = pygcmc.DrudeComplete.calculateEnergy(state1)
    
    # Interaction energy
    interaction_energy = total_energy - 2 * energy1  # Assuming both waters have same self-energy
    
    # Typical water dimer interaction is -20 to -25 kJ/mol
    assert -30 < interaction_energy < -10, \
        f"Water dimer interaction energy {interaction_energy} kJ/mol out of expected range"
    
    pygcmc.DrudeComplete.clear()



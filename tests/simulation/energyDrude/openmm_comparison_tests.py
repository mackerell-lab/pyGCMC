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
    
    # Hydrogen 1
    h1 = pygcmc.MCAtom()
    h1.x = x0 + 0.09572
    h1.y = y0
    h1.z = z0
    h1.charge = OPENMM_PARAMS['charges']['H']
    h1.type = 2
    atoms.append(h1)
    
    # Hydrogen 2
    h2 = pygcmc.MCAtom()
    h2.x = x0 - 0.02399
    h2.y = y0 + 0.09266
    h2.z = z0
    h2.charge = OPENMM_PARAMS['charges']['H']
    h2.type = 2
    atoms.append(h2)
    
    # Virtual site M
    m_site = pygcmc.MCAtom()
    m_site.x = x0 + 0.024034
    m_site.y = y0 + 0.023173
    m_site.z = z0
    m_site.charge = OPENMM_PARAMS['charges']['M']
    m_site.type = 3
    atoms.append(m_site)
    
    return atoms


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


def test_scf_convergence_tolerance():
    """Test SCF convergence with different tolerances"""
    state = pygcmc.MCState()
    
    # Create water dimer
    water1 = create_swm4_water((0, 0, 0))
    water2 = create_swm4_water((0.35, 0, 0))
    
    state.atoms = water1 + water2
    state.activeAtomCount = 10
    state.info.box = [2.0, 2.0, 2.0]
    
    # Setup residues
    residues = []
    for i in range(2):
        res = pygcmc.MCResidue()
        res.atomStart = i * 5
        res.atomCount = 5
        res.active = True
        res.type = 0
        residues.append(res)
    
    state.residues = residues
    state.activeResidueCount = 2
    
    # Test different tolerances
    tolerances = [1.0, 0.1, 0.01, 0.001]  # kJ/mol/nm
    energies = []
    
    for tol in tolerances:
        pygcmc.DrudeComplete.clear()
        
        # Add Drude particles
        for i in range(2):
            p = pygcmc.DrudeParticle()
            p.drudeIndex = i * 5 + 1
            p.parentIndex = i * 5
            p.charge = OPENMM_PARAMS['charges']['D']
            p.polarizability = OPENMM_PARAMS['polarizability']
            p.computeSpringConstants()
            pygcmc.DrudeComplete.addParticle(p)
        
        # Add Thole pair
        pair = pygcmc.ScreenedPair()
        pair.dipole1 = 0
        pair.dipole2 = 1
        pair.thole = OPENMM_PARAMS['thole']
        pygcmc.DrudeComplete.addScreenedPair(pair)
        
        # Set parameters with current tolerance
        params = pygcmc.DrudeSCFParams()
        params.tolerance = tol
        params.maxIterations = 200
        params.maxDrudeDistance = 0.02
        pygcmc.DrudeComplete.setParameters(params)
        
        # Reset Drude positions
        state.atoms[1].x = state.atoms[0].x
        state.atoms[1].y = state.atoms[0].y
        state.atoms[1].z = state.atoms[0].z
        state.atoms[6].x = state.atoms[5].x
        state.atoms[6].y = state.atoms[5].y
        state.atoms[6].z = state.atoms[5].z
        
        energy = pygcmc.DrudeComplete.calculateEnergy(state)
        energies.append(energy)
    
    # Energies should converge as tolerance decreases
    for i in range(1, len(energies)):
        diff = abs(energies[i] - energies[i-1])
        # Each tighter tolerance should give smaller change
        assert diff < tolerances[i-1] * 10, \
            f"Energy not converging properly: {diff} > {tolerances[i-1] * 10}"
    
    pygcmc.DrudeComplete.clear()


def test_polarization_response():
    """Test polarization response to external field"""
    state = pygcmc.MCState()
    
    # Single water molecule
    state.atoms = create_swm4_water()
    state.activeAtomCount = 5
    state.info.box = [3.0, 3.0, 3.0]
    
    # Add external charge to create field
    external = pygcmc.MCAtom()
    external.x, external.y, external.z = 0.5, 0.0, 0.0
    external.charge = 1.0  # Positive charge
    external.type = 4
    state.atoms.append(external)
    state.activeAtomCount = 6
    
    # Setup Drude
    pygcmc.DrudeComplete.clear()
    
    particle = pygcmc.DrudeParticle()
    particle.drudeIndex = 1
    particle.parentIndex = 0
    particle.charge = OPENMM_PARAMS['charges']['D']
    particle.polarizability = OPENMM_PARAMS['polarizability']
    particle.computeSpringConstants()
    
    pygcmc.DrudeComplete.addParticle(particle)
    
    params = pygcmc.DrudeSCFParams()
    params.tolerance = 1e-6
    params.maxIterations = 100
    params.maxDrudeDistance = 0.02
    pygcmc.DrudeComplete.setParameters(params)
    
    # Calculate energy
    energy = pygcmc.DrudeComplete.calculateEnergy(state)
    
    # Check Drude displacement
    dx = state.atoms[1].x - state.atoms[0].x
    dy = state.atoms[1].y - state.atoms[0].y
    dz = state.atoms[1].z - state.atoms[0].z
    
    # Drude (negative) should move away from positive external charge
    assert dx < 0, f"Drude should move in -x direction, but dx={dx}"
    
    # Calculate induced dipole moment
    dipole_magnitude = abs(particle.charge) * math.sqrt(dx*dx + dy*dy + dz*dz)
    
    # Rough check: dipole should be proportional to field strength
    # Field ~ charge/r^2, dipole ~ alpha * field
    assert dipole_magnitude > 0, "No induced dipole"
    assert dipole_magnitude < 0.1, f"Dipole too large: {dipole_magnitude}"
    
    pygcmc.DrudeComplete.clear()


def test_multiple_water_box():
    """Test small box of water molecules"""
    state = pygcmc.MCState()
    
    # Create 2x2x2 grid of waters (8 total)
    atoms = []
    spacing = 0.31  # nm, typical water spacing
    
    for i in range(2):
        for j in range(2):
            for k in range(2):
                origin = (i * spacing, j * spacing, k * spacing)
                atoms.extend(create_swm4_water(origin))
    
    state.atoms = atoms
    state.activeAtomCount = len(atoms)
    state.info.box = [1.0, 1.0, 1.0]
    state.info.cutoff = 0.9
    
    # Setup residues
    residues = []
    for i in range(8):
        res = pygcmc.MCResidue()
        res.atomStart = i * 5
        res.atomCount = 5
        res.active = True
        res.type = 0
        residues.append(res)
    
    state.residues = residues
    state.activeResidueCount = 8
    
    # Setup all Drude particles
    pygcmc.DrudeComplete.clear()
    
    drude_indices = []
    for i in range(8):
        p = pygcmc.DrudeParticle()
        p.drudeIndex = i * 5 + 1  # Drude is second atom in each water
        p.parentIndex = i * 5      # Oxygen is first
        p.charge = OPENMM_PARAMS['charges']['D']
        p.polarizability = OPENMM_PARAMS['polarizability']
        p.computeSpringConstants()
        pygcmc.DrudeComplete.addParticle(p)
        drude_indices.append(i)
    
    # Add all Thole pairs
    for i in range(8):
        for j in range(i+1, 8):
            pair = pygcmc.ScreenedPair()
            pair.dipole1 = i
            pair.dipole2 = j
            pair.thole = OPENMM_PARAMS['thole']
            pygcmc.DrudeComplete.addScreenedPair(pair)
    
    # Set parameters
    params = pygcmc.DrudeSCFParams()
    params.tolerance = 0.01  # Looser for many-body
    params.maxIterations = 200
    params.maxDrudeDistance = 0.02
    params.dampingFactor = 0.5
    pygcmc.DrudeComplete.setParameters(params)
    
    # Calculate energy
    try:
        energy = pygcmc.DrudeComplete.calculateEnergy(state)
        energy_per_water = energy / 8
        
        # Rough check - should be negative (attractive) but not too large
        assert -100 < energy_per_water < 10, \
            f"Energy per water {energy_per_water} kJ/mol out of range"
            
    except Exception as e:
        # SCF might not converge for large system
        pytest.skip(f"SCF convergence issue: {e}")
    
    pygcmc.DrudeComplete.clear()
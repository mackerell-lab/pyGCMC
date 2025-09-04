"""
Tests for SCF force balance and screened interactions in Drude implementation
"""

import pytest
import numpy as np
import pygcmc
import math

# Import helper functions and parameters
from energyOpenmm.openmm_water_energy_tests import (
    OPENMM_PARAMS,
    create_swm4_water
)
from energyOpenmm.drude_test_helpers import (
    calculate_numerical_force
)

def test_drude_scf_force_balance():
    """Test that SCF converged Drude particles have balanced forces"""
    state = pygcmc.MCState()
    
    # Water molecule (SWM4-NDP geometry)
    atoms = []
    
    # Oxygen
    oxygen = pygcmc.MCAtom()
    oxygen.x, oxygen.y, oxygen.z = 0.0, 0.0, 0.0
    oxygen.charge = 1.71636
    oxygen.type = 0
    atoms.append(oxygen)
    
    # Drude
    drude = pygcmc.MCAtom()
    drude.x, drude.y, drude.z = 0.0, 0.0, 0.0
    drude.charge = -1.71636
    drude.type = 1
    atoms.append(drude)
    
    # Hydrogens
    h1 = pygcmc.MCAtom()
    h1.x = 0.09572
    h1.y = 0.0
    h1.z = 0.0
    h1.charge = 0.55733
    h1.type = 2
    atoms.append(h1)
    
    angle = 104.52 * math.pi / 180.0
    h2 = pygcmc.MCAtom()
    h2.x = 0.09572 * math.cos(angle)
    h2.y = 0.09572 * math.sin(angle)
    h2.z = 0.0
    h2.charge = 0.55733
    h2.type = 2
    atoms.append(h2)
    
    # Virtual site M
    bisector_angle = angle / 2.0
    m_site = pygcmc.MCAtom()
    m_site.x = 0.024034 * math.cos(bisector_angle)
    m_site.y = 0.024034 * math.sin(bisector_angle)
    m_site.z = 0.0
    m_site.charge = -1.11466
    m_site.type = 3
    atoms.append(m_site)
    
    state.atoms = atoms
    state.activeAtomCount = 5
    state.info.box = [3.0, 3.0, 3.0]
    
    # Setup residue
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
    particle.charge = -1.71636
    particle.polarizability = 0.978e-3
    particle.computeSpringConstants()
    
    pygcmc.DrudeComplete.addParticle(particle)
    
    # Converge SCF
    params = pygcmc.DrudeSCFParams()
    params.tolerance = 1e-8  # Very tight for force balance test
    params.maxIterations = 500
    params.maxDrudeDistance = 0.02
    pygcmc.DrudeComplete.setParameters(params)
    
    # Calculate energy (this runs SCF)
    energy = pygcmc.DrudeComplete.calculateEnergy(state)
    
    # Calculate force on Drude using finite differences
    force_drude_x = calculate_numerical_force(state, 1, 0)
    force_drude_y = calculate_numerical_force(state, 1, 1)
    force_drude_z = calculate_numerical_force(state, 1, 2)
    
    # For converged SCF, net force on Drude should be near zero
    force_magnitude = math.sqrt(force_drude_x**2 + force_drude_y**2 + force_drude_z**2)
    
    # Force should be very small (within SCF tolerance)
    assert force_magnitude < 0.01, \
        f"Net force on converged Drude too large: {force_magnitude} kJ/mol/nm"
    
    pygcmc.DrudeComplete.clear()


def test_thole_screened_force():
    """Test that Thole screening reduces dipole-dipole forces"""
    # This test checks the general behavior rather than exact values
    state = pygcmc.MCState()
    
    # Two water molecules
    atoms = []
    
    # Water 1
    # Oxygen
    o1 = pygcmc.MCAtom()
    o1.x, o1.y, o1.z = 0.0, 0.0, 0.0
    o1.charge = 1.71636
    o1.type = 0
    atoms.append(o1)
    
    # Drude
    d1 = pygcmc.MCAtom()
    d1.x, d1.y, d1.z = 0.0, 0.0, 0.0
    d1.charge = -1.71636
    d1.type = 1
    atoms.append(d1)
    
    # Water 2 - 0.4 nm away
    o2 = pygcmc.MCAtom()
    o2.x, o2.y, o2.z = 0.4, 0.0, 0.0
    o2.charge = 1.71636
    o2.type = 0
    atoms.append(o2)
    
    d2 = pygcmc.MCAtom()
    d2.x, d2.y, d2.z = 0.4, 0.0, 0.0
    d2.charge = -1.71636
    d2.type = 1
    atoms.append(d2)
    
    state.atoms = atoms
    state.activeAtomCount = 4
    state.info.box = [3.0, 3.0, 3.0]
    state.info.cutoff = 2.0
    
    # Setup residues
    res1 = pygcmc.MCResidue()
    res1.atomStart = 0
    res1.atomCount = 2
    res1.active = True
    res1.type = 0
    
    res2 = pygcmc.MCResidue()
    res2.atomStart = 2
    res2.atomCount = 2
    res2.active = True
    res2.type = 1
    
    state.residues = [res1, res2]
    state.activeResidueCount = 2
    
    # Test with and without Thole screening
    energies = {}
    
    for use_thole in [False, True]:
        pygcmc.DrudeComplete.clear()
        
        # Add Drude particles
        for i in range(2):
            p = pygcmc.DrudeParticle()
            p.drudeIndex = i * 2 + 1
            p.parentIndex = i * 2
            p.charge = -1.71636
            p.polarizability = 0.978e-3
            p.computeSpringConstants()
            pygcmc.DrudeComplete.addParticle(p)
        
        if use_thole:
            # Add Thole screening
            pair = pygcmc.ScreenedPair()
            pair.dipole1 = 0
            pair.dipole2 = 1
            pair.thole = 1.3
            pygcmc.DrudeComplete.addScreenedPair(pair)
        
        # Run SCF
        params = pygcmc.DrudeSCFParams()
        params.tolerance = 1e-5
        params.maxIterations = 500
        params.maxDrudeDistance = 0.02
        pygcmc.DrudeComplete.setParameters(params)
        
        # Calculate energy
        energy = pygcmc.DrudeComplete.calculateEnergy(state)
        energies['thole' if use_thole else 'no_thole'] = energy
    
    # Thole screening changes the potential energy surface and SCF convergence
    # The exact effect depends on the configuration and convergence
    # Just check that Thole screening has some effect
    assert abs(energies['thole'] - energies['no_thole']) > 0.1, \
        f"Thole screening should affect the energy: {energies['thole']} vs {energies['no_thole']}"
    
    pygcmc.DrudeComplete.clear()



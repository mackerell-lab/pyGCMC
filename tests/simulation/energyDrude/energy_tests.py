"""
Energy calculation tests for Drude oscillators
"""

import pytest
import numpy as np
import pygcmc


def test_drude_with_mcstate():
    """Test Drude energy calculation with MCState"""
    # Create a simple MCState with two atoms (parent and Drude)
    state = pygcmc.MCState()
    
    # Parent atom (oxygen)
    parent = pygcmc.MCAtom()
    parent.x = 0.0
    parent.y = 0.0
    parent.z = 0.0
    parent.charge = 1.71636  # Positive charge to balance Drude
    parent.type = 0
    
    # Drude particle
    drude = pygcmc.MCAtom()
    drude.x = 0.0
    drude.y = 0.0
    drude.z = 0.01  # Slightly displaced
    drude.charge = -1.71636
    drude.type = 1
    
    state.atoms = [parent, drude]
    state.activeAtomCount = 2
    
    # Set box dimensions
    state.info.box = [3.0, 3.0, 3.0]
    state.info.cutoff = 1.0
    
    # Setup force field
    state.forcefield.numTotalTypes = 2
    state.forcefield.ljSigma = [0.0, 0.0]
    state.forcefield.ljEps = [0.0, 0.0]
    
    # Clear and setup Drude system
    pygcmc.DrudeComplete.clear()
    
    particle = pygcmc.DrudeParticle()
    particle.drudeIndex = 1
    particle.parentIndex = 0
    particle.charge = -1.71636
    particle.polarizability = 0.0009782237
    particle.computeSpringConstants()
    
    pygcmc.DrudeComplete.addParticle(particle)
    
    # Calculate energy
    energy = pygcmc.DrudeComplete.calculateEnergy(state)
    
    # Energy should be positive (spring is stretched)
    assert energy > 0
    
    # Clean up
    pygcmc.DrudeComplete.clear()


def test_harmonic_energy():
    """Test harmonic spring energy calculation"""
    state = pygcmc.MCState()
    
    # Create parent-Drude pair with known displacement
    parent = pygcmc.MCAtom()
    parent.x = 0.0
    parent.y = 0.0
    parent.z = 0.0
    parent.charge = 0.0
    parent.type = 0
    
    drude = pygcmc.MCAtom()
    displacement = 0.005  # 5 pm displacement
    drude.x = displacement
    drude.y = 0.0
    drude.z = 0.0
    drude.charge = -1.71636
    drude.type = 1
    
    state.atoms = [parent, drude]
    state.activeAtomCount = 2
    state.info.box = [3.0, 3.0, 3.0]
    
    # Setup Drude system
    pygcmc.DrudeComplete.clear()
    
    particle = pygcmc.DrudeParticle()
    particle.drudeIndex = 1
    particle.parentIndex = 0
    particle.charge = -1.71636
    particle.polarizability = 0.0009782237
    particle.computeSpringConstants()
    
    pygcmc.DrudeComplete.addParticle(particle)
    
    # Calculate energy
    energy = pygcmc.DrudeComplete.calculateEnergy(state)
    
    # Expected harmonic energy: 0.5 * k * r^2
    expected_energy = 0.5 * particle.kSpring * displacement * displacement
    
    # Should be close (within numerical precision)
    assert abs(energy - expected_energy) < 1e-6
    
    pygcmc.DrudeComplete.clear()


def test_zero_displacement_energy():
    """Test that energy is zero when Drude is at parent position"""
    state = pygcmc.MCState()
    
    # Create parent and Drude at same position
    parent = pygcmc.MCAtom()
    parent.x = 1.0
    parent.y = 1.0
    parent.z = 1.0
    parent.charge = 0.0
    parent.type = 0
    
    drude = pygcmc.MCAtom()
    drude.x = 1.0  # Same position as parent
    drude.y = 1.0
    drude.z = 1.0
    drude.charge = -1.71636
    drude.type = 1
    
    state.atoms = [parent, drude]
    state.activeAtomCount = 2
    state.info.box = [3.0, 3.0, 3.0]
    
    # Setup Drude system
    pygcmc.DrudeComplete.clear()
    
    particle = pygcmc.DrudeParticle()
    particle.drudeIndex = 1
    particle.parentIndex = 0
    particle.charge = -1.71636
    particle.polarizability = 0.0009782237
    particle.computeSpringConstants()
    
    pygcmc.DrudeComplete.addParticle(particle)
    
    # Calculate energy
    energy = pygcmc.DrudeComplete.calculateEnergy(state)
    
    # Energy should be zero (or very close due to numerical precision)
    assert abs(energy) < 1e-10
    
    pygcmc.DrudeComplete.clear()
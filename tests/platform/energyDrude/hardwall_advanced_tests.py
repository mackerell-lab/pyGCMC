"""
Advanced tests for hard wall constraints in Drude oscillators
"""

import pytest
import numpy as np
import pygcmc
import math


def test_hardwall_with_multiple_drudes():
    """Test hard wall with multiple Drude particles"""
    state = pygcmc.MCState()
    
    # Create two Drude oscillators
    atoms = []
    
    # First oscillator
    parent1 = pygcmc.MCAtom()
    parent1.x, parent1.y, parent1.z = 0.0, 0.0, 0.0
    parent1.charge = 0.0
    parent1.type = 0
    atoms.append(parent1)
    
    drude1 = pygcmc.MCAtom()
    drude1.x, drude1.y, drude1.z = 0.0, 0.0, 0.0
    drude1.charge = -1.0
    drude1.type = 1
    atoms.append(drude1)
    
    # Second oscillator
    parent2 = pygcmc.MCAtom()
    parent2.x, parent2.y, parent2.z = 1.0, 0.0, 0.0
    parent2.charge = 0.0
    parent2.type = 0
    atoms.append(parent2)
    
    drude2 = pygcmc.MCAtom()
    drude2.x, drude2.y, drude2.z = 1.0, 0.0, 0.0
    drude2.charge = -1.0
    drude2.type = 1
    atoms.append(drude2)
    
    # Add external charges to pull Drudes
    external1 = pygcmc.MCAtom()
    external1.x, external1.y, external1.z = -0.2, 0.0, 0.0
    external1.charge = 1.0
    external1.type = 2
    atoms.append(external1)
    
    external2 = pygcmc.MCAtom()
    external2.x, external2.y, external2.z = 1.2, 0.0, 0.0
    external2.charge = 1.0
    external2.type = 2
    atoms.append(external2)
    
    state.atoms = atoms
    state.activeAtomCount = 6
    state.info.box = [5.0, 5.0, 5.0]
    
    # Setup both Drude particles
    pygcmc.DrudeComplete.clear()
    
    for i in range(2):
        particle = pygcmc.DrudeParticle()
        particle.drudeIndex = i * 2 + 1
        particle.parentIndex = i * 2
        particle.charge = -1.0
        particle.polarizability = 0.001
        particle.computeSpringConstants()
        pygcmc.DrudeComplete.addParticle(particle)
    
    # Set SCF parameters with hard wall
    params = pygcmc.DrudeSCFParams()
    params.tolerance = 1e-6
    params.maxIterations = 100
    params.maxDrudeDistance = 0.02  # 0.2 Å
    params.enableHardWall = True
    pygcmc.DrudeComplete.setParameters(params)
    
    # Calculate energy
    energy = pygcmc.DrudeComplete.calculateEnergy(state)
    
    # Check both Drude particles are within hard wall
    for i in range(2):
        parent_idx = i * 2
        drude_idx = i * 2 + 1
        
        dx = state.atoms[drude_idx].x - state.atoms[parent_idx].x
        dy = state.atoms[drude_idx].y - state.atoms[parent_idx].y
        dz = state.atoms[drude_idx].z - state.atoms[parent_idx].z
        dist = math.sqrt(dx*dx + dy*dy + dz*dz)
        
        assert dist <= params.maxDrudeDistance * 1.01, \
            f"Drude {i} exceeds hard wall: {dist} > {params.maxDrudeDistance}"
    
    pygcmc.DrudeComplete.clear()


def test_hardwall_anisotropic_force():
    """Test hard wall with anisotropic external forces"""
    state = pygcmc.MCState()
    
    # Parent at origin
    parent = pygcmc.MCAtom()
    parent.x, parent.y, parent.z = 0.0, 0.0, 0.0
    parent.charge = 0.0
    parent.type = 0
    
    # Drude
    drude = pygcmc.MCAtom()
    drude.x, drude.y, drude.z = 0.0, 0.0, 0.0
    drude.charge = -1.0
    drude.type = 1
    
    # External charges in different directions
    external1 = pygcmc.MCAtom()
    external1.x, external1.y, external1.z = 0.1, 0.0, 0.0  # X direction
    external1.charge = 2.0
    external1.type = 2
    
    external2 = pygcmc.MCAtom()
    external2.x, external2.y, external2.z = 0.0, 0.1, 0.0  # Y direction
    external2.charge = 1.0
    external2.type = 2
    
    state.atoms = [parent, drude, external1, external2]
    state.activeAtomCount = 4
    state.info.box = [5.0, 5.0, 5.0]
    
    # Setup Drude
    pygcmc.DrudeComplete.clear()
    
    particle = pygcmc.DrudeParticle()
    particle.drudeIndex = 1
    particle.parentIndex = 0
    particle.charge = -1.0
    particle.polarizability = 0.001
    particle.computeSpringConstants()
    
    pygcmc.DrudeComplete.addParticle(particle)
    
    params = pygcmc.DrudeSCFParams()
    params.tolerance = 1e-6
    params.maxIterations = 100
    params.maxDrudeDistance = 0.015  # 0.15 Å
    params.enableHardWall = True  # Explicitly enable hard wall
    pygcmc.DrudeComplete.setParameters(params)
    
    energy = pygcmc.DrudeComplete.calculateEnergy(state)
    
    # Check final position
    dx = state.atoms[1].x - state.atoms[0].x
    dy = state.atoms[1].y - state.atoms[0].y
    dz = state.atoms[1].z - state.atoms[0].z
    distance = math.sqrt(dx*dx + dy*dy + dz*dz)
    
    assert distance <= params.maxDrudeDistance * 1.001, \
        f"Distance {distance} exceeds wall"
    
    # Drude should be pulled more in X direction (stronger force)
    assert abs(dx) > abs(dy), \
        "Drude should be displaced more in X direction"
    
    pygcmc.DrudeComplete.clear()
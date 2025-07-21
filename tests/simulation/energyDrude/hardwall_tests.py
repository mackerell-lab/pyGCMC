"""
Tests for hard wall constraints in Drude oscillators
Based on CHARMM implementation with default 0.02 nm (0.2 Å) constraint
"""

import pytest
import numpy as np
import pygcmc
import math


def test_hardwall_constraint_basic():
    """Test that Drude particles respect hard wall constraint"""
    state = pygcmc.MCState()
    
    # Create a parent-Drude pair with large initial displacement
    parent = pygcmc.MCAtom()
    parent.x, parent.y, parent.z = 0.0, 0.0, 0.0
    parent.charge = 0.0
    parent.type = 0
    
    drude = pygcmc.MCAtom()
    # Place Drude far from parent (0.1 nm = 1 Å)
    drude.x, drude.y, drude.z = 0.1, 0.0, 0.0
    drude.charge = -1.0
    drude.type = 1
    
    # Add external charge to pull Drude away
    external = pygcmc.MCAtom()
    external.x, external.y, external.z = 0.2, 0.0, 0.0
    external.charge = 1.0  # Positive to attract negative Drude
    external.type = 2
    
    state.atoms = [parent, drude, external]
    state.activeAtomCount = 3
    state.info.box = [5.0, 5.0, 5.0]
    
    # Setup Drude with SCF optimization
    pygcmc.DrudeComplete.clear()
    
    particle = pygcmc.DrudeParticle()
    particle.drudeIndex = 1
    particle.parentIndex = 0
    particle.charge = -1.0
    particle.polarizability = 0.001
    particle.computeSpringConstants()
    
    pygcmc.DrudeComplete.addParticle(particle)
    
    # Set SCF parameters with hard wall constraint
    params = pygcmc.DrudeSCFParams()
    params.tolerance = 1e-6
    params.maxIterations = 100
    params.maxDrudeDistance = 0.02  # Default hard wall at 0.2 Å
    pygcmc.DrudeComplete.setParameters(params)
    
    # Calculate energy (this will optimize Drude position)
    energy = pygcmc.DrudeComplete.calculateEnergy(state)
    
    # Check that Drude-parent distance is within hard wall
    dx = state.atoms[1].x - state.atoms[0].x
    dy = state.atoms[1].y - state.atoms[0].y
    dz = state.atoms[1].z - state.atoms[0].z
    distance = math.sqrt(dx*dx + dy*dy + dz*dz)
    
    assert distance <= params.maxDrudeDistance * 1.001, \
        f"Drude-parent distance {distance} exceeds hard wall {params.maxDrudeDistance}"
    
    pygcmc.DrudeComplete.clear()


def test_hardwall_different_limits():
    """Test hard wall constraint with different limits"""
    state = pygcmc.MCState()
    
    # Setup system with strong external field
    atoms = []
    
    parent = pygcmc.MCAtom()
    parent.x, parent.y, parent.z = 0.0, 0.0, 0.0
    parent.charge = 0.0
    parent.type = 0
    atoms.append(parent)
    
    drude = pygcmc.MCAtom()
    drude.x, drude.y, drude.z = 0.0, 0.0, 0.0
    drude.charge = -2.0  # Strong charge
    drude.type = 1
    atoms.append(drude)
    
    # Strong external charge
    external = pygcmc.MCAtom()
    external.x, external.y, external.z = 0.1, 0.0, 0.0
    external.charge = 10.0  # Very strong
    external.type = 2
    atoms.append(external)
    
    state.atoms = atoms
    state.activeAtomCount = 3
    state.info.box = [5.0, 5.0, 5.0]
    
    # Test different hard wall distances
    wall_distances = [0.01, 0.02, 0.03, 0.05]  # nm
    measured_distances = []
    
    for max_dist in wall_distances:
        pygcmc.DrudeComplete.clear()
        
        particle = pygcmc.DrudeParticle()
        particle.drudeIndex = 1
        particle.parentIndex = 0
        particle.charge = -2.0
        particle.polarizability = 0.001
        particle.computeSpringConstants()
        pygcmc.DrudeComplete.addParticle(particle)
        
        params = pygcmc.DrudeSCFParams()
        params.tolerance = 1e-6
        params.maxIterations = 100
        params.maxDrudeDistance = max_dist
        pygcmc.DrudeComplete.setParameters(params)
        
        # Reset Drude position
        state.atoms[1].x = 0.0
        state.atoms[1].y = 0.0
        state.atoms[1].z = 0.0
        
        energy = pygcmc.DrudeComplete.calculateEnergy(state)
        
        # Measure actual distance
        dx = state.atoms[1].x - state.atoms[0].x
        dy = state.atoms[1].y - state.atoms[0].y
        dz = state.atoms[1].z - state.atoms[0].z
        distance = math.sqrt(dx*dx + dy*dy + dz*dz)
        measured_distances.append(distance)
        
        # Should not exceed hard wall
        assert distance <= max_dist * 1.001, \
            f"Distance {distance} exceeds wall {max_dist}"
    
    # Distances should increase with wall distance
    for i in range(1, len(measured_distances)):
        assert measured_distances[i] >= measured_distances[i-1], \
            "Larger wall should allow larger distance"
    
    pygcmc.DrudeComplete.clear()


def test_hardwall_energy_consistency():
    """Test that hard wall doesn't cause energy discontinuities"""
    state = pygcmc.MCState()
    
    # Simple system
    parent = pygcmc.MCAtom()
    parent.x, parent.y, parent.z = 0.0, 0.0, 0.0
    parent.charge = 0.0
    parent.type = 0
    
    drude = pygcmc.MCAtom()
    drude.x, drude.y, drude.z = 0.0, 0.0, 0.0
    drude.charge = -1.0
    drude.type = 1
    
    # Moderate external charge - use larger distance to create weaker field
    # With new spring constant (100x smaller), we need much weaker field
    external = pygcmc.MCAtom()
    external.x, external.y, external.z = 5.0, 0.0, 0.0  # 5 nm away for gentle field
    external.charge = 1.0
    external.type = 2
    
    state.atoms = [parent, drude, external]
    state.activeAtomCount = 3
    state.info.box = [5.0, 5.0, 5.0]
    
    # Test energies with different wall distances
    wall_distances = np.linspace(0.01, 0.05, 10)
    energies = []
    
    for max_dist in wall_distances:
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
        params.maxDrudeDistance = max_dist
        pygcmc.DrudeComplete.setParameters(params)
        
        # Reset Drude
        state.atoms[1].x = 0.0
        state.atoms[1].y = 0.0
        state.atoms[1].z = 0.0
        
        energy = pygcmc.DrudeComplete.calculateEnergy(state)
        energies.append(energy)
    
    # Check for smooth energy variation (no jumps)
    for i in range(1, len(energies)):
        energy_change = abs(energies[i] - energies[i-1])
        # Energy should change smoothly
        assert energy_change < 50.0, \
            f"Energy jump too large: {energies[i-1]} -> {energies[i]}"
    
    pygcmc.DrudeComplete.clear()


def test_hardwall_with_multiple_drudes():
    """Test hard wall with multiple Drude particles"""
    state = pygcmc.MCState()
    
    # Create two Drude oscillators
    atoms = []
    
    # First oscillator
    p1 = pygcmc.MCAtom()
    p1.x, p1.y, p1.z = 0.0, 0.0, 0.0
    p1.charge = 0.0
    p1.type = 0
    atoms.append(p1)
    
    d1 = pygcmc.MCAtom()
    d1.x, d1.y, d1.z = 0.0, 0.0, 0.0
    d1.charge = -1.0
    d1.type = 1
    atoms.append(d1)
    
    # Second oscillator
    p2 = pygcmc.MCAtom()
    p2.x, p2.y, p2.z = 0.3, 0.0, 0.0
    p2.charge = 0.0
    p2.type = 0
    atoms.append(p2)
    
    d2 = pygcmc.MCAtom()
    d2.x, d2.y, d2.z = 0.3, 0.0, 0.0
    d2.charge = -1.0
    d2.type = 1
    atoms.append(d2)
    
    state.atoms = atoms
    state.activeAtomCount = 4
    state.info.box = [5.0, 5.0, 5.0]
    
    # Setup Drude system
    pygcmc.DrudeComplete.clear()
    
    for i in range(2):
        particle = pygcmc.DrudeParticle()
        particle.drudeIndex = 2*i + 1
        particle.parentIndex = 2*i
        particle.charge = -1.0
        particle.polarizability = 0.001
        particle.computeSpringConstants()
        pygcmc.DrudeComplete.addParticle(particle)
    
    # Add interaction between Drudes
    pair = pygcmc.ScreenedPair()
    pair.dipole1 = 0
    pair.dipole2 = 1
    pair.thole = 1.3
    pygcmc.DrudeComplete.addScreenedPair(pair)
    
    # Set parameters with hard wall
    params = pygcmc.DrudeSCFParams()
    params.tolerance = 1e-6
    params.maxIterations = 100
    params.maxDrudeDistance = 0.02
    pygcmc.DrudeComplete.setParameters(params)
    
    # Calculate energy
    energy = pygcmc.DrudeComplete.calculateEnergy(state)
    
    # Check both Drudes respect hard wall
    for i in range(2):
        dx = state.atoms[2*i+1].x - state.atoms[2*i].x
        dy = state.atoms[2*i+1].y - state.atoms[2*i].y
        dz = state.atoms[2*i+1].z - state.atoms[2*i].z
        distance = math.sqrt(dx*dx + dy*dy + dz*dz)
        
        assert distance <= params.maxDrudeDistance * 1.001, \
            f"Drude {i} distance {distance} exceeds wall"
    
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
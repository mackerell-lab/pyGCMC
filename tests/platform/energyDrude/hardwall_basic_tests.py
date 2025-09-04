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
    params.enableHardWall = True  # Explicitly enable hard wall
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
        params.enableHardWall = True  # Explicitly enable hard wall
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

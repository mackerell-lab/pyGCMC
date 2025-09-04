#!/usr/bin/env python3
"""
Test hard wall constraint on/off functionality
"""

import pytest
import numpy as np
import pygcmc

def test_hard_wall_default_off():
    """Test that hard wall is OFF by default"""
    params = pygcmc.DrudeSCFParams()
    assert params.enableHardWall == False  # Should be False by default
    assert params.maxDrudeDistance == 0.02  # Default 0.2 Å

def test_hard_wall_can_be_enabled():
    """Test that hard wall can be enabled"""
    params = pygcmc.DrudeSCFParams()
    params.enableHardWall = True
    assert params.enableHardWall == True

def check_hard_wall_off_behavior():
    """Check SCF without hard wall allows larger displacements (helper function)"""
    # This is a helper function, not a test
    # Create a simple system with strong external field
    state = pygcmc.MCState()
    
    # Parent atom
    parent = pygcmc.MCAtom()
    parent.x, parent.y, parent.z = 0.0, 0.0, 0.0
    parent.charge = 0.0
    parent.type = 0
    
    # Drude atom
    drude = pygcmc.MCAtom()
    drude.x, drude.y, drude.z = 0.0, 0.0, 0.0
    drude.charge = -1.0
    drude.type = 1
    
    # External charge to create strong field
    external = pygcmc.MCAtom()
    external.x, external.y, external.z = 0.5, 0.0, 0.0  # 5 Å away
    external.charge = 10.0  # Strong positive charge
    external.type = 2
    
    state.atoms = [parent, drude, external]
    state.activeAtomCount = 3
    state.info.box = [3.0, 3.0, 3.0]
    
    # Setup Drude system
    pygcmc.DrudeComplete.clear()
    
    particle = pygcmc.DrudeParticle()
    particle.drudeIndex = 1
    particle.parentIndex = 0
    particle.charge = -1.0
    particle.polarizability = 0.001  # Larger polarizability for bigger displacement
    particle.aniso1Index = -1
    particle.aniso2Index = -1
    particle.aniso3Index = -1
    particle.aniso4Index = -1
    particle.aniso12 = 1.0
    particle.aniso34 = 1.0
    particle.computeSpringConstants()
    
    pygcmc.DrudeComplete.addParticle(particle)
    
    # Set SCF parameters WITHOUT hard wall
    params = pygcmc.DrudeSCFParams()
    params.tolerance = 0.01
    params.maxIterations = 100
    params.dampingFactor = 0.5
    params.maxDrudeDistance = 0.02  # 0.2 Å
    params.enableHardWall = False  # Disable hard wall
    pygcmc.DrudeComplete.setParameters(params)
    
    # Calculate energy (SCF will optimize Drude position)
    energy = pygcmc.DrudeComplete.calculateEnergy(state)
    
    # Check Drude displacement
    dx = state.atoms[1].x - state.atoms[0].x
    dy = state.atoms[1].y - state.atoms[0].y
    dz = state.atoms[1].z - state.atoms[0].z
    displacement = np.sqrt(dx*dx + dy*dy + dz*dz)
    
    print(f"Without hard wall: displacement = {displacement:.4f} nm")
    
    # Without hard wall, displacement can exceed maxDrudeDistance
    # Strong field should pull Drude beyond 0.02 nm
    assert displacement > 0.02  # Should exceed the "max" distance
    
    pygcmc.DrudeComplete.clear()
    return displacement

def check_hard_wall_on_behavior():
    """Check SCF with hard wall limits displacements (helper function)"""
    # This is a helper function, not a test
    # Same system as above
    state = pygcmc.MCState()
    
    parent = pygcmc.MCAtom()
    parent.x, parent.y, parent.z = 0.0, 0.0, 0.0
    parent.charge = 0.0
    parent.type = 0
    
    drude = pygcmc.MCAtom()
    drude.x, drude.y, drude.z = 0.0, 0.0, 0.0
    drude.charge = -1.0
    drude.type = 1
    
    external = pygcmc.MCAtom()
    external.x, external.y, external.z = 0.5, 0.0, 0.0
    external.charge = 10.0
    external.type = 2
    
    state.atoms = [parent, drude, external]
    state.activeAtomCount = 3
    state.info.box = [3.0, 3.0, 3.0]
    
    # Setup Drude system
    pygcmc.DrudeComplete.clear()
    
    particle = pygcmc.DrudeParticle()
    particle.drudeIndex = 1
    particle.parentIndex = 0
    particle.charge = -1.0
    particle.polarizability = 0.001
    particle.aniso1Index = -1
    particle.aniso2Index = -1
    particle.aniso3Index = -1
    particle.aniso4Index = -1
    particle.aniso12 = 1.0
    particle.aniso34 = 1.0
    particle.computeSpringConstants()
    
    pygcmc.DrudeComplete.addParticle(particle)
    
    # Set SCF parameters WITH hard wall
    params = pygcmc.DrudeSCFParams()
    params.tolerance = 0.01
    params.maxIterations = 100
    params.dampingFactor = 0.5
    params.maxDrudeDistance = 0.02  # 0.2 Å
    params.enableHardWall = True  # Enable hard wall
    pygcmc.DrudeComplete.setParameters(params)
    
    # Calculate energy
    energy = pygcmc.DrudeComplete.calculateEnergy(state)
    
    # Check Drude displacement
    dx = state.atoms[1].x - state.atoms[0].x
    dy = state.atoms[1].y - state.atoms[0].y
    dz = state.atoms[1].z - state.atoms[0].z
    displacement = np.sqrt(dx*dx + dy*dy + dz*dz)
    
    print(f"With hard wall: displacement = {displacement:.4f} nm")
    
    # With hard wall, displacement should be limited
    assert displacement <= 0.02 * 1.01  # Allow 1% tolerance for numerical precision
    
    pygcmc.DrudeComplete.clear()
    return displacement

def test_hard_wall_comparison():
    """Compare displacements with and without hard wall"""
    disp_without_wall = check_hard_wall_off_behavior()
    disp_with_wall = check_hard_wall_on_behavior()
    
    # Without wall should have larger displacement
    assert disp_without_wall > disp_with_wall
    print(f"\nDisplacement ratio: {disp_without_wall/disp_with_wall:.2f}x")

def test_hard_wall_off_behavior():
    """Test SCF without hard wall allows larger displacements"""
    displacement = check_hard_wall_off_behavior()
    assert displacement > 0.02  # Already checked in helper

def test_hard_wall_on_behavior():
    """Test SCF with hard wall limits displacements"""
    displacement = check_hard_wall_on_behavior()
    assert displacement <= 0.02 * 1.01  # Already checked in helper

if __name__ == "__main__":
    test_hard_wall_default_off()
    print("✓ Default hard wall setting: OFF")
    
    test_hard_wall_can_be_enabled()
    print("✓ Hard wall can be enabled")
    
    test_hard_wall_comparison()
    print("✓ Hard wall constraint working correctly")
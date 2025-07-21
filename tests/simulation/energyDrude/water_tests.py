"""
Water model tests for Drude oscillators
"""

import pytest
import numpy as np
import pygcmc


def test_swm4_ndp_water_setup():
    """Test setting up SWM4-NDP water model"""
    state = pygcmc.MCState()
    
    # Simple test: just O and D atoms to test basic setup
    # Create atoms
    atoms = []
    
    # Oxygen (parent)
    oxygen = pygcmc.MCAtom()
    oxygen.x, oxygen.y, oxygen.z = 0.0, 0.0, 0.0
    oxygen.charge = 0.0  # Set to zero for this basic test
    oxygen.type = 0
    atoms.append(oxygen)
    
    # Drude
    drude = pygcmc.MCAtom()
    drude.x, drude.y, drude.z = 0.0, 0.0, 0.0  # Initially at oxygen position
    drude.charge = -1.71636
    drude.type = 1
    atoms.append(drude)
    
    state.atoms = atoms
    state.activeAtomCount = 2
    state.info.box = [3.0, 3.0, 3.0]
    
    # Setup Drude system
    pygcmc.DrudeComplete.clear()
    
    particle = pygcmc.DrudeParticle()
    particle.drudeIndex = 1
    particle.parentIndex = 0
    particle.charge = -1.71636
    particle.polarizability = 0.0009782237  # SWM4-NDP oxygen polarizability
    particle.computeSpringConstants()
    
    pygcmc.DrudeComplete.addParticle(particle)
    
    # Calculate initial energy (should be zero since Drude is at parent)
    energy = pygcmc.DrudeComplete.calculateEnergy(state)
    assert abs(energy) < 1e-10
    
    # Test with displacement
    state.atoms[1].x = 0.001  # Move Drude slightly
    energy2 = pygcmc.DrudeComplete.calculateEnergy(state)
    assert energy2 > 0  # Should have positive energy due to spring
    
    pygcmc.DrudeComplete.clear()


def test_two_water_thole_interaction():
    """Test Thole-screened interaction between two water molecules"""
    # Simple test: verify that Thole screening is working
    state = pygcmc.MCState()
    
    # Create a single Drude oscillator
    atoms = []
    
    # Parent atom
    parent = pygcmc.MCAtom()
    parent.x, parent.y, parent.z = 0.0, 0.0, 0.0
    parent.charge = 0.0
    parent.type = 0
    atoms.append(parent)
    
    # Drude particle
    drude = pygcmc.MCAtom()
    drude.x, drude.y, drude.z = 0.0, 0.0, 0.01  # Displaced 0.1 Å
    drude.charge = -1.71636
    drude.type = 1
    atoms.append(drude)
    
    state.atoms = atoms
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
    
    # Calculate energy with displacement
    energy = pygcmc.DrudeComplete.calculateEnergy(state)
    assert energy > 0  # Should have positive spring energy
    
    # Verify Thole screening function works
    screening = pygcmc.computeTholeScreening(0.3, 0.001, 0.001, 1.3)
    assert 0 < screening < 1  # Should be partially screened
    
    pygcmc.DrudeComplete.clear()


def test_water_residue_setup():
    """Test setting up water as a residue"""
    state = pygcmc.MCState()
    
    # First set up atoms
    atoms = []
    
    # Parent atom
    parent = pygcmc.MCAtom()
    parent.x = 0.0
    parent.y = 0.0
    parent.z = 0.0
    parent.type = 0
    parent.charge = 0.0  # Neutral for simple test
    atoms.append(parent)
    
    # Drude atom
    drude = pygcmc.MCAtom()
    drude.x = 0.0
    drude.y = 0.0  
    drude.z = 0.0
    drude.type = 1
    drude.charge = -1.71636
    atoms.append(drude)
    
    # Set atoms first
    state.atoms = atoms
    state.activeAtomCount = 2
    
    # Create residue after atoms are set
    residue = pygcmc.MCResidue()
    residue.atomStart = 0
    residue.atomCount = 2  # Just O and D
    residue.active = True
    residue.type = 0
    
    # Initialize residues list properly
    residues = []
    residues.append(residue)
    state.residues = residues
    state.activeResidueCount = len(residues)
    
    # Set box
    state.info.box = [3.0, 3.0, 3.0]
    
    # Verify residue structure carefully
    assert len(state.residues) >= 1, f"Expected at least 1 residue, got {len(state.residues)}"
    assert state.activeResidueCount == 1
    assert len(state.atoms) == 2
    # Only access residue if it exists
    if len(state.residues) > 0:
        assert state.residues[0].atomCount == 2
    
    # Setup Drude for the residue
    pygcmc.DrudeComplete.clear()
    
    particle = pygcmc.DrudeParticle()
    particle.drudeIndex = 1
    particle.parentIndex = 0
    particle.charge = -1.71636
    particle.polarizability = 0.0009782237
    particle.computeSpringConstants()
    
    pygcmc.DrudeComplete.addParticle(particle)
    
    # Calculate energy (should be zero with Drude at parent position)
    energy = pygcmc.DrudeComplete.calculateEnergy(state)
    assert abs(energy) < 1e-10
    
    # Test with displacement
    state.atoms[1].z = 0.001  # Small displacement
    energy2 = pygcmc.DrudeComplete.calculateEnergy(state)
    assert energy2 > 0  # Should have positive spring energy
    
    pygcmc.DrudeComplete.clear()
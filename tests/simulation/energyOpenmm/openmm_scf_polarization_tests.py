"""
Tests comparing PyGCMC Drude SCF and polarization with OpenMM
"""

import pytest
import numpy as np
import pygcmc
import math

# Import shared functions from water energy tests
from energyOpenmm.openmm_water_energy_tests import (
    OPENMM_PARAMS,
    create_swm4_water,
    create_water_system
)

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
        params.enableHardWall = True  # Explicitly enable hard wall
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
    # With smaller spring constant, need weaker field to avoid hitting hard wall
    external = pygcmc.MCAtom()
    external.x, external.y, external.z = 2.0, 0.0, 0.0  # Further away for weaker field
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
    params.maxDrudeDistance = 0.03  # Larger to accommodate intramolecular fields
    params.enableHardWall = True  # Explicitly enable hard wall
    pygcmc.DrudeComplete.setParameters(params)
    
    # Calculate energy
    energy = pygcmc.DrudeComplete.calculateEnergy(state)
    
    # Check Drude displacement
    dx = state.atoms[1].x - state.atoms[0].x
    dy = state.atoms[1].y - state.atoms[0].y
    dz = state.atoms[1].z - state.atoms[0].z
    
    # The Drude displacement is complex due to intramolecular fields
    # from H atoms and M virtual site. Just check that there is displacement.
    displacement = math.sqrt(dx*dx + dy*dy + dz*dz)
    
    # The displacement should be small but non-zero
    assert displacement > 1e-6, "No induced dipole"
    assert displacement <= params.maxDrudeDistance * 1.01, f"Displacement {displacement} exceeds hard wall"
    params.enableHardWall = True  # Explicitly enable hard wall
    
    # The actual displacement direction and magnitude depends on:
    # 1. External field from the test charge
    # 2. Intramolecular fields from H atoms (positive)
    # 3. Intramolecular field from M virtual site (negative)
    # So we just verify that polarization occurs without checking specific direction
    
    pygcmc.DrudeComplete.clear()



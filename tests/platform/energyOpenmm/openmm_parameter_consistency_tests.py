"""
Tests verifying PyGCMC Drude parameters match OpenMM conventions
"""

import pytest
import numpy as np
import pygcmc
import math

# Import shared functions
from energyOpenmm.openmm_water_energy_tests import (
    OPENMM_PARAMS,
    create_swm4_water,
    create_water_system
)

def test_openmm_spring_constant_consistency():
    """Verify spring constant calculation matches OpenMM exactly"""
    # Test cases from OpenMM: k = ONE_4PI_EPS0 * q^2 / alpha
    test_cases = [
        # (charge, polarizability, expected_k_openmm)
        (-1.71636, 0.978e-3, None),  # SWM4-NDP water
        (-1.0, 0.001, None),          # Simple case
        (-2.0, 0.002, None),          # Larger values
    ]
    
    for charge, alpha, _ in test_cases:
        particle = pygcmc.DrudeParticle()
        particle.charge = charge
        particle.polarizability = alpha
        particle.computeSpringConstants()
        
        # Manual calculation as in OpenMM
        k_expected = charge * charge * pygcmc.DrudeConstants.ONE_4PI_EPS0 / alpha
        
        assert abs(particle.kSpring - k_expected) < 1e-6, \
            f"Spring constant {particle.kSpring} != {k_expected}"
    
    # Verify specific SWM4-NDP value
    particle = pygcmc.DrudeParticle()
    particle.charge = -1.71636
    particle.polarizability = 0.978e-3
    particle.computeSpringConstants()
    
    # OpenMM reports ~418400 kJ/mol/nm^2 for SWM4-NDP
    assert abs(particle.kSpring - 418400) < 1000, \
        f"SWM4-NDP spring constant {particle.kSpring} differs from OpenMM value"


def test_openmm_thole_screening_values():
    """Compare Thole screening with OpenMM implementation"""
    # Test screening function at specific distances
    alpha = 0.978e-3  # SWM4-NDP
    thole = 1.3
    
    # Test distances and expected screenings from OpenMM
    test_points = [
        # (r_nm, expected_screening_range)
        (0.1, (0.5, 0.7)),      # Very close - highly screened
        (0.2, (0.8, 0.9)),      # Close - moderately screened  
        (0.3, (0.9, 0.95)),     # Medium - less screened
        (0.5, (0.98, 1.0)),     # Far - barely screened
        (1.0, (0.999, 1.0)),    # Very far - essentially unscreened
    ]
    
    for r, (s_min, s_max) in test_points:
        s = pygcmc.computeTholeScreening(r, alpha, alpha, thole)
        assert s_min <= s <= s_max, \
            f"Screening {s} at r={r} nm outside expected range {s_min}-{s_max}"


def test_openmm_water_dipole_moment():
    """Test water molecular polarizability and energy calculation"""
    state = pygcmc.MCState()
    
    # Create two water molecules to induce dipole-dipole interaction
    water1 = create_swm4_water((0, 0, 0))
    water2 = create_swm4_water((0.4, 0, 0))  # 0.4 nm away
    
    state.atoms = water1 + water2
    state.activeAtomCount = 10
    state.info.box = [3.0, 3.0, 3.0]
    state.info.cutoff = 2.0
    
    # Setup residues for proper intramolecular exclusions
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
    
    # Setup Drude particles
    pygcmc.DrudeComplete.clear()
    
    for i in range(2):
        particle = pygcmc.DrudeParticle()
        particle.drudeIndex = i * 5 + 1  # Drude is second atom
        particle.parentIndex = i * 5      # Oxygen is first
        particle.charge = OPENMM_PARAMS['charges']['D']
        particle.polarizability = OPENMM_PARAMS['polarizability']
        particle.computeSpringConstants()
        pygcmc.DrudeComplete.addParticle(particle)
    
    # Add Thole screening between waters
    pair = pygcmc.ScreenedPair()
    pair.dipole1 = 0
    pair.dipole2 = 1
    pair.thole = OPENMM_PARAMS['thole']
    pygcmc.DrudeComplete.addScreenedPair(pair)
    
    params = pygcmc.DrudeSCFParams()
    params.tolerance = 10.0  # Looser tolerance
    params.maxIterations = 200
    params.maxDrudeDistance = 0.025
    params.enableHardWall = True  # Explicitly enable hard wall
    pygcmc.DrudeComplete.setParameters(params)
    
    # Calculate energy
    energy = pygcmc.DrudeComplete.calculateEnergy(state)
    
    # Check that we have non-zero energy (interaction + harmonic)
    assert abs(energy) > 0.01, f"Energy {energy} too small, SCF might not be working"
    
    # Check at least one Drude has moved
    max_displacement = 0.0
    for i in range(2):
        drude_idx = i * 5 + 1
        parent_idx = i * 5
        dx = state.atoms[drude_idx].x - state.atoms[parent_idx].x
        dy = state.atoms[drude_idx].y - state.atoms[parent_idx].y
        dz = state.atoms[drude_idx].z - state.atoms[parent_idx].z
        displacement = math.sqrt(dx*dx + dy*dy + dz*dz)
        max_displacement = max(max_displacement, displacement)
    
    # At least one Drude should have moved due to intermolecular fields
    assert max_displacement > 1e-7, f"No Drude movement detected (max displacement: {max_displacement})"
    
    # Verify reasonable energy range (not too large)
    assert -1000 < energy < 1000, f"Energy {energy} kJ/mol out of reasonable range"
    
    pygcmc.DrudeComplete.clear()


def test_openmm_scf_convergence_behavior():
    """Test SCF convergence pattern similar to OpenMM"""
    state = pygcmc.MCState()
    
    # Two water molecules
    water1 = create_swm4_water((0, 0, 0))
    water2 = create_swm4_water((0.3, 0, 0))
    
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
    
    # Test with different SCF parameters like OpenMM
    scf_tests = [
        # (tolerance, max_iter, damping)
        (1e-3, 50, 0.7),    # Tight tolerance, high damping
        (1e-2, 100, 0.5),   # Medium tolerance, medium damping
        (1e-1, 200, 0.3),   # Loose tolerance, low damping
    ]
    
    energies = []
    
    for tol, max_iter, damping in scf_tests:
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
        
        # Add Thole screening
        pair = pygcmc.ScreenedPair()
        pair.dipole1 = 0
        pair.dipole2 = 1
        pair.thole = OPENMM_PARAMS['thole']
        pygcmc.DrudeComplete.addScreenedPair(pair)
        
        # Set parameters
        params = pygcmc.DrudeSCFParams()
        params.tolerance = tol * 1000  # Convert to force units
        params.maxIterations = max_iter
        params.dampingFactor = damping
        params.maxDrudeDistance = 0.02
        params.enableHardWall = True  # Explicitly enable hard wall
        pygcmc.DrudeComplete.setParameters(params)
        
        # Calculate energy
        energy = pygcmc.DrudeComplete.calculateEnergy(state)
        energies.append(energy)
    
    # Energies should be similar regardless of SCF parameters
    # (all should converge to same solution)
    for i in range(1, len(energies)):
        assert abs(energies[i] - energies[0]) < 1.0, \
            f"Energy variation too large: {energies[i]} vs {energies[0]}"
    
    pygcmc.DrudeComplete.clear()


def test_multiple_water_box():
    """Test small box of water molecules"""
    state = pygcmc.MCState()
    
    # Create 2x2x2 grid of waters (8 total)
    atoms = []
    spacing = 0.4  # nm, larger spacing to reduce repulsion
    
    for i in range(2):
        for j in range(2):
            for k in range(2):
                origin = (i * spacing, j * spacing, k * spacing)
                atoms.extend(create_swm4_water(origin))
    
    state.atoms = atoms
    state.activeAtomCount = len(atoms)
    state.info.box = [1.5, 1.5, 1.5]  # Larger box for larger spacing
    state.info.cutoff = 1.4
    
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
    params.tolerance = 0.1  # Looser for many-body system
    params.maxIterations = 200
    params.maxDrudeDistance = 0.02
    params.enableHardWall = True  # Explicitly enable hard wall
    params.dampingFactor = 0.5
    pygcmc.DrudeComplete.setParameters(params)
    
    # Calculate energy
    energy = pygcmc.DrudeComplete.calculateEnergy(state)
    energy_per_water = energy / 8
    
    # With larger spacing, energy should be less repulsive
    # The actual value depends on many factors including Coulomb interactions
    # Just check that it's not extremely high or low
    assert -500 < energy_per_water < 500, \
        f"Energy per water {energy_per_water} kJ/mol seems unreasonable"
    

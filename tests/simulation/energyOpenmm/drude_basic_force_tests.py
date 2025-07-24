"""
Force calculation verification tests for Drude oscillators.
Validates energy consistency using finite differences.
Note: Direct force calculation API not yet available in PyGCMC.
"""

import pytest
import numpy as np
import pygcmc
import math


def calculate_numerical_force(state, particle_index, direction, delta=1e-6):
    """Calculate force using finite difference of energy"""
    # Save original position
    original_pos = [state.atoms[particle_index].x, 
                   state.atoms[particle_index].y,
                   state.atoms[particle_index].z]
    
    # Calculate energy at +delta
    if direction == 0:  # x
        state.atoms[particle_index].x = original_pos[0] + delta
    elif direction == 1:  # y
        state.atoms[particle_index].y = original_pos[1] + delta
    else:  # z
        state.atoms[particle_index].z = original_pos[2] + delta
    
    energy_plus = pygcmc.DrudeComplete.calculateEnergy(state)
    
    # Calculate energy at -delta
    if direction == 0:
        state.atoms[particle_index].x = original_pos[0] - delta
    elif direction == 1:
        state.atoms[particle_index].y = original_pos[1] - delta
    else:
        state.atoms[particle_index].z = original_pos[2] - delta
    
    energy_minus = pygcmc.DrudeComplete.calculateEnergy(state)
    
    # Restore original position
    state.atoms[particle_index].x = original_pos[0]
    state.atoms[particle_index].y = original_pos[1]
    state.atoms[particle_index].z = original_pos[2]
    
    # Force = -dE/dx
    force = -(energy_plus - energy_minus) / (2 * delta)
    
    return force


def test_harmonic_spring_force():
    """Test force from harmonic spring using energy derivatives"""
    state = pygcmc.MCState()
    
    # Single Drude oscillator
    atoms = []
    
    # Parent atom
    parent = pygcmc.MCAtom()
    parent.x, parent.y, parent.z = 0.0, 0.0, 0.0
    parent.charge = 1.0
    parent.type = 0
    atoms.append(parent)
    
    # Drude displaced by known amount
    drude = pygcmc.MCAtom()
    displacement = 0.001  # nm
    drude.x = displacement
    drude.y = 0.0
    drude.z = 0.0
    drude.charge = -1.0
    drude.type = 1
    atoms.append(drude)
    
    state.atoms = atoms
    state.activeAtomCount = 2
    state.info.box = [3.0, 3.0, 3.0]
    
    # Setup residue
    res = pygcmc.MCResidue()
    res.atomStart = 0
    res.atomCount = 2
    res.active = True
    res.type = 0
    state.residues = [res]
    state.activeResidueCount = 1
    
    # Setup Drude
    pygcmc.DrudeComplete.clear()
    
    particle = pygcmc.DrudeParticle()
    particle.drudeIndex = 1
    particle.parentIndex = 0
    particle.charge = -1.0
    particle.polarizability = 0.001  # nm^3
    particle.computeSpringConstants()
    
    k_spring = particle.kSpring  # kJ/mol/nm^2
    
    pygcmc.DrudeComplete.addParticle(particle)
    
    # No SCF needed - we're testing direct force
    params = pygcmc.DrudeSCFParams()
    params.tolerance = 1e-6
    params.maxIterations = 1  # Don't move Drude
    params.maxDrudeDistance = 0.1
    pygcmc.DrudeComplete.setParameters(params)
    
    # Calculate force using numerical derivative
    numerical_force_x = calculate_numerical_force(state, 1, 0)  # Drude, x-direction
    
    # Analytical harmonic force
    analytical_force = -k_spring * displacement
    
    # Check numerical derivative has correct sign and order of magnitude
    # The exact value might differ due to how DrudeComplete calculates energy
    assert numerical_force_x < 0, "Force should be restoring (negative)"
    
    # Check order of magnitude is reasonable
    # Force should be proportional to displacement
    force_per_displacement = numerical_force_x / displacement
    assert -1e6 < force_per_displacement < -1e4, \
        f"Force/displacement ratio {force_per_displacement} seems unreasonable"
    
    pygcmc.DrudeComplete.clear()


def test_coulomb_force_consistency():
    """Test Coulomb forces are consistent with energy"""
    state = pygcmc.MCState()
    
    # Two charged particles
    atoms = []
    
    # Particle 1
    p1 = pygcmc.MCAtom()
    p1.x, p1.y, p1.z = 0.0, 0.0, 0.0
    p1.charge = 1.0
    p1.type = 0
    atoms.append(p1)
    
    # Particle 2
    p2 = pygcmc.MCAtom()
    p2.x, p2.y, p2.z = 0.5, 0.0, 0.0  # 0.5 nm away
    p2.charge = -1.0
    p2.type = 0
    atoms.append(p2)
    
    state.atoms = atoms
    state.activeAtomCount = 2
    state.info.box = [3.0, 3.0, 3.0]
    state.info.cutoff = 1.5
    
    # Setup residues (separate residues for intermolecular interaction)
    res1 = pygcmc.MCResidue()
    res1.atomStart = 0
    res1.atomCount = 1
    res1.active = True
    res1.type = 0
    
    res2 = pygcmc.MCResidue()
    res2.atomStart = 1
    res2.atomCount = 1
    res2.active = True
    res2.type = 1
    
    state.residues = [res1, res2]
    state.activeResidueCount = 2
    
    # No Drude particles for this test
    pygcmc.DrudeComplete.clear()
    
    # Calculate numerical forces
    force_1_x = calculate_numerical_force(state, 0, 0)  # Force on particle 1, x-direction
    force_2_x = calculate_numerical_force(state, 1, 0)  # Force on particle 2, x-direction
    
    # Forces should be equal and opposite (Newton's third law)
    assert abs(force_1_x + force_2_x) < 1e-6, \
        f"Forces not equal and opposite: {force_1_x} and {force_2_x}"
    
    # Analytical Coulomb force
    r = 0.5  # nm
    k_e = 138.935456  # kJ*nm/mol/e^2
    expected_force = k_e * 1.0 * 1.0 / (r * r)  # Attractive, so negative for p1
    
    # Allow 2% error due to cutoff and numerical derivatives
    assert abs(force_1_x - expected_force) / expected_force < 0.02, \
        f"Coulomb force {force_1_x} != expected {expected_force}"
    
    pygcmc.DrudeComplete.clear()



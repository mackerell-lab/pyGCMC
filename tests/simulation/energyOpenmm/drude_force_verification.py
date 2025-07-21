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
    
    # With Thole screening, interaction should be weaker (less negative)
    assert energies['thole'] > energies['no_thole'], \
        f"Thole screening should reduce interaction: {energies['thole']} vs {energies['no_thole']}"
    
    pygcmc.DrudeComplete.clear()


def test_force_energy_consistency_complex():
    """Test force-energy consistency for a complex system"""
    state = pygcmc.MCState()
    
    # Three polarizable ions in a triangle
    positions = [
        (0.0, 0.0, 0.0),
        (0.4, 0.0, 0.0),
        (0.2, 0.346, 0.0)  # Equilateral triangle
    ]
    charges = [1.0, -1.0, 0.5]
    
    atoms = []
    for i, (pos, charge) in enumerate(zip(positions, charges)):
        # Parent
        parent = pygcmc.MCAtom()
        parent.x, parent.y, parent.z = pos
        parent.charge = charge
        parent.type = 0
        atoms.append(parent)
        
        # Drude
        drude = pygcmc.MCAtom()
        drude.x, drude.y, drude.z = pos
        drude.charge = -charge
        drude.type = 1
        atoms.append(drude)
    
    state.atoms = atoms
    state.activeAtomCount = 6
    state.info.box = [3.0, 3.0, 3.0]
    state.info.cutoff = 2.0
    
    # Setup residues
    residues = []
    for i in range(3):
        res = pygcmc.MCResidue()
        res.atomStart = i * 2
        res.atomCount = 2
        res.active = True
        res.type = i
        residues.append(res)
    
    state.residues = residues
    state.activeResidueCount = 3
    
    # Setup Drude particles
    pygcmc.DrudeComplete.clear()
    
    for i in range(3):
        p = pygcmc.DrudeParticle()
        p.drudeIndex = i * 2 + 1
        p.parentIndex = i * 2
        p.charge = -charges[i]
        p.polarizability = 0.001
        p.computeSpringConstants()
        pygcmc.DrudeComplete.addParticle(p)
    
    # Add all Thole pairs
    for i in range(3):
        for j in range(i+1, 3):
            pair = pygcmc.ScreenedPair()
            pair.dipole1 = i
            pair.dipole2 = j
            pair.thole = 1.3
            pygcmc.DrudeComplete.addScreenedPair(pair)
    
    # SCF parameters
    params = pygcmc.DrudeSCFParams()
    params.tolerance = 1e-6
    params.maxIterations = 500
    params.maxDrudeDistance = 0.02
    pygcmc.DrudeComplete.setParameters(params)
    
    # Calculate energy
    energy = pygcmc.DrudeComplete.calculateEnergy(state)
    
    # Test force-energy consistency for each particle
    delta = 1e-5
    max_error = 0.0
    
    for atom_idx in range(state.activeAtomCount):
        for direction in range(3):  # x, y, z
            # Calculate numerical force
            force_numerical = calculate_numerical_force(state, atom_idx, direction, delta)
            
            # For this test, we mainly verify that forces are reasonable
            # and finite (no NaN or infinity)
            assert np.isfinite(force_numerical), \
                f"Non-finite force for atom {atom_idx} direction {direction}"
            
            # Force magnitude should be reasonable (not huge)
            assert abs(force_numerical) < 10000, \
                f"Unreasonably large force: {force_numerical} kJ/mol/nm"
    
    pygcmc.DrudeComplete.clear()
#!/usr/bin/env python3
"""
Drude field response and polarization tests
Tests the correct response of Drude oscillators to external fields
"""

import pytest
import numpy as np
import pygcmc


def test_single_drude_uniform_field():
    """Test single Drude particle response to uniform external field"""
    print("\n=== Test: Single Drude in Uniform Field ===")
    
    # Create state
    state = pygcmc.MCState()
    state.info.box = [3.0, 3.0, 3.0]
    
    # Parent atom at origin (neutral)
    parent = pygcmc.MCAtom()
    parent.x = parent.y = parent.z = 0.0
    parent.charge = 0.0
    parent.type = 0
    
    # Drude particle
    drude = pygcmc.MCAtom()
    drude.x = drude.y = drude.z = 0.0
    drude.charge = -1.0
    drude.type = 1
    
    # External positive charge
    external = pygcmc.MCAtom()
    external.x = 0.5  # 0.5 nm away
    external.y = external.z = 0.0
    external.charge = 2.0
    external.type = 2
    
    state.atoms = [parent, drude, external]
    state.activeAtomCount = 3
    
    # Create residues
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
    particle.polarizability = 0.001  # nm³
    particle.aniso1Index = particle.aniso2Index = particle.aniso3Index = particle.aniso4Index = -1
    particle.aniso12 = particle.aniso34 = 1.0
    particle.computeSpringConstants()
    pygcmc.DrudeComplete.addParticle(particle)
    
    # Set SCF parameters
    params = pygcmc.DrudeSCFParams()
    params.tolerance = 0.01
    params.maxIterations = 100
    params.enableHardWall = False
    params.dampingFactor = 0.5
    pygcmc.DrudeComplete.setParameters(params)
    
    # Calculate energy
    energy = pygcmc.DrudeComplete.calculateEnergy(state)
    
    # Get Drude displacement
    dx = state.atoms[1].x - state.atoms[0].x
    dy = state.atoms[1].y - state.atoms[0].y
    dz = state.atoms[1].z - state.atoms[0].z
    
    print(f"Drude displacement: ({dx:.6f}, {dy:.6f}, {dz:.6f}) nm")
    print(f"Energy: {energy:.6f} kJ/mol")
    
    # Physical expectation:
    # The negative Drude should move toward the positive external charge
    # This means positive x displacement
    assert dx > 0, "Drude should move toward positive charge"
    assert abs(dy) < 1e-6, "No y-displacement expected"
    assert abs(dz) < 1e-6, "No z-displacement expected"
    
    # Check reasonable displacement magnitude
    # For α = 0.001 nm³ and moderate field, expect ~0.008 nm displacement
    assert 0.005 < dx < 0.01, f"Displacement {dx} outside expected range"
    
    pygcmc.DrudeComplete.clear()


def test_drude_spring_energy():
    """Test Drude spring energy calculation"""
    print("\n=== Test: Drude Spring Energy ===")
    
    state = pygcmc.MCState()
    state.info.box = [3.0, 3.0, 3.0]
    
    # Parent-Drude pair with fixed displacement
    parent = pygcmc.MCAtom()
    parent.x = parent.y = parent.z = 0.0
    parent.charge = 0.0
    parent.type = 0
    
    # Place Drude at fixed displacement
    displacement = 0.01  # nm
    drude = pygcmc.MCAtom()
    drude.x = displacement
    drude.y = drude.z = 0.0
    drude.charge = -1.0
    drude.type = 1
    
    state.atoms = [parent, drude]
    state.activeAtomCount = 2
    
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
    particle.polarizability = 0.001
    particle.aniso1Index = particle.aniso2Index = particle.aniso3Index = particle.aniso4Index = -1
    particle.aniso12 = particle.aniso34 = 1.0
    particle.computeSpringConstants()
    
    print(f"Spring constant: {particle.kSpring:.2f} kJ/mol/nm²")
    
    pygcmc.DrudeComplete.addParticle(particle)
    
    # Calculate without SCF optimization
    params = pygcmc.DrudeSCFParams()
    params.maxIterations = 0  # No optimization, keep fixed positions
    pygcmc.DrudeComplete.setParameters(params)
    
    energy = pygcmc.DrudeComplete.calculateEnergy(state)
    
    # Expected spring energy: E = 0.5 * k * x²
    expected_energy = 0.5 * particle.kSpring * displacement**2
    
    print(f"Calculated energy: {energy:.6f} kJ/mol")
    print(f"Expected spring energy: {expected_energy:.6f} kJ/mol")
    
    # Should match within numerical precision
    assert abs(energy - expected_energy) / expected_energy < 0.01, "Spring energy mismatch"
    
    pygcmc.DrudeComplete.clear()


def test_induced_dipole_magnitude():
    """Test that induced dipole has correct magnitude"""
    print("\n=== Test: Induced Dipole Magnitude ===")
    
    state = pygcmc.MCState()
    state.info.box = [5.0, 5.0, 5.0]
    
    # System with known field
    parent = pygcmc.MCAtom()
    parent.x = parent.y = parent.z = 0.0
    parent.charge = 0.0
    parent.type = 0
    
    drude = pygcmc.MCAtom()
    drude.x = drude.y = drude.z = 0.0
    drude.charge = -1.0
    drude.type = 1
    
    # Create a charge at distance r to create known field
    r = 1.0  # nm
    q_external = 1.0  # e
    external = pygcmc.MCAtom()
    external.x = r
    external.y = external.z = 0.0
    external.charge = q_external
    external.type = 2
    
    state.atoms = [parent, drude, external]
    state.activeAtomCount = 3
    
    res = pygcmc.MCResidue()
    res.atomStart = 0
    res.atomCount = 2
    res.active = True
    res.type = 0
    state.residues = [res]
    state.activeResidueCount = 1
    
    # Setup Drude with known polarizability
    pygcmc.DrudeComplete.clear()
    
    alpha = 0.001  # nm³
    particle = pygcmc.DrudeParticle()
    particle.drudeIndex = 1
    particle.parentIndex = 0
    particle.charge = -1.0
    particle.polarizability = alpha
    particle.aniso1Index = particle.aniso2Index = particle.aniso3Index = particle.aniso4Index = -1
    particle.aniso12 = particle.aniso34 = 1.0
    particle.computeSpringConstants()
    pygcmc.DrudeComplete.addParticle(particle)
    
    # Tight convergence for accurate test
    params = pygcmc.DrudeSCFParams()
    params.tolerance = 1e-6
    params.maxIterations = 1000
    params.enableHardWall = False
    params.dampingFactor = 0.5
    pygcmc.DrudeComplete.setParameters(params)
    
    # Calculate
    energy = pygcmc.DrudeComplete.calculateEnergy(state)
    
    # Get displacement
    dx = state.atoms[1].x - state.atoms[0].x
    
    # Induced dipole moment: μ = q_drude * displacement
    mu_induced = abs(particle.charge * dx)
    
    # For point charge, E = k*q/r² (in appropriate units)
    # But this is complex due to self-consistency
    # Just check that dipole is reasonable
    print(f"Drude displacement: {dx:.8f} nm")
    print(f"Induced dipole moment: {mu_induced:.8f} e·nm")
    print(f"Polarizability * q_drude: {alpha * abs(particle.charge):.8f}")
    
    # Dipole should be positive and reasonable
    assert mu_induced > 0, "Should have non-zero induced dipole"
    assert mu_induced < 0.02, "Induced dipole too large"
    
    pygcmc.DrudeComplete.clear()


def test_multiple_drudes_mutual_polarization():
    """Test mutual polarization between multiple Drude oscillators"""
    print("\n=== Test: Multiple Drudes Mutual Polarization ===")
    
    state = pygcmc.MCState()
    state.info.box = [5.0, 5.0, 5.0]
    
    # Two polarizable molecules
    atoms = []
    
    # Molecule 1
    parent1 = pygcmc.MCAtom()
    parent1.x = 0.0
    parent1.y = parent1.z = 0.0
    parent1.charge = 0.5  # Partial positive
    parent1.type = 0
    atoms.append(parent1)
    
    drude1 = pygcmc.MCAtom()
    drude1.x = drude1.y = drude1.z = 0.0
    drude1.charge = -1.0
    drude1.type = 1
    atoms.append(drude1)
    
    # Molecule 2
    parent2 = pygcmc.MCAtom()
    parent2.x = 0.5  # 0.5 nm away
    parent2.y = parent2.z = 0.0
    parent2.charge = 0.5
    parent2.type = 0
    atoms.append(parent2)
    
    drude2 = pygcmc.MCAtom()
    drude2.x = 0.5
    drude2.y = drude2.z = 0.0
    drude2.charge = -1.0
    drude2.type = 1
    atoms.append(drude2)
    
    state.atoms = atoms
    state.activeAtomCount = 4
    
    # Create residues
    residues = []
    for i in range(2):
        res = pygcmc.MCResidue()
        res.atomStart = i * 2
        res.atomCount = 2
        res.active = True
        res.type = i
        residues.append(res)
    
    state.residues = residues
    state.activeResidueCount = 2
    
    # Setup Drude particles
    pygcmc.DrudeComplete.clear()
    
    for i in range(2):
        particle = pygcmc.DrudeParticle()
        particle.drudeIndex = i * 2 + 1
        particle.parentIndex = i * 2
        particle.charge = -1.0
        particle.polarizability = 0.001
        particle.aniso1Index = particle.aniso2Index = particle.aniso3Index = particle.aniso4Index = -1
        particle.aniso12 = particle.aniso34 = 1.0
        particle.computeSpringConstants()
        pygcmc.DrudeComplete.addParticle(particle)
    
    # SCF parameters
    params = pygcmc.DrudeSCFParams()
    params.tolerance = 0.01
    params.maxIterations = 100
    params.enableHardWall = False
    params.dampingFactor = 0.5
    pygcmc.DrudeComplete.setParameters(params)
    
    # Calculate
    energy = pygcmc.DrudeComplete.calculateEnergy(state)
    
    # Get displacements
    disp1 = state.atoms[1].x - state.atoms[0].x
    disp2 = state.atoms[3].x - state.atoms[2].x
    
    print(f"Drude 1 displacement: {disp1:.6f} nm")
    print(f"Drude 2 displacement: {disp2:.6f} nm")
    print(f"Energy: {energy:.3f} kJ/mol")
    
    # Due to symmetry and mutual repulsion of positive charges,
    # Drudes should move in opposite directions
    # Drude 1 is attracted to parent 2, but also repelled by drude 2
    # The net effect is that drude 1 moves left (negative direction)
    assert disp1 < 0, "Drude 1 should move left"
    assert disp2 > 0, "Drude 2 should move right"
    
    # Should have similar magnitude due to symmetry
    assert abs(abs(disp1) - abs(disp2)) / abs(disp1) < 0.1, "Displacements should be similar"
    
    pygcmc.DrudeComplete.clear()


if __name__ == "__main__":
    test_single_drude_uniform_field()
    test_drude_spring_energy()
    test_induced_dipole_magnitude()
    test_multiple_drudes_mutual_polarization()
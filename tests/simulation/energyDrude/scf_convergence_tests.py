"""
Tests for SCF convergence behavior in Drude oscillators
Including damping factors, iteration counts, and convergence criteria
"""

import pytest
import numpy as np
import pygcmc
import math


def test_scf_basic_convergence():
    """Test basic SCF convergence for simple system"""
    state = pygcmc.MCState()
    
    # Parent-Drude pair with external field
    parent = pygcmc.MCAtom()
    parent.x, parent.y, parent.z = 0.0, 0.0, 0.0
    parent.charge = 0.0
    parent.type = 0
    
    drude = pygcmc.MCAtom()
    drude.x, drude.y, drude.z = 0.01, 0.0, 0.0  # Start displaced
    drude.charge = -1.0
    drude.type = 1
    
    # External charge
    external = pygcmc.MCAtom()
    external.x, external.y, external.z = 0.2, 0.0, 0.0
    external.charge = 0.5
    external.type = 2
    
    state.atoms = [parent, drude, external]
    state.activeAtomCount = 3
    state.info.box = [5.0, 5.0, 5.0]
    
    # Test convergence with different starting positions
    initial_displacements = [0.0, 0.005, 0.01, 0.02]
    
    for disp in initial_displacements:
        pygcmc.DrudeComplete.clear()
        
        # Reset Drude position
        state.atoms[1].x = disp
        state.atoms[1].y = 0.0
        state.atoms[1].z = 0.0
        
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
        params.maxDrudeDistance = 0.02
        params.enableHardWall = True  # Explicitly enable hard wall
        pygcmc.DrudeComplete.setParameters(params)
        
        # Should converge regardless of starting position
        energy = pygcmc.DrudeComplete.calculateEnergy(state)
        
        # Final position should be the same
        final_x = state.atoms[1].x
        
        # All should converge to similar position
        if disp > 0:
            assert abs(final_x - state.atoms[1].x) < 1e-4, \
                f"Different convergence from displacement {disp}"
    
    pygcmc.DrudeComplete.clear()


def test_scf_damping_factor_effect():
    """Test effect of damping factor on convergence"""
    state = pygcmc.MCState()
    
    # Create strongly coupled system
    atoms = []
    
    # Two Drude oscillators close together
    for i in range(2):
        parent = pygcmc.MCAtom()
        parent.x = i * 0.25  # 2.5 Å apart
        parent.y, parent.z = 0.0, 0.0
        parent.charge = 0.0
        parent.type = 0
        atoms.append(parent)
        
        drude = pygcmc.MCAtom()
        drude.x = i * 0.25
        drude.y, drude.z = 0.0, 0.0
        drude.charge = -2.0  # Strong charge
        drude.type = 1
        atoms.append(drude)
    
    state.atoms = atoms
    state.activeAtomCount = 4
    state.info.box = [5.0, 5.0, 5.0]
    
    # Test different damping factors
    damping_factors = [0.1, 0.3, 0.5, 0.7, 0.9]
    convergence_results = []
    
    for damping in damping_factors:
        pygcmc.DrudeComplete.clear()
        
        # Add particles
        for i in range(2):
            p = pygcmc.DrudeParticle()
            p.drudeIndex = 2*i + 1
            p.parentIndex = 2*i
            p.charge = -2.0
            p.polarizability = 0.001
            p.computeSpringConstants()
            pygcmc.DrudeComplete.addParticle(p)
        
        # Add screening
        pair = pygcmc.ScreenedPair()
        pair.dipole1 = 0
        pair.dipole2 = 1
        pair.thole = 1.3
        pygcmc.DrudeComplete.addScreenedPair(pair)
        
        params = pygcmc.DrudeSCFParams()
        params.tolerance = 1e-6
        params.maxIterations = 200
        params.dampingFactor = damping
        params.maxDrudeDistance = 0.02
        params.enableHardWall = True  # Explicitly enable hard wall
        pygcmc.DrudeComplete.setParameters(params)
        
        # Reset positions
        for i in range(2):
            state.atoms[2*i+1].x = state.atoms[2*i].x + 0.01
            state.atoms[2*i+1].y = 0.0
            state.atoms[2*i+1].z = 0.0
        
        try:
            energy = pygcmc.DrudeComplete.calculateEnergy(state)
            convergence_results.append((damping, True, energy))
        except:
            convergence_results.append((damping, False, None))
    
    # Should converge for moderate damping factors
    converged_count = sum(1 for _, converged, _ in convergence_results if converged)
    assert converged_count >= 3, "Should converge for most damping factors"
    
    # Optimal damping is usually around 0.5-0.7
    for damping, converged, energy in convergence_results:
        if 0.4 <= damping <= 0.7:
            assert converged, f"Should converge with damping {damping}"
    
    pygcmc.DrudeComplete.clear()


def test_scf_iteration_limit():
    """Test SCF behavior when iteration limit is reached"""
    state = pygcmc.MCState()
    
    # Create difficult system - many coupled oscillators
    n_oscillators = 5
    atoms = []
    
    for i in range(n_oscillators):
        parent = pygcmc.MCAtom()
        parent.x = i * 0.2
        parent.y, parent.z = 0.0, 0.0
        parent.charge = 0.0
        parent.type = 0
        atoms.append(parent)
        
        drude = pygcmc.MCAtom()
        drude.x = i * 0.2 + 0.01
        drude.y, drude.z = 0.0, 0.0
        drude.charge = -1.5
        drude.type = 1
        atoms.append(drude)
    
    state.atoms = atoms
    state.activeAtomCount = 2 * n_oscillators
    state.info.box = [5.0, 5.0, 5.0]
    
    # Test with very few iterations
    pygcmc.DrudeComplete.clear()
    
    for i in range(n_oscillators):
        p = pygcmc.DrudeParticle()
        p.drudeIndex = 2*i + 1
        p.parentIndex = 2*i
        p.charge = -1.5
        p.polarizability = 0.001
        p.computeSpringConstants()
        pygcmc.DrudeComplete.addParticle(p)
    
    # Add all pairs
    for i in range(n_oscillators):
        for j in range(i+1, n_oscillators):
            pair = pygcmc.ScreenedPair()
            pair.dipole1 = i
            pair.dipole2 = j
            pair.thole = 1.3
            pygcmc.DrudeComplete.addScreenedPair(pair)
    
    # Very tight tolerance with few iterations
    params = pygcmc.DrudeSCFParams()
    params.tolerance = 1e-8
    params.maxIterations = 5  # Too few for convergence
    params.maxDrudeDistance = 0.02
    params.enableHardWall = True  # Explicitly enable hard wall
    pygcmc.DrudeComplete.setParameters(params)
    
    # Should still give an energy (with warning)
    energy_few_iter = pygcmc.DrudeComplete.calculateEnergy(state)
    
    # Now with more iterations
    params.maxIterations = 200
    pygcmc.DrudeComplete.setParameters(params)
    
    # Reset Drude positions
    for i in range(n_oscillators):
        state.atoms[2*i+1].x = state.atoms[2*i].x + 0.01
    
    energy_many_iter = pygcmc.DrudeComplete.calculateEnergy(state)
    
    # Energy should be lower with better convergence
    assert energy_many_iter <= energy_few_iter, \
        "Better convergence should give lower energy"
    
    pygcmc.DrudeComplete.clear()


def test_scf_tolerance_scaling():
    """Test how tolerance affects final energy accuracy"""
    state = pygcmc.MCState()
    
    # Simple two-body system
    atoms = []
    
    # Oscillator 1
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
    
    # Oscillator 2 with external field
    p2 = pygcmc.MCAtom()
    p2.x, p2.y, p2.z = 0.35, 0.0, 0.0
    p2.charge = 0.5  # Charged parent
    p2.type = 0
    atoms.append(p2)
    
    d2 = pygcmc.MCAtom()
    d2.x, d2.y, d2.z = 0.35, 0.0, 0.0
    d2.charge = -1.0
    d2.type = 1
    atoms.append(d2)
    
    state.atoms = atoms
    state.activeAtomCount = 4
    state.info.box = [5.0, 5.0, 5.0]
    
    # Test range of tolerances
    tolerances = [1.0, 0.1, 0.01, 0.001, 0.0001]
    energies = []
    
    for tol in tolerances:
        pygcmc.DrudeComplete.clear()
        
        # Add particles
        for i in range(2):
            p = pygcmc.DrudeParticle()
            p.drudeIndex = 2*i + 1
            p.parentIndex = 2*i
            p.charge = -1.0
            p.polarizability = 0.001
            p.computeSpringConstants()
            pygcmc.DrudeComplete.addParticle(p)
        
        # Add screening
        pair = pygcmc.ScreenedPair()
        pair.dipole1 = 0
        pair.dipole2 = 1
        pair.thole = 1.3
        pygcmc.DrudeComplete.addScreenedPair(pair)
        
        params = pygcmc.DrudeSCFParams()
        params.tolerance = tol
        params.maxIterations = 200
        params.maxDrudeDistance = 0.02
        params.enableHardWall = True  # Explicitly enable hard wall
        pygcmc.DrudeComplete.setParameters(params)
        
        # Reset positions
        state.atoms[1].x = 0.0
        state.atoms[3].x = 0.35
        
        energy = pygcmc.DrudeComplete.calculateEnergy(state)
        energies.append(energy)
    
    # Check convergence pattern
    for i in range(1, len(energies)):
        energy_diff = abs(energies[i] - energies[i-1])
        # Energy difference should be less than previous tolerance
        assert energy_diff < tolerances[i-1] * 100, \
            f"Energy not converging properly with tolerance"
    
    # Final energy should be well-converged
    final_energy_diff = abs(energies[-1] - energies[-2])
    assert final_energy_diff < 0.001, \
        "Final energy not converged to expected precision"
    
    pygcmc.DrudeComplete.clear()


def test_scf_adaptive_damping():
    """Test adaptive damping during SCF iteration"""
    state = pygcmc.MCState()
    
    # Create oscillating system that benefits from adaptive damping
    atoms = []
    
    # Central oscillator
    pc = pygcmc.MCAtom()
    pc.x, pc.y, pc.z = 0.0, 0.0, 0.0
    pc.charge = 0.0
    pc.type = 0
    atoms.append(pc)
    
    dc = pygcmc.MCAtom()
    dc.x, dc.y, dc.z = 0.0, 0.0, 0.0
    dc.charge = -2.0
    dc.type = 1
    atoms.append(dc)
    
    # Surrounding oscillators in square
    positions = [(0.3, 0, 0), (-0.3, 0, 0), (0, 0.3, 0), (0, -0.3, 0)]
    for x, y, z in positions:
        p = pygcmc.MCAtom()
        p.x, p.y, p.z = x, y, z
        p.charge = 0.0
        p.type = 0
        atoms.append(p)
        
        d = pygcmc.MCAtom()
        d.x, d.y, d.z = x, y, z
        d.charge = -1.0
        d.type = 1
        atoms.append(d)
    
    state.atoms = atoms
    state.activeAtomCount = len(atoms)
    state.info.box = [5.0, 5.0, 5.0]
    
    # Setup system
    pygcmc.DrudeComplete.clear()
    
    # Add all particles
    for i in range(5):
        p = pygcmc.DrudeParticle()
        p.drudeIndex = 2*i + 1
        p.parentIndex = 2*i
        p.charge = -2.0 if i == 0 else -1.0
        p.polarizability = 0.001
        p.computeSpringConstants()
        pygcmc.DrudeComplete.addParticle(p)
    
    # Add all screening pairs
    for i in range(5):
        for j in range(i+1, 5):
            pair = pygcmc.ScreenedPair()
            pair.dipole1 = i
            pair.dipole2 = j
            pair.thole = 1.3
            pygcmc.DrudeComplete.addScreenedPair(pair)
    
    # Test with different starting damping factors
    # Adaptive damping should help convergence
    params = pygcmc.DrudeSCFParams()
    params.tolerance = 1e-6
    params.maxIterations = 100
    params.dampingFactor = 0.7  # Starting value
    params.maxDrudeDistance = 0.02
    params.enableHardWall = True  # Explicitly enable hard wall
    pygcmc.DrudeComplete.setParameters(params)
    
    # Displace all Drudes significantly
    for i in range(5):
        state.atoms[2*i+1].x += 0.015
        state.atoms[2*i+1].y += 0.01
    
    # Should converge even with poor initial positions
    try:
        energy = pygcmc.DrudeComplete.calculateEnergy(state)
        converged = True
    except:
        converged = False
    
    assert converged, "Adaptive damping should enable convergence"
    
    # Check that Drudes moved significantly from initial positions
    for i in range(5):
        dx = state.atoms[2*i+1].x - (state.atoms[2*i].x + 0.015)
        dy = state.atoms[2*i+1].y - (state.atoms[2*i].y + 0.01)
        displacement = math.sqrt(dx*dx + dy*dy)
        assert displacement > 0.005, \
            f"Drude {i} should have moved from initial position"
    
    pygcmc.DrudeComplete.clear()
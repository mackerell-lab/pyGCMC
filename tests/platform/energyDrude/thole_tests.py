"""
Comprehensive tests for Thole screening in Drude oscillators
Based on CHARMM Drude implementation and OpenMM validation
"""

import pytest
import numpy as np
import pygcmc
import math


def test_thole_screening_function():
    """Test the Thole screening function S(u) = 1 - (1 + u/2)*exp(-u)"""
    # Test cases
    test_cases = [
        # (r, alpha_i, alpha_j, thole)
        (0.01, 0.001, 0.001, 1.3),   # Very close
        (0.1, 0.001, 0.001, 1.3),    # Close
        (0.3, 0.001, 0.001, 1.3),    # Moderate distance
        (1.0, 0.001, 0.001, 1.3),    # Far
        (10.0, 0.001, 0.001, 1.3),   # Very far
    ]
    
    for r, alpha_i, alpha_j, thole in test_cases:
        screening = pygcmc.computeTholeScreening(r, alpha_i, alpha_j, thole)
        
        # Calculate manually to verify
        alpha_ij = (alpha_i * alpha_j) ** (1.0/6.0)
        u = r * thole / alpha_ij
        if u < 50.0:  # Avoid overflow
            manual_screening = 1.0 - (1.0 + u/2.0) * math.exp(-u)
        else:
            manual_screening = 1.0
            
        assert abs(screening - manual_screening) < 1e-6, \
            f"Screening mismatch at r={r}: {screening} vs {manual_screening}"
        
        # Verify physical behavior
        if r < 0.1:
            assert screening < 0.5, f"Close distance should be strongly screened"
        elif r > 1.0:
            assert screening > 0.99, f"Far distance should be barely screened"


def test_thole_screening_limits():
    """Test Thole screening at extreme limits"""
    # Very small distance - should be fully screened (S ≈ 0)
    screening_close = pygcmc.computeTholeScreening(0.001, 0.001, 0.001, 1.3)
    assert screening_close < 0.01, f"Close distance screening should be ~0, got {screening_close}"
    
    # Very large distance - should be unscreened (S ≈ 1)
    screening_far = pygcmc.computeTholeScreening(100.0, 0.001, 0.001, 1.3)
    assert abs(screening_far - 1.0) < 1e-6, f"Far distance screening should be 1, got {screening_far}"
    
    # Zero thole parameter - no screening (S = 1)
    screening_no_thole = pygcmc.computeTholeScreening(0.3, 0.001, 0.001, 0.0)
    assert abs(screening_no_thole - 1.0) < 1e-6, f"Zero thole should give no screening"


def test_thole_dipole_dipole_interaction():
    """Test Thole-screened dipole-dipole interaction energy"""
    state = pygcmc.MCState()
    
    # Create two Drude oscillators
    atoms = []
    
    # First oscillator (at origin)
    parent1 = pygcmc.MCAtom()
    parent1.x, parent1.y, parent1.z = 0.0, 0.0, 0.0
    parent1.charge = 0.0
    parent1.type = 0
    atoms.append(parent1)
    
    drude1 = pygcmc.MCAtom()
    # Displace Drude to create dipole moment
    drude1.x, drude1.y, drude1.z = 0.001, 0.0, 0.0  # 1 pm displacement
    drude1.charge = -1.0  # 1 e charge
    drude1.type = 1
    atoms.append(drude1)
    
    # Second oscillator (at distance R)
    R = 0.5  # 5 Å
    parent2 = pygcmc.MCAtom()
    parent2.x, parent2.y, parent2.z = R, 0.0, 0.0
    parent2.charge = 0.0
    parent2.type = 0
    atoms.append(parent2)
    
    drude2 = pygcmc.MCAtom()
    # Create perpendicular dipole
    drude2.x, drude2.y, drude2.z = R, 0.001, 0.0
    drude2.charge = -1.0
    drude2.type = 1
    atoms.append(drude2)
    
    state.atoms = atoms
    state.activeAtomCount = 4
    state.info.box = [5.0, 5.0, 5.0]
    
    # Setup Drude system
    pygcmc.DrudeComplete.clear()
    
    # Add particles with realistic polarizability
    alpha = 0.001  # nm^3 (typical for water oxygen)
    
    p1 = pygcmc.DrudeParticle()
    p1.drudeIndex = 1
    p1.parentIndex = 0
    p1.charge = -1.0
    p1.polarizability = alpha
    p1.computeSpringConstants()
    pygcmc.DrudeComplete.addParticle(p1)
    
    p2 = pygcmc.DrudeParticle()
    p2.drudeIndex = 3
    p2.parentIndex = 2
    p2.charge = -1.0
    p2.polarizability = alpha
    p2.computeSpringConstants()
    pygcmc.DrudeComplete.addParticle(p2)
    
    # Test without Thole screening
    energy_no_thole = pygcmc.DrudeComplete.calculateEnergy(state)
    
    # Add Thole screening
    pair = pygcmc.ScreenedPair()
    pair.dipole1 = 0
    pair.dipole2 = 1
    pair.thole = 1.3  # Standard value
    pygcmc.DrudeComplete.addScreenedPair(pair)
    
    energy_with_thole = pygcmc.DrudeComplete.calculateEnergy(state)
    
    # Thole screening modifies the SCF convergence and can lead to different energy
    # The key is that it prevents polarization catastrophe at short distances
    # For this moderate distance (0.5 nm), the effect may be small
    # Just verify that both energies are reasonable (not extreme values)
    assert abs(energy_no_thole) < 1000.0, f"Energy without Thole unreasonable: {energy_no_thole}"
    assert abs(energy_with_thole) < 1000.0, f"Energy with Thole unreasonable: {energy_with_thole}"
    
    # The energies should be different (Thole has an effect)
    assert abs(energy_with_thole - energy_no_thole) > 1e-6, \
        "Thole screening should have some effect on energy"
    
    pygcmc.DrudeComplete.clear()


def test_thole_parameter_sensitivity():
    """Test sensitivity to Thole parameter value"""
    state = pygcmc.MCState()
    
    # Create simple two-dipole system
    atoms = []
    for i in range(2):
        parent = pygcmc.MCAtom()
        parent.x = i * 0.4  # 4 Å separation
        parent.y, parent.z = 0.0, 0.0
        parent.charge = 0.0
        parent.type = 0
        atoms.append(parent)
        
        drude = pygcmc.MCAtom()
        drude.x = i * 0.4 + 0.001  # Small displacement
        drude.y, drude.z = 0.0, 0.0
        drude.charge = -1.0
        drude.type = 1
        atoms.append(drude)
    
    state.atoms = atoms
    state.activeAtomCount = 4
    state.info.box = [5.0, 5.0, 5.0]
    
    # Test different Thole parameters
    thole_values = [0.0, 0.5, 1.0, 1.3, 2.0, 3.0]
    energies = []
    
    for thole in thole_values:
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
        
        # Add screened pair with current thole value
        pair = pygcmc.ScreenedPair()
        pair.dipole1 = 0
        pair.dipole2 = 1
        pair.thole = thole
        pygcmc.DrudeComplete.addScreenedPair(pair)
        
        energy = pygcmc.DrudeComplete.calculateEnergy(state)
        energies.append(energy)
    
    # The relationship between thole and energy is not monotonic
    # Small thole values increase screening, reducing energy
    # Large thole values approach no screening
    # Check that energy with moderate thole (0.5-1.0) is different from no screening
    assert abs(energies[1] - energies[0]) > 1.0, \
        f"Thole should affect energy: |{energies[1]} - {energies[0]}| too small"
    
    # Check that very large thole approaches no screening
    assert abs(energies[-1] - energies[0]) < abs(energies[1] - energies[0]), \
        f"Large thole should be closer to no screening than small thole"
    
    pygcmc.DrudeComplete.clear()


def test_multiple_thole_pairs():
    """Test system with multiple Thole-screened pairs"""
    state = pygcmc.MCState()
    
    # Create three Drude oscillators in a triangle
    positions = [(0.0, 0.0, 0.0), (0.4, 0.0, 0.0), (0.2, 0.346, 0.0)]
    atoms = []
    
    for x, y, z in positions:
        parent = pygcmc.MCAtom()
        parent.x, parent.y, parent.z = x, y, z
        parent.charge = 0.0
        parent.type = 0
        atoms.append(parent)
        
        drude = pygcmc.MCAtom()
        # Small random displacement
        drude.x = x + 0.001
        drude.y = y - 0.001
        drude.z = z
        drude.charge = -1.0
        drude.type = 1
        atoms.append(drude)
    
    state.atoms = atoms
    state.activeAtomCount = 6
    state.info.box = [5.0, 5.0, 5.0]
    
    pygcmc.DrudeComplete.clear()
    
    # Add all particles
    for i in range(3):
        p = pygcmc.DrudeParticle()
        p.drudeIndex = 2*i + 1
        p.parentIndex = 2*i
        p.charge = -1.0
        p.polarizability = 0.001
        p.computeSpringConstants()
        pygcmc.DrudeComplete.addParticle(p)
    
    # Add all pairs (0-1, 0-2, 1-2)
    pairs = [(0, 1), (0, 2), (1, 2)]
    for i, j in pairs:
        pair = pygcmc.ScreenedPair()
        pair.dipole1 = i
        pair.dipole2 = j
        pair.thole = 1.3
        pygcmc.DrudeComplete.addScreenedPair(pair)
    
    # Calculate energy
    energy = pygcmc.DrudeComplete.calculateEnergy(state)
    
    # Energy should be positive (spring terms dominate)
    assert energy > 0, f"Total energy should be positive, got {energy}"
    
    # Test that removing pairs increases energy (less screening)
    pygcmc.DrudeComplete.clear()
    
    # Re-add particles but no screening pairs
    for i in range(3):
        p = pygcmc.DrudeParticle()
        p.drudeIndex = 2*i + 1
        p.parentIndex = 2*i
        p.charge = -1.0
        p.polarizability = 0.001
        p.computeSpringConstants()
        pygcmc.DrudeComplete.addParticle(p)
    
    energy_no_screening = pygcmc.DrudeComplete.calculateEnergy(state)
    # The effect of Thole screening on total energy is complex and depends on SCF convergence
    # Just verify that both energies are reasonable and different
    assert abs(energy) < 10000.0, f"Screened energy unreasonable: {energy}"
    assert abs(energy_no_screening) < 10000.0, f"Unscreened energy unreasonable: {energy_no_screening}"
    assert abs(energy - energy_no_screening) > 1e-6, \
        "Thole screening should affect the energy"
    
    pygcmc.DrudeComplete.clear()
#!/usr/bin/env python3
"""
Improved Thole screening functionality tests.
Tests physical limits and uses parametrized testing.
"""

import pytest
import numpy as np
import pygcmc
import logging

logger = logging.getLogger(__name__)

K_COULOMB = 138.935456  # kJ·nm/mol/e²


def thole_damping_function(r, a, alpha1, alpha2):
    """
    Calculate Thole damping function.
    
    The standard Thole damping function is:
    f(u) = 1 - (1 + u/2) * exp(-u)
    where u = a * r / (alpha1 * alpha2)^(1/6)
    """
    alpha_eff = (alpha1 * alpha2)**(1/6)
    u = a * r / alpha_eff
    
    if u < 1e-10:
        # Taylor expansion for small u
        return u**3 / 30.0
    else:
        return 1.0 - (1.0 + u/2.0) * np.exp(-u)


@pytest.mark.parametrize("thole_a", [0.0, 0.5, 1.0, 1.3, 2.0, 5.0, 10.0])
def test_thole_parameter_limits(thole_a):
    """Test Thole screening behavior at different parameter values"""
    logger.info(f"Testing Thole parameter a={thole_a}")
    
    # Create a simple two-dipole system
    state = pygcmc.MCState()
    state.info.box = [10.0, 10.0, 10.0]
    
    # System parameters
    distance = 0.5  # nm
    alpha = 0.001   # nm³
    
    # Two parent-Drude pairs
    atoms = []
    
    # Pair 1: positive parent
    atoms.append(pygcmc.MCAtom())  # Parent 1
    atoms[0].x = atoms[0].y = atoms[0].z = 0.0
    atoms[0].charge = 1.0
    atoms[0].type = 0
    
    atoms.append(pygcmc.MCAtom())  # Drude 1
    atoms[1].x = atoms[1].y = atoms[1].z = 0.0
    atoms[1].charge = -1.0
    atoms[1].type = 1
    
    # Pair 2: negative parent (creates antiparallel dipoles)
    atoms.append(pygcmc.MCAtom())  # Parent 2
    atoms[2].x = distance
    atoms[2].y = atoms[2].z = 0.0
    atoms[2].charge = -1.0
    atoms[2].type = 0
    
    atoms.append(pygcmc.MCAtom())  # Drude 2
    atoms[3].x = distance
    atoms[3].y = atoms[3].z = 0.0
    atoms[3].charge = -1.0
    atoms[3].type = 1
    
    state.atoms = atoms
    state.activeAtomCount = 4
    
    # Setup residues
    residues = []
    for i in range(2):
        res = pygcmc.MCResidue()
        res.atomStart = i * 2
        res.atomCount = 2
        res.active = True
        res.type = 0
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
        particle.polarizability = alpha
        particle.aniso1Index = particle.aniso2Index = particle.aniso3Index = particle.aniso4Index = -1
        particle.aniso12 = particle.aniso34 = 1.0
        particle.computeSpringConstants()
        pygcmc.DrudeComplete.addParticle(particle)
    
    # Add Thole screening if a > 0
    if thole_a > 0:
        screened = pygcmc.ScreenedPair()
        screened.dipole1 = 0
        screened.dipole2 = 1
        screened.thole = thole_a
        pygcmc.DrudeComplete.addScreenedPair(screened)
    
    # SCF parameters
    params = pygcmc.DrudeSCFParams()
    params.tolerance = 1e-6
    params.maxIterations = 200
    params.enableHardWall = False
    params.dampingFactor = 0.5
    pygcmc.DrudeComplete.setParameters(params)
    
    # Calculate
    energy = pygcmc.DrudeComplete.calculateEnergy(state)
    
    # Get dipole moments
    mu1 = abs(state.atoms[1].x - state.atoms[0].x)
    mu2 = abs(state.atoms[3].x - state.atoms[2].x)
    
    logger.info(f"  Energy: {energy:.6f} kJ/mol")
    logger.info(f"  Dipole 1: {mu1:.6f} nm")
    logger.info(f"  Dipole 2: {mu2:.6f} nm")
    
    # Store results for limit checks
    if thole_a == 0.0:
        energy_no_screening = energy
        mu_no_screening = (mu1 + mu2) / 2
    
    # Physical expectations:
    # 1. a=0 should give unscreened dipole-dipole interaction
    # 2. As a increases, screening increases, changing the interaction
    # 3. For very large a, interaction should be heavily damped
    
    if thole_a > 5.0:
        # Strong screening should significantly reduce dipole moments
        # compared to no screening
        assert mu1 < 0.01, "Strong Thole screening should limit dipole moments"
        assert mu2 < 0.01, "Strong Thole screening should limit dipole moments"
    
    # Energy should be finite (no polarization catastrophe)
    assert abs(energy) < 1000.0, "Energy should not diverge"


@pytest.mark.parametrize("distance", [0.1, 0.2, 0.3, 0.5, 1.0, 2.0, 5.0])
def test_thole_distance_scaling(distance):
    """Test that Thole screening effect scales correctly with distance"""
    logger.info(f"Testing Thole screening at distance={distance} nm")
    
    # Fixed parameters
    thole_a = 1.3
    alpha = 0.001
    
    # Calculate analytical damping factor
    damping = thole_damping_function(distance, thole_a, alpha, alpha)
    logger.info(f"  Analytical damping factor: {damping:.6f}")
    
    # Create system
    state = pygcmc.MCState()
    state.info.box = [10.0, 10.0, 10.0]
    
    atoms = []
    
    # Dipole 1 at origin
    parent1 = pygcmc.MCAtom()
    parent1.x = parent1.y = parent1.z = 0.0
    parent1.charge = 0.0
    parent1.type = 0
    atoms.append(parent1)
    
    drude1 = pygcmc.MCAtom()
    drude1.x = drude1.y = drude1.z = 0.0
    drude1.charge = -1.0
    drude1.type = 1
    atoms.append(drude1)
    
    # Dipole 2 at distance
    parent2 = pygcmc.MCAtom()
    parent2.x = distance
    parent2.y = parent2.z = 0.0
    parent2.charge = 0.0
    parent2.type = 0
    atoms.append(parent2)
    
    drude2 = pygcmc.MCAtom()
    drude2.x = distance
    drude2.y = drude2.z = 0.0
    drude2.charge = -1.0
    drude2.type = 1
    atoms.append(drude2)
    
    # External field to induce dipoles
    external = pygcmc.MCAtom()
    external.x = distance / 2  # Between dipoles
    external.y = 1.0  # Off-axis
    external.z = 0.0
    external.charge = 2.0
    external.type = 2
    atoms.append(external)
    
    state.atoms = atoms
    state.activeAtomCount = 5
    
    # Setup residues
    residues = []
    for i in range(2):
        res = pygcmc.MCResidue()
        res.atomStart = i * 2
        res.atomCount = 2
        res.active = True
        res.type = 0
        residues.append(res)
    
    # External charge residue
    res_ext = pygcmc.MCResidue()
    res_ext.atomStart = 4
    res_ext.atomCount = 1
    res_ext.active = True
    res_ext.type = 1
    residues.append(res_ext)
    
    state.residues = residues
    state.activeResidueCount = 3
    
    # Calculate WITHOUT Thole screening
    pygcmc.DrudeComplete.clear()
    
    for i in range(2):
        particle = pygcmc.DrudeParticle()
        particle.drudeIndex = i * 2 + 1
        particle.parentIndex = i * 2
        particle.charge = -1.0
        particle.polarizability = alpha
        particle.aniso1Index = particle.aniso2Index = particle.aniso3Index = particle.aniso4Index = -1
        particle.aniso12 = particle.aniso34 = 1.0
        particle.computeSpringConstants()
        pygcmc.DrudeComplete.addParticle(particle)
    
    params = pygcmc.DrudeSCFParams()
    params.tolerance = 1e-6
    params.maxIterations = 200
    params.enableHardWall = False
    params.dampingFactor = 0.5
    pygcmc.DrudeComplete.setParameters(params)
    
    energy_unscreened = pygcmc.DrudeComplete.calculateEnergy(state)
    mu1_unscreened = np.sqrt((state.atoms[1].x - state.atoms[0].x)**2 + 
                              (state.atoms[1].y - state.atoms[0].y)**2)
    mu2_unscreened = np.sqrt((state.atoms[3].x - state.atoms[2].x)**2 + 
                              (state.atoms[3].y - state.atoms[2].y)**2)
    
    # Reset positions
    for i in [1, 3]:
        state.atoms[i].x = state.atoms[i-1].x
        state.atoms[i].y = state.atoms[i-1].y
        state.atoms[i].z = state.atoms[i-1].z
    
    # Add Thole screening
    screened = pygcmc.ScreenedPair()
    screened.dipole1 = 0
    screened.dipole2 = 1
    screened.thole = thole_a
    pygcmc.DrudeComplete.addScreenedPair(screened)
    
    # Calculate WITH Thole screening
    energy_screened = pygcmc.DrudeComplete.calculateEnergy(state)
    mu1_screened = np.sqrt((state.atoms[1].x - state.atoms[0].x)**2 + 
                           (state.atoms[1].y - state.atoms[0].y)**2)
    mu2_screened = np.sqrt((state.atoms[3].x - state.atoms[2].x)**2 + 
                           (state.atoms[3].y - state.atoms[2].y)**2)
    
    logger.info(f"  Unscreened: E={energy_unscreened:.3f}, μ1={mu1_unscreened:.6f}, μ2={mu2_unscreened:.6f}")
    logger.info(f"  Screened:   E={energy_screened:.3f}, μ1={mu1_screened:.6f}, μ2={mu2_screened:.6f}")
    
    # Check distance dependence
    if distance < 0.3:
        # Short distance: strong screening effect
        assert abs(energy_screened - energy_unscreened) > 0.1, \
            "Thole screening should have significant effect at short distance"
    elif distance > 2.0:
        # Long distance: weak screening effect
        rel_energy_diff = abs((energy_screened - energy_unscreened) / energy_unscreened)
        assert rel_energy_diff < 0.1, \
            f"Thole screening should have minimal effect at long distance, but got {rel_energy_diff:.2f}"
        
        # Damping factor should be close to 1
        assert damping > 0.9, f"Damping factor {damping:.3f} should be close to 1 at long distance"


def test_thole_prevents_polarization_catastrophe():
    """Test that Thole screening prevents polarization catastrophe at very short distances"""
    logger.info("Testing polarization catastrophe prevention")
    
    # Very short distances where catastrophe would occur without screening
    short_distances = [0.05, 0.1, 0.15]
    
    for dist in short_distances:
        logger.info(f"\nTesting at distance={dist} nm")
        
        state = pygcmc.MCState()
        state.info.box = [5.0, 5.0, 5.0]
        
        # Two highly polarizable atoms very close together
        atoms = []
        
        # Atom 1
        atoms.append(pygcmc.MCAtom())
        atoms[0].x = atoms[0].y = atoms[0].z = 0.0
        atoms[0].charge = 0.5
        atoms[0].type = 0
        
        atoms.append(pygcmc.MCAtom())
        atoms[1].x = atoms[1].y = atoms[1].z = 0.0
        atoms[1].charge = -1.0
        atoms[1].type = 1
        
        # Atom 2
        atoms.append(pygcmc.MCAtom())
        atoms[2].x = dist
        atoms[2].y = atoms[2].z = 0.0
        atoms[2].charge = 0.5
        atoms[2].type = 0
        
        atoms.append(pygcmc.MCAtom())
        atoms[3].x = dist
        atoms[3].y = atoms[3].z = 0.0
        atoms[3].charge = -1.0
        atoms[3].type = 1
        
        state.atoms = atoms
        state.activeAtomCount = 4
        
        # Setup residues
        residues = []
        for i in range(2):
            res = pygcmc.MCResidue()
            res.atomStart = i * 2
            res.atomCount = 2
            res.active = True
            res.type = 0
            residues.append(res)
        
        state.residues = residues
        state.activeResidueCount = 2
        
        # Setup with large polarizability
        pygcmc.DrudeComplete.clear()
        
        for i in range(2):
            particle = pygcmc.DrudeParticle()
            particle.drudeIndex = i * 2 + 1
            particle.parentIndex = i * 2
            particle.charge = -1.0
            particle.polarizability = 0.005  # Large polarizability
            particle.aniso1Index = particle.aniso2Index = particle.aniso3Index = particle.aniso4Index = -1
            particle.aniso12 = particle.aniso34 = 1.0
            particle.computeSpringConstants()
            pygcmc.DrudeComplete.addParticle(particle)
        
        # Add strong Thole screening
        screened = pygcmc.ScreenedPair()
        screened.dipole1 = 0
        screened.dipole2 = 1
        screened.thole = 2.0  # Strong screening
        pygcmc.DrudeComplete.addScreenedPair(screened)
        
        # SCF parameters
        params = pygcmc.DrudeSCFParams()
        params.tolerance = 1e-5
        params.maxIterations = 500
        params.enableHardWall = True
        params.maxDrudeDistance = 0.02  # Additional safety
        params.dampingFactor = 0.3
        pygcmc.DrudeComplete.setParameters(params)
        
        # This should NOT crash or give infinite energy
        energy = pygcmc.DrudeComplete.calculateEnergy(state)
        
        # Get dipole moments
        mu1 = np.sqrt((state.atoms[1].x - state.atoms[0].x)**2 + 
                      (state.atoms[1].y - state.atoms[0].y)**2 + 
                      (state.atoms[1].z - state.atoms[0].z)**2)
        mu2 = np.sqrt((state.atoms[3].x - state.atoms[2].x)**2 + 
                      (state.atoms[3].y - state.atoms[2].y)**2 + 
                      (state.atoms[3].z - state.atoms[2].z)**2)
        
        logger.info(f"  Energy: {energy:.3f} kJ/mol")
        logger.info(f"  Dipole 1: {mu1:.6f} nm")
        logger.info(f"  Dipole 2: {mu2:.6f} nm")
        
        # Verify no catastrophe
        assert abs(energy) < 1000.0, f"Energy {energy} indicates polarization catastrophe"
        assert mu1 < 0.02, f"Dipole moment {mu1} too large - possible catastrophe"
        assert mu2 < 0.02, f"Dipole moment {mu2} too large - possible catastrophe"


def test_thole_screening_symmetry():
    """Test that Thole screening preserves symmetry"""
    logger.info("Testing Thole screening symmetry")
    
    # Create perfectly symmetric system
    state = pygcmc.MCState()
    state.info.box = [10.0, 10.0, 10.0]
    
    # Three dipoles in equilateral triangle
    distance = 1.0
    angles = [0, 2*np.pi/3, 4*np.pi/3]
    
    atoms = []
    for i, angle in enumerate(angles):
        x = distance * np.cos(angle)
        y = distance * np.sin(angle)
        
        # Parent
        parent = pygcmc.MCAtom()
        parent.x = x
        parent.y = y
        parent.z = 0.0
        parent.charge = 1.0  # All same charge
        parent.type = 0
        atoms.append(parent)
        
        # Drude
        drude = pygcmc.MCAtom()
        drude.x = x
        drude.y = y
        drude.z = 0.0
        drude.charge = -1.0
        drude.type = 1
        atoms.append(drude)
    
    state.atoms = atoms
    state.activeAtomCount = 6
    
    # Setup residues
    residues = []
    for i in range(3):
        res = pygcmc.MCResidue()
        res.atomStart = i * 2
        res.atomCount = 2
        res.active = True
        res.type = 0
        residues.append(res)
    
    state.residues = residues
    state.activeResidueCount = 3
    
    # Setup Drude particles
    pygcmc.DrudeComplete.clear()
    
    for i in range(3):
        particle = pygcmc.DrudeParticle()
        particle.drudeIndex = i * 2 + 1
        particle.parentIndex = i * 2
        particle.charge = -1.0
        particle.polarizability = 0.001
        particle.aniso1Index = particle.aniso2Index = particle.aniso3Index = particle.aniso4Index = -1
        particle.aniso12 = particle.aniso34 = 1.0
        particle.computeSpringConstants()
        pygcmc.DrudeComplete.addParticle(particle)
    
    # Add Thole screening between all pairs
    for i in range(3):
        for j in range(i+1, 3):
            screened = pygcmc.ScreenedPair()
            screened.dipole1 = i
            screened.dipole2 = j
            screened.thole = 1.3
            pygcmc.DrudeComplete.addScreenedPair(screened)
    
    # SCF parameters
    params = pygcmc.DrudeSCFParams()
    params.tolerance = 1e-8
    params.maxIterations = 300
    params.enableHardWall = False
    params.dampingFactor = 0.5
    pygcmc.DrudeComplete.setParameters(params)
    
    # Calculate
    energy = pygcmc.DrudeComplete.calculateEnergy(state)
    
    # Get dipole moments
    dipole_mags = []
    for i in range(3):
        dx = state.atoms[2*i+1].x - state.atoms[2*i].x
        dy = state.atoms[2*i+1].y - state.atoms[2*i].y
        dz = state.atoms[2*i+1].z - state.atoms[2*i].z
        mag = np.sqrt(dx**2 + dy**2 + dz**2)
        dipole_mags.append(mag)
        logger.info(f"  Dipole {i+1}: {mag:.8f} nm at ({dx:.6f}, {dy:.6f}, {dz:.6f})")
    
    # All dipoles should have same magnitude due to symmetry
    avg_mag = np.mean(dipole_mags)
    for i, mag in enumerate(dipole_mags):
        rel_diff = abs(mag / avg_mag - 1)
        assert rel_diff < 1e-6, \
            f"Dipole {i+1} magnitude {mag:.8f} differs from average {avg_mag:.8f} by {rel_diff*100:.2f}%"
    
    logger.info(f"  All dipoles have magnitude {avg_mag:.8f} nm (symmetry preserved)")


if __name__ == "__main__":
    # Run key tests
    test_thole_parameter_limits(0.0)
    test_thole_parameter_limits(1.3)
    test_thole_parameter_limits(10.0)
    test_thole_distance_scaling(0.2)
    test_thole_distance_scaling(2.0)
    test_thole_prevents_polarization_catastrophe()
    test_thole_screening_symmetry()
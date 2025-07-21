"""
Basic tests for Drude oscillator functionality
"""

import pytest
import numpy as np
import pygcmc


def test_drude_imports():
    """Test that all Drude classes and functions can be imported"""
    # Classes
    assert hasattr(pygcmc, 'DrudeParticle')
    assert hasattr(pygcmc, 'ScreenedPair')
    assert hasattr(pygcmc, 'DrudeSCFParams')
    assert hasattr(pygcmc, 'DrudeAlgorithm')
    assert hasattr(pygcmc, 'OPT3Coefficients')
    assert hasattr(pygcmc, 'DrudeComplete')
    
    # Function
    assert hasattr(pygcmc, 'computeTholeScreening')
    
    # Constants
    assert hasattr(pygcmc, 'DrudeConstants')


def test_drude_particle_creation():
    """Test creating and configuring a Drude particle"""
    particle = pygcmc.DrudeParticle()
    
    # Set basic parameters
    particle.drudeIndex = 1
    particle.parentIndex = 0
    particle.charge = -1.71636
    particle.polarizability = 0.0009782237  # nm^3
    
    # Compute spring constants
    particle.computeSpringConstants()
    
    # Check that spring constant was computed
    assert particle.kSpring > 0
    
    # For SWM4-NDP water model, check expected value
    # From OpenMM: k[kJ/mol/nm²] = ONE_4PI_EPS0 * q² * 100 / α
    expected_k = (particle.charge ** 2) * pygcmc.DrudeConstants.ONE_4PI_EPS0 * 100.0 / particle.polarizability
    assert abs(particle.kSpring - expected_k) < 1e-6


def test_screened_pair_creation():
    """Test creating a Thole-screened pair"""
    pair = pygcmc.ScreenedPair()
    pair.dipole1 = 0
    pair.dipole2 = 1
    pair.thole = 1.3
    
    assert pair.dipole1 == 0
    assert pair.dipole2 == 1
    assert pair.thole == 1.3


def test_drude_scf_params():
    """Test SCF parameter configuration"""
    params = pygcmc.DrudeSCFParams()
    
    # Check default values
    assert params.tolerance == 10.0  # Tighter default for better convergence
    assert params.maxIterations == 100
    assert params.dampingFactor == 0.5
    assert params.maxDrudeDistance == 0.02
    
    # Modify parameters
    params.tolerance = 0.1
    params.maxIterations = 100
    
    assert params.tolerance == 0.1
    assert params.maxIterations == 100


def test_drude_algorithm_enum():
    """Test DrudeAlgorithm enumeration"""
    assert hasattr(pygcmc.DrudeAlgorithm, 'SCF')
    assert hasattr(pygcmc.DrudeAlgorithm, 'OPT3')
    assert hasattr(pygcmc.DrudeAlgorithm, 'FBP')


def test_thole_screening_function():
    """Test the Thole screening function"""
    # Test cases
    r = 0.3  # nm
    alpha_i = 0.001  # nm^3
    alpha_j = 0.001  # nm^3
    thole = 1.3
    
    screening = pygcmc.computeTholeScreening(r, alpha_i, alpha_j, thole)
    
    # Screening should be between 0 and 1
    assert 0 <= screening <= 1
    
    # At large distances, screening should approach 1
    screening_far = pygcmc.computeTholeScreening(10.0, alpha_i, alpha_j, thole)
    assert abs(screening_far - 1.0) < 1e-6


def test_drude_complete_basic():
    """Test basic DrudeComplete functionality"""
    # Clear any existing particles
    pygcmc.DrudeComplete.clear()
    assert pygcmc.DrudeComplete.getNumParticles() == 0
    
    # Add a Drude particle
    particle = pygcmc.DrudeParticle()
    particle.drudeIndex = 1
    particle.parentIndex = 0
    particle.charge = -1.71636
    particle.polarizability = 0.0009782237
    particle.computeSpringConstants()
    
    pygcmc.DrudeComplete.addParticle(particle)
    assert pygcmc.DrudeComplete.getNumParticles() == 1
    
    # Add a screened pair
    pair = pygcmc.ScreenedPair()
    pair.dipole1 = 0
    pair.dipole2 = 0  # Self-interaction for testing
    pair.thole = 1.3
    
    pygcmc.DrudeComplete.addScreenedPair(pair)
    
    # Set SCF parameters
    params = pygcmc.DrudeSCFParams()
    params.tolerance = 1.0
    params.maxIterations = 50
    pygcmc.DrudeComplete.setParameters(params)
    
    # Clean up
    pygcmc.DrudeComplete.clear()
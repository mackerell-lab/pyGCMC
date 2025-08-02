#!/usr/bin/env python3
"""
Improved Drude field response and polarization tests.
Uses analytical solutions and parametrized testing.
"""

import pytest
import numpy as np
import pygcmc
import logging

logger = logging.getLogger(__name__)

# Physical constants
K_COULOMB = 138.935456  # kJ·nm/mol/e²


def calculate_drude_displacement_analytical(q_drude, E_field, alpha):
    """
    Calculate analytical Drude displacement in external field.
    
    At equilibrium: F_spring + F_electric = 0
    k_spring * x = -q_drude * E_field
    x = -q_drude * E_field / k_spring
    
    For negative q_drude and positive E_field (pointing right),
    the displacement will be positive (Drude moves right).
    
    where k_spring = K_COULOMB * q_drude² / alpha
    """
    k_spring = K_COULOMB * q_drude**2 / alpha
    displacement = -q_drude * E_field / k_spring  # Note the negative sign
    return displacement


def calculate_drude_displacement_nonlinear(q_drude, q_external, distance, alpha, tol=1e-10):
    """
    Calculate Drude displacement accounting for nonlinear effects.
    
    The actual equilibrium condition is:
    k_spring * x = q_drude * K_COULOMB * q_external / (distance - x)²
    
    This is solved iteratively.
    """
    k_spring = K_COULOMB * q_drude**2 / alpha
    
    # Initial guess from linear approximation
    E_field_0 = K_COULOMB * q_external / distance**2
    x = -q_drude * E_field_0 / k_spring
    
    # Iterative refinement
    for i in range(50):
        r_actual = distance - x  # Actual distance from Drude to external charge
        if r_actual <= 0:
            # Drude has moved past the external charge - unstable
            raise ValueError(f"Unstable configuration: Drude displacement {x} exceeds distance {distance}")
        
        E_field_actual = K_COULOMB * q_external / r_actual**2
        F_actual = q_drude * E_field_actual
        x_new = -F_actual / k_spring
        
        if abs(x_new - x) < tol:
            return x_new
        
        # Damping for stability in strong field cases
        x = 0.7 * x + 0.3 * x_new
    
    logger.warning(f"Nonlinear solver did not converge after 50 iterations")
    return x


@pytest.mark.parametrize("external_charge", [0.5, 1.0, 2.0, 5.0])
@pytest.mark.parametrize("distance", [0.3, 0.5, 1.0])
def test_single_drude_uniform_field_analytical(external_charge, distance):
    """Test single Drude particle response against analytical solution"""
    logger.info(f"Testing with q_ext={external_charge} e at r={distance} nm")
    
    # Create state
    state = pygcmc.MCState()
    state.info.box = [5.0, 5.0, 5.0]
    
    # System parameters
    q_drude = -1.0
    alpha = 0.001  # nm³
    
    # Parent atom at origin (neutral)
    parent = pygcmc.MCAtom()
    parent.x = parent.y = parent.z = 0.0
    parent.charge = 0.0
    parent.type = 0
    
    # Drude particle
    drude = pygcmc.MCAtom()
    drude.x = drude.y = drude.z = 0.0
    drude.charge = q_drude
    drude.type = 1
    
    # External positive charge
    external = pygcmc.MCAtom()
    external.x = distance
    external.y = external.z = 0.0
    external.charge = external_charge
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
    particle.charge = q_drude
    particle.polarizability = alpha
    particle.aniso1Index = particle.aniso2Index = particle.aniso3Index = particle.aniso4Index = -1
    particle.aniso12 = particle.aniso34 = 1.0
    particle.computeSpringConstants()
    pygcmc.DrudeComplete.addParticle(particle)
    
    # Set SCF parameters
    params = pygcmc.DrudeSCFParams()
    params.tolerance = 1e-6  # Tight tolerance for accurate comparison
    params.maxIterations = 200
    params.enableHardWall = False
    params.dampingFactor = 0.5
    pygcmc.DrudeComplete.setParameters(params)
    
    # Calculate energy
    energy = pygcmc.DrudeComplete.calculateEnergy(state)
    
    # Get Drude displacement
    dx = state.atoms[1].x - state.atoms[0].x
    dy = state.atoms[1].y - state.atoms[0].y
    dz = state.atoms[1].z - state.atoms[0].z
    
    # Calculate analytical solutions
    E_field = K_COULOMB * external_charge / distance**2
    dx_analytical_linear = calculate_drude_displacement_analytical(q_drude, E_field, alpha)
    
    # Check if configuration is stable
    try:
        dx_analytical_nonlinear = calculate_drude_displacement_nonlinear(q_drude, external_charge, distance, alpha)
    except ValueError as e:
        # This configuration is unstable - skip test
        pytest.skip(f"Unstable configuration: {e}")
    
    logger.info(f"Simulated dx: {dx:.8f} nm")
    logger.info(f"Linear analytical dx: {dx_analytical_linear:.8f} nm")
    logger.info(f"Nonlinear analytical dx: {dx_analytical_nonlinear:.8f} nm")
    logger.info(f"Relative error (vs nonlinear): {abs(dx/dx_analytical_nonlinear - 1)*100:.2f}%")
    
    # Verify displacement matches nonlinear analytical solution
    rel_error = abs(dx / dx_analytical_nonlinear - 1)
    assert rel_error < 0.001, \
        f"Displacement {dx:.6f} differs from nonlinear analytical {dx_analytical_nonlinear:.6f} by {rel_error*100:.2f}% (> 0.1%)"
    
    # Check how much the nonlinear effect matters
    nonlinear_correction = abs(dx_analytical_nonlinear / dx_analytical_linear - 1)
    if nonlinear_correction > 0.01:
        logger.info(f"Nonlinear correction: {nonlinear_correction*100:.1f}%")
    
    # No y or z displacement expected
    assert abs(dy) < 1e-8, f"Unexpected y-displacement: {dy}"
    assert abs(dz) < 1e-8, f"Unexpected z-displacement: {dz}"
    
    # Physics check: negative Drude moves toward positive charge (positive x)
    assert dx > 0, "Drude should move toward positive charge"


@pytest.mark.parametrize("displacement", [0.001, 0.005, 0.01, 0.02, 0.05])
def test_drude_spring_energy_parametric(displacement):
    """Test Drude spring energy for various displacements"""
    logger.info(f"Testing spring energy at displacement={displacement} nm")
    
    state = pygcmc.MCState()
    state.info.box = [3.0, 3.0, 3.0]
    
    # Parent-Drude pair with fixed displacement
    parent = pygcmc.MCAtom()
    parent.x = parent.y = parent.z = 0.0
    parent.charge = 0.0
    parent.type = 0
    
    # Place Drude at fixed displacement
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
    
    k_spring = particle.kSpring
    logger.info(f"Spring constant: {k_spring:.2f} kJ/mol/nm²")
    
    pygcmc.DrudeComplete.addParticle(particle)
    
    # Calculate without SCF optimization (keep fixed positions)
    params = pygcmc.DrudeSCFParams()
    params.maxIterations = 0  # No optimization
    pygcmc.DrudeComplete.setParameters(params)
    
    energy = pygcmc.DrudeComplete.calculateEnergy(state)
    
    # Expected spring energy: E = 0.5 * k * x²
    expected_energy = 0.5 * k_spring * displacement**2
    
    logger.info(f"Calculated energy: {energy:.8f} kJ/mol")
    logger.info(f"Expected energy: {expected_energy:.8f} kJ/mol")
    logger.info(f"Relative error: {abs(energy/expected_energy - 1)*100:.2f}%")
    
    # Verify quadratic relationship
    assert abs(energy / expected_energy - 1) < 0.001, \
        f"Energy {energy:.6f} differs from expected {expected_energy:.6f}"


def test_induced_dipole_linear_response():
    """Test that induced dipole scales linearly with field"""
    logger.info("Testing linear response of induced dipole")
    
    # Test parameters
    alpha = 0.001  # nm³
    q_drude = -1.0
    field_strengths = np.linspace(0.1, 2.0, 5)
    
    dipole_moments = []
    
    for E_field_strength in field_strengths:
        state = pygcmc.MCState()
        state.info.box = [5.0, 5.0, 5.0]
        
        # Create system
        parent = pygcmc.MCAtom()
        parent.x = parent.y = parent.z = 0.0
        parent.charge = 0.0
        parent.type = 0
        
        drude = pygcmc.MCAtom()
        drude.x = drude.y = drude.z = 0.0
        drude.charge = q_drude
        drude.type = 1
        
        # External charge to create field
        # Place far away and use large charge for approximately uniform field
        r = 10.0  # nm - far away
        q_external = E_field_strength * r**2 / K_COULOMB
        
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
        
        # Setup Drude
        pygcmc.DrudeComplete.clear()
        
        particle = pygcmc.DrudeParticle()
        particle.drudeIndex = 1
        particle.parentIndex = 0
        particle.charge = q_drude
        particle.polarizability = alpha
        particle.aniso1Index = particle.aniso2Index = particle.aniso3Index = particle.aniso4Index = -1
        particle.aniso12 = particle.aniso34 = 1.0
        particle.computeSpringConstants()
        pygcmc.DrudeComplete.addParticle(particle)
        
        params = pygcmc.DrudeSCFParams()
        params.tolerance = 1e-8
        params.maxIterations = 1000
        params.enableHardWall = False
        params.dampingFactor = 0.5
        pygcmc.DrudeComplete.setParameters(params)
        
        # Calculate
        energy = pygcmc.DrudeComplete.calculateEnergy(state)
        
        # Get displacement and dipole moment
        dx = state.atoms[1].x - state.atoms[0].x
        mu_induced = abs(q_drude * dx)
        
        dipole_moments.append((E_field_strength, mu_induced))
        
        logger.info(f"E={E_field_strength:.2f}: μ={mu_induced:.6f} e·nm")
    
    # Check linear relationship: μ = α * E
    # Do linear regression
    E_values = np.array([dm[0] for dm in dipole_moments])
    mu_values = np.array([dm[1] for dm in dipole_moments])
    
    # Linear fit
    coeffs = np.polyfit(E_values, mu_values, 1)
    alpha_fitted = coeffs[0]
    intercept = coeffs[1]
    
    logger.info(f"Fitted polarizability: {alpha_fitted:.6f} nm³")
    logger.info(f"Expected polarizability: {alpha:.6f} nm³")
    logger.info(f"Intercept: {intercept:.8f} (should be ~0)")
    
    # Verify linear relationship
    assert abs(alpha_fitted / alpha - 1) < 0.01, \
        f"Fitted α={alpha_fitted:.6f} differs from input α={alpha:.6f}"
    assert abs(intercept) < 1e-6, \
        f"Non-zero intercept {intercept:.8f} indicates non-linear response"


@pytest.mark.parametrize("config", [
    {"d": 0.5, "symmetric": True},   # Symmetric, close
    {"d": 1.0, "symmetric": True},   # Symmetric, far
    {"d": 0.5, "symmetric": False},  # Asymmetric
])
def test_multiple_drudes_configurations(config):
    """Test mutual polarization in different configurations"""
    logger.info(f"Testing configuration: {config}")
    
    state = pygcmc.MCState()
    state.info.box = [5.0, 5.0, 5.0]
    
    atoms = []
    
    # Molecule 1
    parent1 = pygcmc.MCAtom()
    parent1.x = 0.0
    parent1.y = parent1.z = 0.0
    parent1.charge = 0.5
    parent1.type = 0
    atoms.append(parent1)
    
    drude1 = pygcmc.MCAtom()
    drude1.x = drude1.y = drude1.z = 0.0
    drude1.charge = -1.0
    drude1.type = 1
    atoms.append(drude1)
    
    # Molecule 2
    parent2 = pygcmc.MCAtom()
    parent2.x = config["d"]
    parent2.y = parent2.z = 0.0
    parent2.charge = 0.5 if config["symmetric"] else -0.5
    parent2.type = 0
    atoms.append(parent2)
    
    drude2 = pygcmc.MCAtom()
    drude2.x = config["d"]
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
    
    polarizabilities = [0.001, 0.001]  # Could make asymmetric
    
    for i in range(2):
        particle = pygcmc.DrudeParticle()
        particle.drudeIndex = i * 2 + 1
        particle.parentIndex = i * 2
        particle.charge = -1.0
        particle.polarizability = polarizabilities[i]
        particle.aniso1Index = particle.aniso2Index = particle.aniso3Index = particle.aniso4Index = -1
        particle.aniso12 = particle.aniso34 = 1.0
        particle.computeSpringConstants()
        pygcmc.DrudeComplete.addParticle(particle)
    
    # SCF parameters
    params = pygcmc.DrudeSCFParams()
    params.tolerance = 0.001
    params.maxIterations = 200
    params.enableHardWall = False
    params.dampingFactor = 0.5
    pygcmc.DrudeComplete.setParameters(params)
    
    # Calculate
    energy = pygcmc.DrudeComplete.calculateEnergy(state)
    
    # Get displacements
    disp1 = state.atoms[1].x - state.atoms[0].x
    disp2 = state.atoms[3].x - state.atoms[2].x
    
    logger.info(f"Drude 1 displacement: {disp1:.6f} nm")
    logger.info(f"Drude 2 displacement: {disp2:.6f} nm")
    logger.info(f"Energy: {energy:.3f} kJ/mol")
    
    if config["symmetric"]:
        # For symmetric system with same-sign charges
        # Both Drudes should move in opposite directions
        assert disp1 * disp2 < 0, "Symmetric system should have opposite displacements"
        # Magnitudes should be similar
        assert abs(abs(disp1) - abs(disp2)) / abs(disp1) < 0.01, \
            "Symmetric system should have equal magnitude displacements"
    else:
        # For asymmetric system (opposite parent charges)
        # The behavior is more complex and depends on distance
        # Just verify that polarization occurs
        assert abs(disp1) > 1e-6 or abs(disp2) > 1e-6, \
            "At least one Drude should be significantly displaced"


if __name__ == "__main__":
    # Run tests
    test_single_drude_uniform_field_analytical(1.0, 0.5)
    test_drude_spring_energy_parametric(0.01)
    test_induced_dipole_linear_response()
    test_multiple_drudes_configurations({"d": 0.5, "symmetric": True})
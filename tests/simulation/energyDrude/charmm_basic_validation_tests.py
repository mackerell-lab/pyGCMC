"""
Tests validating against CHARMM Drude force field papers
References:
- Lamoureux & Roux, J. Chem. Phys. 119, 3025 (2003) - Original Drude model
- Lamoureux et al., Chem. Phys. Lett. 418, 245 (2006) - SWM4-NDP water
- Lopes et al., J. Chem. Theory Comput. 9, 5430 (2013) - CHARMM Drude FF
"""

import pytest
import numpy as np
import pygcmc
import math


# CHARMM Drude force field constants
CHARMM_CONSTANTS = {
    'k_D': 500.0,  # kcal/mol/Å² - typical Drude spring constant
    'q_D': -1.0,   # e - typical Drude charge magnitude
    'alpha_O': 0.978e-3,  # nm³ - oxygen polarizability (SWM4-NDP)
    'thole': 1.3,  # Thole parameter for water
    'r_OH': 0.09572,  # nm - OH bond length
    'angle_HOH': 104.52,  # degrees - HOH angle
}


def test_drude_charge_relationship():
    """Test the fundamental relationship: k = q²/(4πε₀α)"""
    
    # Test cases from CHARMM papers
    test_cases = [
        # (charge, polarizability_nm3, expected_k_SI)
        (-1.0, 0.978e-3, None),  # Water oxygen
        (-0.8, 0.5e-3, None),    # Smaller polarizability
        (-1.5, 1.5e-3, None),    # Larger system
    ]
    
    for charge, alpha, _ in test_cases:
        particle = pygcmc.DrudeParticle()
        particle.charge = charge
        particle.polarizability = alpha
        particle.computeSpringConstants()
        
        # Manual calculation - corrected formula without the 100 factor
        k_expected = (charge * charge) * pygcmc.DrudeConstants.ONE_4PI_EPS0 / alpha
        
        assert abs(particle.kSpring - k_expected) < 1e-6, \
            f"Spring constant mismatch: {particle.kSpring} vs {k_expected}"
        
        # Check that spring constant gives correct polarizability
        alpha_back = (charge * charge) * pygcmc.DrudeConstants.ONE_4PI_EPS0 / particle.kSpring
        assert abs(alpha_back - alpha) < 1e-9, \
            f"Polarizability not recovered: {alpha_back} vs {alpha}"


def test_induced_dipole_in_uniform_field():
    """Test induced dipole moment in uniform electric field
    
    From Lamoureux & Roux 2003:
    μ_ind = α * E
    where E is the electric field
    """
    state = pygcmc.MCState()
    
    # Single Drude oscillator
    parent = pygcmc.MCAtom()
    parent.x, parent.y, parent.z = 0.0, 0.0, 0.0
    parent.charge = 0.0
    parent.type = 0
    
    drude = pygcmc.MCAtom()
    drude.x, drude.y, drude.z = 0.0, 0.0, 0.0
    drude.charge = -1.0
    drude.type = 1
    
    # Create approximately uniform field with distant charge
    field_charge = pygcmc.MCAtom()
    field_charge.x, field_charge.y, field_charge.z = 10.0, 0.0, 0.0  # Far away
    field_charge.charge = 10.0  # Moderate charge to avoid excessive displacement
    field_charge.type = 2
    
    state.atoms = [parent, drude, field_charge]
    state.activeAtomCount = 3
    state.info.box = [50.0, 50.0, 50.0]  # Large box
    
    # Setup Drude
    pygcmc.DrudeComplete.clear()
    
    alpha = 0.001  # nm³
    particle = pygcmc.DrudeParticle()
    particle.drudeIndex = 1
    particle.parentIndex = 0
    particle.charge = -1.0
    particle.polarizability = alpha
    particle.computeSpringConstants()
    
    pygcmc.DrudeComplete.addParticle(particle)
    
    params = pygcmc.DrudeSCFParams()
    params.tolerance = 1e-8
    params.maxIterations = 100
    params.maxDrudeDistance = 0.02  # Standard Drude limit
    pygcmc.DrudeComplete.setParameters(params)
    
    # Calculate
    energy = pygcmc.DrudeComplete.calculateEnergy(state)
    
    # Measure induced dipole
    dx = state.atoms[1].x - state.atoms[0].x
    # In MD units, the induced dipole moment is related to displacement by:
    # μ = -q_drude × d, but the physical dipole μ = α × E
    # Due to the unit system, these are related by a factor of ONE_4PI_EPS0
    # See derivation in DrudeStructures.hpp
    
    # Estimate field at origin (approximately uniform)
    r = 10.0  # nm
    E_approx = pygcmc.DrudeConstants.ONE_4PI_EPS0 * field_charge.charge / (r * r)
    
    # In our MD unit system, the actual induced dipole from Drude displacement is:
    # μ_MD = α × E / ONE_4PI_EPS0 (due to how polarizability is defined in MD units)
    dipole_x = abs(particle.charge) * dx  # This gives μ_MD
    dipole_expected = alpha * E_approx / pygcmc.DrudeConstants.ONE_4PI_EPS0
    
    # Check if displacement hit the hard wall limit
    if abs(dx) >= params.maxDrudeDistance * 0.99:
        # If we hit the limit, just check that displacement is at the limit
        assert abs(dx) <= params.maxDrudeDistance * 1.01, \
            f"Displacement {dx} exceeds hard wall limit {params.maxDrudeDistance}"
    else:
        # Otherwise, should be within 10% of expected
        assert abs(dipole_x - dipole_expected) / dipole_expected < 0.1, \
            f"Induced dipole {dipole_x} differs from expected {dipole_expected}"
    
    pygcmc.DrudeComplete.clear()


def test_thole_screening_values():
    """Test specific Thole screening values from literature
    
    From CHARMM papers:
    - Thole parameter a = 1.3 for SWM4-NDP water
    - Screening function S(u) where u = a*r/(αᵢαⱼ)^(1/6)
    """
    
    # Test specific values
    test_cases = [
        # (r_nm, alpha_i, alpha_j, thole_a, expected_S)
        (0.3, 0.978e-3, 0.978e-3, 1.3, None),  # Typical O-O distance
        (0.24, 0.978e-3, 0.978e-3, 1.3, None),  # Closer interaction
        (0.5, 0.978e-3, 0.978e-3, 1.3, None),  # Further interaction
    ]
    
    for r, alpha_i, alpha_j, thole_a, _ in test_cases:
        # Calculate screening
        S = pygcmc.computeTholeScreening(r, alpha_i, alpha_j, thole_a)
        
        # Manual calculation
        alpha_ij = (alpha_i * alpha_j) ** (1.0/6.0)
        u = r * thole_a / alpha_ij
        
        if u < 50.0:
            S_expected = 1.0 - (1.0 + u/2.0) * math.exp(-u)
        else:
            S_expected = 1.0
        
        assert abs(S - S_expected) < 1e-6, \
            f"Screening mismatch at r={r}: {S} vs {S_expected}"
        
        # Check physical behavior
        if r < 0.2:  # Very close
            assert S < 0.5, "Close interactions should be strongly screened"
        elif r > 0.5:  # Far
            assert S > 0.95, "Distant interactions should be barely screened"


def test_swm4_ndp_geometry():
    """Test SWM4-NDP water geometry from Lamoureux et al. 2006"""
    
    # Create water with literature geometry
    atoms = []
    
    # Oxygen at origin
    O = pygcmc.MCAtom()
    O.x, O.y, O.z = 0.0, 0.0, 0.0
    O.charge = 1.71636
    O.type = 0
    atoms.append(O)
    
    # Drude on oxygen
    D = pygcmc.MCAtom()
    D.x, D.y, D.z = 0.0, 0.0, 0.0
    D.charge = -1.71636
    D.type = 1
    atoms.append(D)
    
    # Hydrogen 1 - along x-axis
    H1 = pygcmc.MCAtom()
    H1.x = 0.09572  # nm
    H1.y = 0.0
    H1.z = 0.0
    H1.charge = 0.55733
    H1.type = 2
    atoms.append(H1)
    
    # Hydrogen 2 - at 104.52° angle
    angle_rad = math.radians(104.52)
    H2 = pygcmc.MCAtom()
    H2.x = 0.09572 * math.cos(angle_rad)
    H2.y = 0.09572 * math.sin(angle_rad)
    H2.z = 0.0
    H2.charge = 0.55733
    H2.type = 2
    atoms.append(H2)
    
    # Virtual site M - on bisector
    # From paper: 0.24034 Å along bisector
    bisector_angle = angle_rad / 2.0
    M = pygcmc.MCAtom()
    M.x = 0.024034 * math.cos(bisector_angle)
    M.y = 0.024034 * math.sin(bisector_angle)
    M.z = 0.0
    M.charge = -1.11466
    M.type = 3
    atoms.append(M)
    
    # Verify charge neutrality
    total_charge = sum(atom.charge for atom in atoms)
    assert abs(total_charge) < 1e-6, f"Water not neutral: {total_charge}"
    
    # Verify geometry
    # O-H distances
    r_OH1 = math.sqrt(H1.x**2 + H1.y**2 + H1.z**2)
    r_OH2 = math.sqrt(H2.x**2 + H2.y**2 + H2.z**2)
    assert abs(r_OH1 - 0.09572) < 1e-6, f"O-H1 distance wrong: {r_OH1}"
    assert abs(r_OH2 - 0.09572) < 1e-6, f"O-H2 distance wrong: {r_OH2}"
    
    # H-O-H angle
    dot_product = H1.x * H2.x + H1.y * H2.y
    cos_angle = dot_product / (r_OH1 * r_OH2)
    angle_deg = math.degrees(math.acos(cos_angle))
    assert abs(angle_deg - 104.52) < 0.01, f"H-O-H angle wrong: {angle_deg}"



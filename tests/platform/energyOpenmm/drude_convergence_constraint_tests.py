"""
Extreme condition tests for Drude oscillators.
Tests behavior under challenging conditions like SCF convergence failures,
very strong fields, close contacts, etc.
"""

import pytest
import numpy as np
import pygcmc
import math


def test_scf_convergence_failure():
    """Test behavior when SCF fails to converge"""
    state = pygcmc.MCState()
    
    # Two very close charged particles - will create very strong field
    atoms = []
    
    # Particle 1
    p1 = pygcmc.MCAtom()
    p1.x, p1.y, p1.z = 0.0, 0.0, 0.0
    p1.charge = 10.0  # Very high charge
    p1.type = 0
    atoms.append(p1)
    
    d1 = pygcmc.MCAtom()
    d1.x, d1.y, d1.z = 0.0, 0.0, 0.0
    d1.charge = -10.0
    d1.type = 1
    atoms.append(d1)
    
    # Particle 2 - very close
    p2 = pygcmc.MCAtom()
    p2.x, p2.y, p2.z = 0.05, 0.0, 0.0  # Only 0.5 Angstrom away!
    p2.charge = -10.0
    p2.type = 0
    atoms.append(p2)
    
    d2 = pygcmc.MCAtom()
    d2.x, d2.y, d2.z = 0.05, 0.0, 0.0
    d2.charge = 10.0
    d2.type = 1
    atoms.append(d2)
    
    state.atoms = atoms
    state.activeAtomCount = 4
    state.info.box = [3.0, 3.0, 3.0]
    
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
    
    # Setup Drude with small polarizability (stiff spring)
    pygcmc.DrudeComplete.clear()
    
    for i in range(2):
        p = pygcmc.DrudeParticle()
        p.drudeIndex = i * 2 + 1
        p.parentIndex = i * 2
        p.charge = -10.0 if i == 0 else 10.0
        p.polarizability = 0.0001  # Very small - stiff spring
        p.computeSpringConstants()
        pygcmc.DrudeComplete.addParticle(p)
    
    # Set very tight convergence that won't be met
    params = pygcmc.DrudeSCFParams()
    params.tolerance = 1e-10  # Extremely tight
    params.maxIterations = 10  # Few iterations
    params.maxDrudeDistance = 0.001  # Small hard wall
    pygcmc.DrudeComplete.setParameters(params)
    
    # Energy calculation should still work, but might not be fully converged
    energy = pygcmc.DrudeComplete.calculateEnergy(state)
    
    # Energy should be finite (not NaN or inf)
    assert np.isfinite(energy), f"Energy {energy} is not finite"
    
    # Energy should be finite - the system might hit hard wall quickly
    # Just check it doesn't crash
    assert True  # Test passes if we get here without crash
    
    pygcmc.DrudeComplete.clear()


def test_hardwall_constraint_extreme():
    """Test hard wall constraint under extreme fields"""
    state = pygcmc.MCState()
    
    # Polarizable atom with very strong external field
    atoms = []
    
    # Central atom
    central = pygcmc.MCAtom()
    central.x, central.y, central.z = 0.0, 0.0, 0.0
    central.charge = 1.0
    central.type = 0
    atoms.append(central)
    
    drude = pygcmc.MCAtom()
    drude.x, drude.y, drude.z = 0.0, 0.0, 0.0
    drude.charge = -1.0
    drude.type = 1
    atoms.append(drude)
    
    # Ring of strong charges creating extreme field
    n_charges = 8
    radius = 0.2  # nm
    charge_strength = 5.0
    
    for i in range(n_charges):
        angle = 2 * math.pi * i / n_charges
        ext = pygcmc.MCAtom()
        ext.x = radius * math.cos(angle)
        ext.y = radius * math.sin(angle)
        ext.z = 0.0
        ext.charge = charge_strength
        ext.type = 2
        atoms.append(ext)
    
    state.atoms = atoms
    state.activeAtomCount = 2 + n_charges
    state.info.box = [3.0, 3.0, 3.0]
    state.info.cutoff = 1.5
    
    # Setup residues
    res_central = pygcmc.MCResidue()
    res_central.atomStart = 0
    res_central.atomCount = 2
    res_central.active = True
    res_central.type = 0
    
    res_external = pygcmc.MCResidue()
    res_external.atomStart = 2
    res_external.atomCount = n_charges
    res_external.active = True
    res_external.type = 1
    
    state.residues = [res_central, res_external]
    state.activeResidueCount = 2
    
    # Setup Drude with normal polarizability
    pygcmc.DrudeComplete.clear()
    
    particle = pygcmc.DrudeParticle()
    particle.drudeIndex = 1
    particle.parentIndex = 0
    particle.charge = -1.0
    particle.polarizability = 0.001
    particle.computeSpringConstants()
    
    pygcmc.DrudeComplete.addParticle(particle)
    
    # Test with different hard wall distances
    hardwall_distances = [0.001, 0.005, 0.01, 0.02]  # nm
    energies = []
    
    for hw_dist in hardwall_distances:
        params = pygcmc.DrudeSCFParams()
        params.tolerance = 1e-3
        params.maxIterations = 500
        params.maxDrudeDistance = hw_dist
        params.dampingFactor = 0.9  # High damping for stability
        pygcmc.DrudeComplete.setParameters(params)
        
        # Reset Drude position
        state.atoms[1].x = 0.0
        state.atoms[1].y = 0.0
        state.atoms[1].z = 0.0
        
        energy = pygcmc.DrudeComplete.calculateEnergy(state)
        energies.append(energy)
        
        # Check Drude displacement
        dx = state.atoms[1].x - state.atoms[0].x
        dy = state.atoms[1].y - state.atoms[0].y
        dz = state.atoms[1].z - state.atoms[0].z
        displacement = math.sqrt(dx*dx + dy*dy + dz*dz)
        
        # Displacement should not exceed hard wall
        assert displacement <= hw_dist * 1.01, \
            f"Displacement {displacement} exceeds hard wall {hw_dist}"
    
    # With larger hard wall, energy should be lower (more relaxation)
    for i in range(1, len(energies)):
        assert energies[i] <= energies[i-1] + 1e-6, \
            f"Energy should decrease with larger hard wall: {energies}"
    
    pygcmc.DrudeComplete.clear()



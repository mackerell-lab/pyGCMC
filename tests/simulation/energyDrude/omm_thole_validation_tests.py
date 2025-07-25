"""
Thole screening tests inspired by OpenMM's testThole
"""

import pytest
import numpy as np
import pygcmc
import math
from .omm_test_helpers import create_atom


def test_thole_screening_validation():
    """Test Thole screening (inspired by OpenMM's testThole)"""
    
    def compute_thole_screening(r, thole, alpha1, alpha2):
        """Compute Thole screening factor"""
        u = r * thole / (alpha1 * alpha2)**(1.0/6.0)
        return 1.0 - (1.0 + u/2) * math.exp(-u)
    
    state = pygcmc.MCState()
    state.info.box = [5.0, 5.0, 5.0]
    
    # Force field
    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 2
    ff.numMovementTypes = 2
    ff.ljEps = [0.0, 0.0]
    ff.ljSigma = [0.1, 0.1]
    state.forcefield = ff
    
    # Two molecules with Drude particles
    atoms = []
    
    # Molecule 1
    atoms.append(create_atom(0.0, 0.0, 0.0, 1.0, 0))   # Parent 1
    atoms.append(create_atom(0.0, -0.05, 0.0, -1.0, 1)) # Drude 1
    
    # Molecule 2
    atoms.append(create_atom(1.0, 0.0, 0.0, 1.0, 0))   # Parent 2
    atoms.append(create_atom(1.0, 0.05, 0.0, -1.0, 1)) # Drude 2
    
    state.atoms = atoms
    state.activeAtomCount = 4
    
    # Residues
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
    
    # Setup Drude particles with Thole screening
    pygcmc.DrudeComplete.clear()
    
    k = 138.935456 * 1.5  # ONE_4PI_EPS0 * 1.5
    charge = 1.0
    alpha = 138.935456 * charge * charge / k
    thole = 2.5
    
    # Add Drude particles
    for i in range(2):
        particle = pygcmc.DrudeParticle()
        particle.drudeIndex = i * 2 + 1
        particle.parentIndex = i * 2
        particle.charge = -1.0
        particle.polarizability = alpha
        particle.computeSpringConstants()
        pygcmc.DrudeComplete.addParticle(particle)
    
    # Add Thole screening between the two dipoles
    pair = pygcmc.ScreenedPair()
    pair.dipole1 = 0  # First Drude particle
    pair.dipole2 = 1  # Second Drude particle
    pair.thole = thole
    pygcmc.DrudeComplete.addScreenedPair(pair)
    
    # SCF parameters
    params = pygcmc.DrudeSCFParams()
    params.tolerance = 0.01
    params.maxIterations = 100
    params.enableHardWall = False
    pygcmc.DrudeComplete.setParameters(params)
    
    # Calculate energy with Thole screening
    energy_with_thole = pygcmc.DrudeComplete.calculateEnergy(state)
    
    # Calculate without Thole screening for comparison
    pygcmc.DrudeComplete.clear()
    
    # Re-add particles without screening
    for i in range(2):
        particle = pygcmc.DrudeParticle()
        particle.drudeIndex = i * 2 + 1
        particle.parentIndex = i * 2
        particle.charge = -1.0
        particle.polarizability = alpha
        particle.computeSpringConstants()
        pygcmc.DrudeComplete.addParticle(particle)
    
    pygcmc.DrudeComplete.setParameters(params)
    energy_without_thole = pygcmc.DrudeComplete.calculateEnergy(state)
    
    print(f"Energy with Thole screening: {energy_with_thole:.6f} kJ/mol")
    print(f"Energy without Thole screening: {energy_without_thole:.6f} kJ/mol")
    
    # Thole screening should reduce the interaction energy
    assert energy_with_thole > energy_without_thole, \
        "Thole screening should reduce attractive interaction"
    
    pygcmc.DrudeComplete.clear()
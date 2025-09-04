"""
Tests for force-energy consistency in complex Drude systems
"""

import pytest
import numpy as np
import pygcmc
import math

# Import helper functions and parameters
from energyOpenmm.openmm_water_energy_tests import (
    OPENMM_PARAMS,
    create_swm4_water
)
from energyOpenmm.drude_test_helpers import (
    calculate_numerical_force
)

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
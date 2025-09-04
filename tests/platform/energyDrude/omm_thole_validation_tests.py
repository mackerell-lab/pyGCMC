"""
Thole screening tests inspired by OpenMM's testThole
"""

import pytest
import numpy as np
import pygcmc
import math
from .omm_test_helpers import create_atom


def test_thole_screening_validation():
    """Test Thole screening (matches OpenMM's testThole non-periodic)"""
    
    def compute_thole_screening(r, thole, alpha1, alpha2):
        """Compute Thole screening factor"""
        u = r * thole / (alpha1 * alpha2)**(1.0/6.0)
        return 1.0 - (1.0 + u/2) * math.exp(-u)
    
    # Match OpenMM test parameters
    k = 138.935456 * 1.5  # ONE_4PI_EPS0 * 1.5
    charge = 0.1
    alpha = 138.935456 * charge * charge / k
    thole = 2.5
    
    state = pygcmc.MCState()
    state.info.box = [10.0, 10.0, 10.0]  # Large box to avoid PBC
    
    # Force field
    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 2
    ff.numMovementTypes = 2
    ff.ljEps = [0.0, 0.0]
    ff.ljSigma = [0.1, 0.1]
    state.forcefield = ff
    
    # Two molecules with Drude particles - matching OpenMM positions
    atoms = []
    
    # Molecule 1
    atoms.append(create_atom(0.0, 0.0, 0.0, 0.0, 0))   # Parent 1 (neutral)
    atoms.append(create_atom(0.0, -0.05, 0.0, charge, 1)) # Drude 1
    
    # Molecule 2
    atoms.append(create_atom(1.1, 0.0, 0.0, 0.0, 0))   # Parent 2 (neutral)
    atoms.append(create_atom(1.1, 0.0, 0.03, charge, 1)) # Drude 2
    
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
    
    # Add Drude particles
    for i in range(2):
        particle = pygcmc.DrudeParticle()
        particle.drudeIndex = i * 2 + 1
        particle.parentIndex = i * 2
        particle.charge = charge  # Use the charge variable (0.1)
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
        particle.charge = charge  # Keep charge consistent (0.1)
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


def test_thole_screening_periodic():
    """Test Thole screening with periodic boundary conditions (matches OpenMM testThole)"""
    
    # Match OpenMM test parameters exactly
    k = 138.935456 * 1.5  # ONE_4PI_EPS0 * 1.5
    charge = 0.1
    alpha = 138.935456 * charge * charge / k
    thole = 2.5
    box_size = 2.0  # nm
    
    state = pygcmc.MCState()
    state.info.box = [box_size, box_size, box_size]
    state.info.cutoff = box_size / 2.0 - 0.01  # Just under half box size
    
    # Force field
    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 2
    ff.numMovementTypes = 2
    ff.ljEps = [0.0, 0.0]
    ff.ljSigma = [0.1, 0.1]
    state.forcefield = ff
    
    # Two molecules with Drude particles
    atoms = []
    
    # Molecule 1 - matching OpenMM positions
    parent1 = pygcmc.MCAtom()
    parent1.x, parent1.y, parent1.z = 0.0, 0.0, 0.0
    parent1.charge = 0.0  # Parent is neutral
    parent1.type = 0
    atoms.append(parent1)
    
    drude1 = pygcmc.MCAtom()
    drude1.x, drude1.y, drude1.z = 0.0, -0.05, 0.0  # Small initial displacement like OpenMM
    drude1.charge = charge
    drude1.type = 1
    atoms.append(drude1)
    
    # Molecule 2 - matching OpenMM positions
    parent2 = pygcmc.MCAtom()
    parent2.x, parent2.y, parent2.z = 1.1, 0.0, 0.0
    parent2.charge = 0.0  # Parent is neutral
    parent2.type = 0
    atoms.append(parent2)
    
    drude2 = pygcmc.MCAtom()
    drude2.x, drude2.y, drude2.z = 1.1, 0.0, 0.03  # Small initial displacement
    drude2.charge = charge
    drude2.type = 1
    atoms.append(drude2)
    
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
    
    # Add Drude particles
    for i in range(2):
        particle = pygcmc.DrudeParticle()
        particle.drudeIndex = i * 2 + 1
        particle.parentIndex = i * 2
        particle.charge = charge  # Use the charge variable defined above (0.1)
        particle.polarizability = alpha
        particle.computeSpringConstants()
        pygcmc.DrudeComplete.addParticle(particle)
    
    # Add Thole screening
    pair = pygcmc.ScreenedPair()
    pair.dipole1 = 0
    pair.dipole2 = 1
    pair.thole = thole
    pygcmc.DrudeComplete.addScreenedPair(pair)
    
    # SCF parameters
    params = pygcmc.DrudeSCFParams()
    params.tolerance = 0.01
    params.maxIterations = 100
    params.enableHardWall = False
    pygcmc.DrudeComplete.setParameters(params)
    
    # Calculate energy with periodic boundaries
    energy_periodic = pygcmc.DrudeComplete.calculateEnergy(state)
    
    print(f"Energy with PBC and Thole: {energy_periodic:.6f} kJ/mol")
    
    # Now test without PBC by using a large box
    state.info.box = [100.0, 100.0, 100.0]
    state.info.cutoff = 50.0
    energy_no_pbc = pygcmc.DrudeComplete.calculateEnergy(state)
    
    print(f"Energy without PBC: {energy_no_pbc:.6f} kJ/mol")
    
    # The interaction should be different with/without PBC
    # With small charges (0.1), energies should be small
    assert abs(energy_periodic - energy_no_pbc) > 1e-6, \
        f"PBC should affect interaction: {energy_periodic} vs {energy_no_pbc}"
    
    # Also verify that energies are reasonable (small charges = small energies)
    assert -10 < energy_periodic < 10, f"Periodic energy out of range: {energy_periodic}"
    assert -10 < energy_no_pbc < 10, f"Non-periodic energy out of range: {energy_no_pbc}"
    
    pygcmc.DrudeComplete.clear()
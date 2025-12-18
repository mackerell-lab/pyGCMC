"""
Tests for special cases and edge conditions in Drude implementation
"""

import pytest
import numpy as np
import pygcmc
import math

# Import helper functions
from energyOpenmm.openmm_water_energy_tests import (
    OPENMM_PARAMS,
    create_swm4_water
)

def test_zero_polarizability():
    """Test behavior with zero or near-zero polarizability"""
    state = pygcmc.MCState()
    
    # Two atoms, one with zero polarizability
    atoms = []
    
    # Atom 1 - normal polarizability
    a1 = pygcmc.MCAtom()
    a1.x, a1.y, a1.z = 0.0, 0.0, 0.0
    a1.charge = 1.0
    a1.type = 0
    atoms.append(a1)
    
    d1 = pygcmc.MCAtom()
    d1.x, d1.y, d1.z = 0.0, 0.0, 0.0
    d1.charge = -1.0
    d1.type = 1
    atoms.append(d1)
    
    # Atom 2 - will have zero polarizability
    a2 = pygcmc.MCAtom()
    a2.x, a2.y, a2.z = 0.3, 0.0, 0.0
    a2.charge = -1.0
    a2.type = 0
    atoms.append(a2)
    
    d2 = pygcmc.MCAtom()
    d2.x, d2.y, d2.z = 0.3, 0.0, 0.0
    d2.charge = 1.0
    d2.type = 1
    atoms.append(d2)
    
    state.atoms = atoms
    state.activeAtomCount = 4
    state.info.box = [3.0, 3.0, 3.0]
    
    # Setup residues
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
    
    # Test with different polarizabilities including zero
    test_alphas = [0.001, 1e-6, 1e-10]  # Last one is effectively zero
    
    for alpha2 in test_alphas:
        pygcmc.DrudeComplete.clear()
        
        # Particle 1 - normal
        p1 = pygcmc.DrudeParticle()
        p1.drudeIndex = 1
        p1.parentIndex = 0
        p1.charge = -1.0
        p1.polarizability = 0.001
        p1.computeSpringConstants()
        pygcmc.DrudeComplete.addParticle(p1)
        
        # Particle 2 - varying polarizability
        p2 = pygcmc.DrudeParticle()
        p2.drudeIndex = 3
        p2.parentIndex = 2
        p2.charge = 1.0
        p2.polarizability = alpha2
        p2.computeSpringConstants()
        
        # Check spring constant scales inversely with alpha
        # k ~ q^2/alpha, so should be large for small alpha
        if alpha2 < 1e-8:
            assert p2.kSpring > 1e12, \
                f"Spring constant {p2.kSpring} too small for tiny alpha={alpha2}"
        
        pygcmc.DrudeComplete.addParticle(p2)
        
        # SCF parameters
        params = pygcmc.DrudeSCFParams()
        params.tolerance = 1e-6
        params.maxIterations = 100
        params.maxDrudeDistance = 0.02
        pygcmc.DrudeComplete.setParameters(params)
        
        # Should still calculate energy
        energy = pygcmc.DrudeComplete.calculateEnergy(state)
        assert np.isfinite(energy), f"Energy not finite for alpha={alpha2}"
        
        # Drude 2 should barely move with tiny polarizability
        dx = state.atoms[3].x - state.atoms[2].x
        dy = state.atoms[3].y - state.atoms[2].y
        dz = state.atoms[3].z - state.atoms[2].z
        displacement = math.sqrt(dx*dx + dy*dy + dz*dz)
        
        if alpha2 < 1e-8:
            assert displacement < 1e-8, \
                f"Zero-polarizability Drude moved: {displacement} nm"
    
    pygcmc.DrudeComplete.clear()


def test_many_body_polarization():
    """Test many-body polarization effects in dense system"""
    state = pygcmc.MCState()
    
    # Create a small cubic lattice of polarizable ions
    n_side = 3  # 3x3x3 = 27 ions
    spacing = 0.25  # nm
    
    atoms = []
    charge_pattern = [1.0, -1.0]  # Alternating charges
    
    idx = 0
    for i in range(n_side):
        for j in range(n_side):
            for k in range(n_side):
                x = i * spacing
                y = j * spacing
                z = k * spacing
                
                # Parent (permanent charge)
                parent = pygcmc.MCAtom()
                parent.x, parent.y, parent.z = x, y, z
                parent_charge = charge_pattern[idx % 2]
                parent.charge = parent_charge
                parent.type = 0
                atoms.append(parent)
                
                # Drude
                # Use a smaller magnitude Drude charge so each ion has a net charge and
                # a non-zero local electric field (avoids the trivial "neutral pairs at
                # same position => zero field" fixed point).
                drude = pygcmc.MCAtom()
                drude.x, drude.y, drude.z = x, y, z
                drude.charge = -0.2 * parent_charge
                drude.type = 1
                atoms.append(drude)
                
                idx += 1
    
    state.atoms = atoms
    state.activeAtomCount = len(atoms)
    state.info.box = [2.0, 2.0, 2.0]
    state.info.cutoff = 1.0
    
    # Setup residues
    residues = []
    n_ions = n_side ** 3
    for i in range(n_ions):
        res = pygcmc.MCResidue()
        res.atomStart = i * 2
        res.atomCount = 2
        res.active = True
        res.type = 0
        residues.append(res)
    
    state.residues = residues
    state.activeResidueCount = n_ions
    
    # Setup all Drude particles
    pygcmc.DrudeComplete.clear()
    
    for i in range(n_ions):
        p = pygcmc.DrudeParticle()
        p.drudeIndex = i * 2 + 1
        p.parentIndex = i * 2
        p.charge = -0.2 * charge_pattern[i % 2]
        p.polarizability = 0.001
        p.computeSpringConstants()
        pygcmc.DrudeComplete.addParticle(p)
    
    # Add Thole screening between all pairs
    for i in range(n_ions):
        for j in range(i+1, n_ions):
            pair = pygcmc.ScreenedPair()
            pair.dipole1 = i
            pair.dipole2 = j
            pair.thole = 1.3
            pygcmc.DrudeComplete.addScreenedPair(pair)
    
    # SCF parameters - may need many iterations for many-body
    params = pygcmc.DrudeSCFParams()
    params.tolerance = 0.1  # Looser tolerance for many-body
    params.maxIterations = 1000  # Many iterations
    params.maxDrudeDistance = 0.02
    params.dampingFactor = 0.5  # Conservative damping
    pygcmc.DrudeComplete.setParameters(params)
    
    # Calculate energy
    energy = pygcmc.DrudeComplete.calculateEnergy(state)
    
    # Check that energy is reasonable
    energy_per_ion = energy / n_ions
    assert -1000 < energy_per_ion < 1000, \
        f"Energy per ion {energy_per_ion} kJ/mol unreasonable"
    
    # Check that at least some Drudes have moved
    max_displacement = 0.0
    for i in range(n_ions):
        drude_idx = i * 2 + 1
        parent_idx = i * 2
        dx = state.atoms[drude_idx].x - state.atoms[parent_idx].x
        dy = state.atoms[drude_idx].y - state.atoms[parent_idx].y
        dz = state.atoms[drude_idx].z - state.atoms[parent_idx].z
        displacement = math.sqrt(dx*dx + dy*dy + dz*dz)
        max_displacement = max(max_displacement, displacement)
    
    assert max_displacement > 1e-6, \
        "No Drude displacement in many-body system"
    
    pygcmc.DrudeComplete.clear()


def test_pathological_geometry():
    """Test with pathological molecular geometry"""
    state = pygcmc.MCState()
    
    # Create a molecule with atoms in a straight line
    # This can cause numerical issues with polarization
    atoms = []
    
    # Five atoms in a line, alternating charges
    positions = [0.0, 0.1, 0.2, 0.3, 0.4]  # nm
    charges = [1.0, -2.0, 3.0, -2.0, 1.0]
    
    for i, (x, q) in enumerate(zip(positions, charges)):
        # Parent
        parent = pygcmc.MCAtom()
        parent.x, parent.y, parent.z = x, 0.0, 0.0
        parent.charge = q
        parent.type = 0
        atoms.append(parent)
        
        # Drude only on middle three atoms
        if 0 < i < 4:
            drude = pygcmc.MCAtom()
            drude.x, drude.y, drude.z = x, 0.0, 0.0
            drude.charge = -q
            drude.type = 1
            atoms.append(drude)
    
    state.atoms = atoms
    state.activeAtomCount = len(atoms)
    state.info.box = [3.0, 3.0, 3.0]
    
    # All atoms in one residue
    res = pygcmc.MCResidue()
    res.atomStart = 0
    res.atomCount = len(atoms)
    res.active = True
    res.type = 0
    state.residues = [res]
    state.activeResidueCount = 1
    
    # Setup Drude particles
    pygcmc.DrudeComplete.clear()
    
    drude_atom_indices = [2, 4, 6]  # Indices of Drude atoms
    parent_indices = [1, 3, 5]  # Their parents
    drude_particle_idx = 0
    
    for drude_idx, parent_idx in zip(drude_atom_indices, parent_indices):
        p = pygcmc.DrudeParticle()
        p.drudeIndex = drude_idx
        p.parentIndex = parent_idx
        p.charge = state.atoms[drude_idx].charge
        p.polarizability = 0.001
        p.computeSpringConstants()
        pygcmc.DrudeComplete.addParticle(p)
        drude_particle_idx += 1
    
    # SCF with high damping for stability
    params = pygcmc.DrudeSCFParams()
    params.tolerance = 1e-4
    params.maxIterations = 500
    params.maxDrudeDistance = 0.01
    params.dampingFactor = 0.8  # High damping
    pygcmc.DrudeComplete.setParameters(params)
    
    # Should handle pathological geometry
    energy = pygcmc.DrudeComplete.calculateEnergy(state)
    assert np.isfinite(energy), "Energy not finite for linear molecule"
    
    # Check that middle Drude (highest charge) moves most
    displacements = []
    for drude_idx, parent_idx in zip(drude_atom_indices, parent_indices):
        dx = state.atoms[drude_idx].x - state.atoms[parent_idx].x
        dy = state.atoms[drude_idx].y - state.atoms[parent_idx].y
        dz = state.atoms[drude_idx].z - state.atoms[parent_idx].z
        displacement = math.sqrt(dx*dx + dy*dy + dz*dz)
        displacements.append(displacement)
    
    # Middle atom has highest charge magnitude
    assert displacements[1] >= displacements[0] and displacements[1] >= displacements[2], \
        f"Middle Drude should move most: {displacements}"
    
    pygcmc.DrudeComplete.clear()

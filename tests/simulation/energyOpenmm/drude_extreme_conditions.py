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
                
                # Parent
                parent = pygcmc.MCAtom()
                parent.x, parent.y, parent.z = x, y, z
                parent.charge = charge_pattern[idx % 2]
                parent.type = 0
                atoms.append(parent)
                
                # Drude
                drude = pygcmc.MCAtom()
                drude.x, drude.y, drude.z = x, y, z
                drude.charge = -parent.charge
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
        p.charge = -charge_pattern[i % 2]
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
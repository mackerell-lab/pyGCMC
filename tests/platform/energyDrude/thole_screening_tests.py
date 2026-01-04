#!/usr/bin/env python3
"""
Thole screening functionality tests
Tests the damping of dipole-dipole interactions
"""

import pytest
import numpy as np
import pygcmc


def test_thole_screening_effect():
    """Test that Thole screening reduces dipole-dipole interaction"""
    print("\n=== Test: Thole Screening Effect ===")
    
    state = pygcmc.MCState()
    state.info.box = [5.0, 5.0, 5.0]
    
    # Two parent-Drude pairs
    atoms = []
    
    # Pair 1
    parent1 = pygcmc.MCAtom()
    parent1.x = parent1.y = parent1.z = 0.0
    parent1.charge = 1.0
    parent1.type = 0
    atoms.append(parent1)
    
    drude1 = pygcmc.MCAtom()
    drude1.x = drude1.y = drude1.z = 0.0
    drude1.charge = -1.0
    drude1.type = 1
    atoms.append(drude1)
    
    # Pair 2
    parent2 = pygcmc.MCAtom()
    parent2.x = 0.3  # 0.3 nm away
    parent2.y = parent2.z = 0.0
    parent2.charge = 1.0
    parent2.type = 0
    atoms.append(parent2)
    
    drude2 = pygcmc.MCAtom()
    drude2.x = 0.3
    drude2.y = drude2.z = 0.0
    drude2.charge = -1.0
    drude2.type = 1
    atoms.append(drude2)

    # Add a small external charge near dipole 1 to break symmetry and induce polarization.
    # This avoids relying on unphysical "spontaneous polarization" in a perfectly symmetric setup.
    external = pygcmc.MCAtom()
    external.x = -0.2
    external.y = external.z = 0.0
    external.charge = 0.4
    external.type = 0
    atoms.append(external)
    
    state.atoms = atoms
    state.activeAtomCount = 5
    
    # Three residues (two dipoles + external)
    residues = []
    for i in range(2):
        res = pygcmc.MCResidue()
        res.atomStart = i * 2
        res.atomCount = 2
        res.active = True
        res.type = 0
        residues.append(res)

    res_ext = pygcmc.MCResidue()
    res_ext.atomStart = 4
    res_ext.atomCount = 1
    res_ext.active = True
    res_ext.type = 1
    residues.append(res_ext)
    
    state.residues = residues
    state.activeResidueCount = 3
    
    # Setup Drude particles
    pygcmc.DrudeComplete.clear()
    
    for i in range(2):
        particle = pygcmc.DrudeParticle()
        particle.drudeIndex = i * 2 + 1
        particle.parentIndex = i * 2
        particle.charge = -1.0
        particle.polarizability = 0.001
        particle.aniso1Index = particle.aniso2Index = particle.aniso3Index = particle.aniso4Index = -1
        particle.aniso12 = particle.aniso34 = 1.0
        particle.computeSpringConstants()
        pygcmc.DrudeComplete.addParticle(particle)
    
    # SCF parameters
    params = pygcmc.DrudeSCFParams()
    params.tolerance = 0.01
    params.maxIterations = 100
    params.enableHardWall = False
    params.dampingFactor = 0.5
    params.includeCoulombEnergy = True
    pygcmc.DrudeComplete.setParameters(params)
    
    # Calculate WITHOUT Thole screening
    energy_no_thole = pygcmc.DrudeComplete.calculateEnergy(state)
    
    # Save displacements
    disp1_no_thole = np.sqrt((state.atoms[1].x - state.atoms[0].x)**2 + 
                             (state.atoms[1].y - state.atoms[0].y)**2 + 
                             (state.atoms[1].z - state.atoms[0].z)**2)
    disp2_no_thole = np.sqrt((state.atoms[3].x - state.atoms[2].x)**2 + 
                             (state.atoms[3].y - state.atoms[2].y)**2 + 
                             (state.atoms[3].z - state.atoms[2].z)**2)
    
    # Reset Drude positions
    for i in [1, 3]:  # Drude indices
        state.atoms[i].x = state.atoms[i-1].x
        state.atoms[i].y = state.atoms[i-1].y
        state.atoms[i].z = state.atoms[i-1].z
    
    # Add Thole screening
    screened = pygcmc.ScreenedPair()
    screened.dipole1 = 0  # First Drude particle in DrudeComplete
    screened.dipole2 = 1  # Second Drude particle
    screened.thole = 1.3  # Typical value
    pygcmc.DrudeComplete.addScreenedPair(screened)
    
    # Calculate WITH Thole screening
    energy_with_thole = pygcmc.DrudeComplete.calculateEnergy(state)
    
    # Save displacements
    disp1_with_thole = np.sqrt((state.atoms[1].x - state.atoms[0].x)**2 + 
                               (state.atoms[1].y - state.atoms[0].y)**2 + 
                               (state.atoms[1].z - state.atoms[0].z)**2)
    disp2_with_thole = np.sqrt((state.atoms[3].x - state.atoms[2].x)**2 + 
                               (state.atoms[3].y - state.atoms[2].y)**2 + 
                               (state.atoms[3].z - state.atoms[2].z)**2)
    
    print(f"\nWithout Thole screening:")
    print(f"  Energy: {energy_no_thole:.6f} kJ/mol")
    print(f"  Drude 1 displacement: {disp1_no_thole:.6f} nm")
    print(f"  Drude 2 displacement: {disp2_no_thole:.6f} nm")
    
    print(f"\nWith Thole screening (a = {screened.thole}):")
    print(f"  Energy: {energy_with_thole:.6f} kJ/mol")
    print(f"  Drude 1 displacement: {disp1_with_thole:.6f} nm")
    print(f"  Drude 2 displacement: {disp2_with_thole:.6f} nm")
    
    # Thole screening reduces dipole-dipole coupling, so the *remote* dipole response
    # (dipole 2) should be reduced when screening is enabled.
    assert disp2_no_thole > 1e-6, "Dipole 2 should respond to induced polarization"
    assert disp2_with_thole < disp2_no_thole, "Thole screening should reduce induced coupling"
    
    pygcmc.DrudeComplete.clear()


def test_thole_parameter_range():
    """Test different Thole parameter values"""
    print("\n=== Test: Thole Parameter Range ===")
    
    # Same system as before
    state = pygcmc.MCState()
    state.info.box = [5.0, 5.0, 5.0]
    
    # Two dipoles at 0.4 nm distance
    atoms = []
    for i in range(2):
        x = i * 0.4
        # Parent
        parent = pygcmc.MCAtom()
        parent.x = x
        parent.y = parent.z = 0.0
        parent.charge = 0.5
        parent.type = 0
        atoms.append(parent)
        
        # Drude
        drude = pygcmc.MCAtom()
        drude.x = x
        drude.y = drude.z = 0.0
        drude.charge = -1.0
        drude.type = 1
        atoms.append(drude)
    
    state.atoms = atoms
    state.activeAtomCount = 4
    
    # Residues
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
    
    # Test different Thole parameters
    thole_values = [0.0, 0.5, 1.0, 1.3, 2.0, 3.0]
    results = []
    
    for thole in thole_values:
        # Setup Drude particles
        pygcmc.DrudeComplete.clear()
        
        for i in range(2):
            particle = pygcmc.DrudeParticle()
            particle.drudeIndex = i * 2 + 1
            particle.parentIndex = i * 2
            particle.charge = -1.0
            particle.polarizability = 0.001
            particle.aniso1Index = particle.aniso2Index = particle.aniso3Index = particle.aniso4Index = -1
            particle.aniso12 = particle.aniso34 = 1.0
            particle.computeSpringConstants()
            pygcmc.DrudeComplete.addParticle(particle)
        
        # Add screening if thole > 0
        if thole > 0:
            screened = pygcmc.ScreenedPair()
            screened.dipole1 = 0
            screened.dipole2 = 1
            screened.thole = thole
            pygcmc.DrudeComplete.addScreenedPair(screened)
        
        # Set parameters
        params = pygcmc.DrudeSCFParams()
        params.tolerance = 0.01
        params.maxIterations = 100
        params.enableHardWall = False
        params.dampingFactor = 0.5
        pygcmc.DrudeComplete.setParameters(params)
        
        # Reset positions
        for i in range(4):
            if i % 2 == 1:  # Drude
                state.atoms[i].x = state.atoms[i-1].x
                state.atoms[i].y = state.atoms[i-1].y
                state.atoms[i].z = state.atoms[i-1].z
        
        # Calculate
        energy = pygcmc.DrudeComplete.calculateEnergy(state)
        
        # Get average displacement
        disp1 = np.sqrt((state.atoms[1].x - state.atoms[0].x)**2 + 
                        (state.atoms[1].y - state.atoms[0].y)**2 + 
                        (state.atoms[1].z - state.atoms[0].z)**2)
        disp2 = np.sqrt((state.atoms[3].x - state.atoms[2].x)**2 + 
                        (state.atoms[3].y - state.atoms[2].y)**2 + 
                        (state.atoms[3].z - state.atoms[2].z)**2)
        avg_disp = (disp1 + disp2) / 2
        
        results.append((thole, energy, avg_disp))
    
    print(f"\n{'Thole':>6} {'Energy':>12} {'Avg Disp':>12}")
    print("-" * 32)
    
    for thole, energy, disp in results:
        print(f"{thole:6.1f} {energy:12.6f} {disp:12.6f}")
    
    # Check that Thole has an effect
    # The energy should change monotonically with Thole parameter
    # But displacement may not be monotonic due to complex interactions
    no_thole_disp = results[0][2]  # Displacement with no Thole
    
    # At least some Thole values should significantly change displacement
    significant_changes = 0
    # The exact displacement change depends on the chosen parameters and SCF tolerance,
    # but Thole screening must produce a detectable effect above numerical noise.
    min_disp_change = max(2e-5, 0.01 * no_thole_disp)
    for i in range(1, len(results)):
        if abs(results[i][2] - no_thole_disp) > min_disp_change:
            significant_changes += 1
    
    assert significant_changes >= 2, "Thole should affect displacement for multiple values"
    
    pygcmc.DrudeComplete.clear()


def test_thole_distance_dependence():
    """Test that Thole screening effect depends on distance"""
    print("\n=== Test: Thole Distance Dependence ===")
    
    # Test at different distances
    distances = [0.2, 0.3, 0.4, 0.5, 0.7, 1.0]
    thole = 1.3
    
    screening_effects = []
    
    for dist in distances:
        state = pygcmc.MCState()
        state.info.box = [5.0, 5.0, 5.0]
        
        # Two dipoles
        atoms = []
        
        # Dipole 1
        parent1 = pygcmc.MCAtom()
        parent1.x = parent1.y = parent1.z = 0.0
        parent1.charge = 0.5
        parent1.type = 0
        atoms.append(parent1)
        
        drude1 = pygcmc.MCAtom()
        drude1.x = drude1.y = drude1.z = 0.0
        drude1.charge = -1.0
        drude1.type = 1
        atoms.append(drude1)
        
        # Dipole 2
        parent2 = pygcmc.MCAtom()
        parent2.x = dist
        parent2.y = parent2.z = 0.0
        parent2.charge = 0.5
        parent2.type = 0
        atoms.append(parent2)
        
        drude2 = pygcmc.MCAtom()
        drude2.x = dist
        drude2.y = drude2.z = 0.0
        drude2.charge = -1.0
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
            res.type = 0
            residues.append(res)
        
        state.residues = residues
        state.activeResidueCount = 2
        
        # Calculate without Thole
        pygcmc.DrudeComplete.clear()
        
        for i in range(2):
            particle = pygcmc.DrudeParticle()
            particle.drudeIndex = i * 2 + 1
            particle.parentIndex = i * 2
            particle.charge = -1.0
            particle.polarizability = 0.001
            particle.aniso1Index = particle.aniso2Index = particle.aniso3Index = particle.aniso4Index = -1
            particle.aniso12 = particle.aniso34 = 1.0
            particle.computeSpringConstants()
            pygcmc.DrudeComplete.addParticle(particle)
        
        params = pygcmc.DrudeSCFParams()
        params.tolerance = 0.01
        params.maxIterations = 100
        params.enableHardWall = False
        params.dampingFactor = 0.5
        pygcmc.DrudeComplete.setParameters(params)
        
        energy_no_thole = pygcmc.DrudeComplete.calculateEnergy(state)
        disp_no_thole = state.atoms[1].x - state.atoms[0].x
        
        # Reset and calculate with Thole
        state.atoms[1].x = state.atoms[1].y = state.atoms[1].z = 0.0
        state.atoms[3].x = dist
        state.atoms[3].y = state.atoms[3].z = 0.0
        
        # Add Thole screening
        screened = pygcmc.ScreenedPair()
        screened.dipole1 = 0
        screened.dipole2 = 1
        screened.thole = thole
        pygcmc.DrudeComplete.addScreenedPair(screened)
        
        energy_with_thole = pygcmc.DrudeComplete.calculateEnergy(state)
        disp_with_thole = state.atoms[1].x - state.atoms[0].x
        
        # Calculate screening effect
        if abs(disp_no_thole) > 1e-8:
            effect = 1.0 - abs(disp_with_thole / disp_no_thole)
        else:
            effect = 0.0
        
        screening_effects.append((dist, effect))
    
    print(f"\n{'Distance':>10} {'Screening Effect':>20}")
    print("-" * 32)
    
    for dist, effect in screening_effects:
        print(f"{dist:10.1f} {effect:20.3f}")
    
    # Check that screening effect varies with distance
    # The effect should be different at different distances
    effects = [e[1] for e in screening_effects]
    effect_range = max(effects) - min(effects)
    
    assert effect_range > 0.01, "Screening effect should vary with distance"
    
    # At very long distances, the effect should be different from short distances
    short_dist_effect = abs(screening_effects[0][1])
    long_dist_effect = abs(screening_effects[-1][1])
    assert abs(short_dist_effect - long_dist_effect) > 0.01, \
        "Screening effect should differ between short and long distances"
    
    pygcmc.DrudeComplete.clear()


if __name__ == "__main__":
    test_thole_screening_effect()
    test_thole_parameter_range()
    test_thole_distance_dependence()

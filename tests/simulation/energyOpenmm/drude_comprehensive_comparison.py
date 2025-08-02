#!/usr/bin/env python3
"""
Comprehensive Drude-OpenMM comparison tests
Validates numerical consistency between PyGCMC and OpenMM implementations
"""

import pytest
import numpy as np
try:
    import openmm
    import openmm.unit as unit
    HAS_OPENMM = True
except ImportError:
    HAS_OPENMM = False

import pygcmc


@pytest.mark.skipif(not HAS_OPENMM, reason="OpenMM not available")
def test_single_drude_external_field():
    """Test single Drude in external field comparing PyGCMC vs OpenMM"""
    print("\n=== Single Drude in External Field ===")
    
    # System setup
    positions = [
        [0.0, 0.0, 0.0],   # Parent (neutral)
        [0.0, 0.0, 0.0],   # Drude
        [0.5, 0.0, 0.0],   # External charge
    ]
    
    # PyGCMC setup
    state = pygcmc.MCState()
    state.info.box = [3.0, 3.0, 3.0]
    
    # Create atoms
    parent = pygcmc.MCAtom()
    parent.x = parent.y = parent.z = 0.0
    parent.charge = 0.0
    parent.type = 0
    
    drude = pygcmc.MCAtom()
    drude.x = drude.y = drude.z = 0.0
    drude.charge = -1.0
    drude.type = 1
    
    external = pygcmc.MCAtom()
    external.x = 0.5
    external.y = external.z = 0.0
    external.charge = 2.0
    external.type = 2
    
    state.atoms = [parent, drude, external]
    state.activeAtomCount = 3
    
    # Setup residues
    res1 = pygcmc.MCResidue()
    res1.atomStart = 0
    res1.atomCount = 2
    res1.active = True
    res1.type = 0
    state.residues = [res1]
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
    pygcmc.DrudeComplete.addParticle(particle)
    
    # Set tight tolerance for comparison
    params = pygcmc.DrudeSCFParams()
    params.tolerance = 1e-5  # Match OpenMM
    params.maxIterations = 1000
    params.enableHardWall = False
    params.dampingFactor = 0.5
    pygcmc.DrudeComplete.setParameters(params)
    
    # Calculate with PyGCMC
    energy_pygcmc = pygcmc.DrudeComplete.calculateEnergy(state)
    drude_disp_pygcmc = state.atoms[1].x - state.atoms[0].x
    
    # OpenMM setup
    system = openmm.System()
    system.addParticle(16.0 * unit.dalton)  # Parent
    system.addParticle(0.4 * unit.dalton)   # Drude
    system.addParticle(1.0 * unit.dalton)   # External
    
    # Drude force
    drude_force = openmm.DrudeForce()
    drude_force.addParticle(1, 0, -1, -1, -1, -1.0, 0.001, 1.0, 1.0)
    system.addForce(drude_force)
    
    # Nonbonded force
    nonbonded = openmm.NonbondedForce()
    nonbonded.setNonbondedMethod(openmm.NonbondedForce.NoCutoff)
    nonbonded.addParticle(0.0, 0.1 * unit.nanometer, 0.0)
    nonbonded.addParticle(-1.0 * unit.elementary_charge, 0.1 * unit.nanometer, 0.0)
    nonbonded.addParticle(2.0 * unit.elementary_charge, 0.1 * unit.nanometer, 0.0)
    nonbonded.addException(0, 1, 0.0, 1.0, 0.0)  # Exclude parent-drude
    system.addForce(nonbonded)
    
    # Create context
    integrator = openmm.DrudeSCFIntegrator(0.001 * unit.picoseconds)
    integrator.setMinimizationErrorTolerance(1e-5)
    context = openmm.Context(system, integrator)
    context.setPositions(positions * unit.nanometer)
    
    # Optimize
    integrator.step(1)
    
    # Get results
    state_omm = context.getState(getPositions=True, getEnergy=True)
    pos_omm = state_omm.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
    energy_openmm = state_omm.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    drude_disp_openmm = pos_omm[1][0] - pos_omm[0][0]
    
    print(f"\nDrude displacement:")
    print(f"  PyGCMC: {drude_disp_pygcmc:.8f} nm")
    print(f"  OpenMM: {drude_disp_openmm:.8f} nm")
    print(f"  Ratio: {drude_disp_pygcmc/drude_disp_openmm:.6f}")
    
    # We expect ~0.5% difference as documented
    assert 0.990 < drude_disp_pygcmc/drude_disp_openmm < 1.010, \
        "Displacement ratio outside expected 1% tolerance"
    
    pygcmc.DrudeComplete.clear()


@pytest.mark.skipif(not HAS_OPENMM, reason="OpenMM not available")
def test_water_dimer_thole_screening():
    """Test water dimer with Thole screening"""
    print("\n=== Water Dimer with Thole Screening ===")
    
    # Two water-like molecules
    positions = [
        [0.0, 0.0, 0.0],   # O1
        [0.0, 0.0, 0.0],   # D1
        [0.3, 0.0, 0.0],   # O2
        [0.3, 0.0, 0.0],   # D2
    ]
    
    # PyGCMC setup
    state = pygcmc.MCState()
    state.info.box = [3.0, 3.0, 3.0]
    
    # Create atoms
    atoms = []
    charges = [-0.662, -1.338, -0.662, -1.338]  # SWM4-NDP-like
    
    for i, (pos, charge) in enumerate(zip(positions, charges)):
        atom = pygcmc.MCAtom()
        atom.x, atom.y, atom.z = pos
        atom.charge = charge
        atom.type = i % 2
        atoms.append(atom)
    
    state.atoms = atoms
    state.activeAtomCount = 4
    
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
    
    # Setup Drude particles
    pygcmc.DrudeComplete.clear()
    
    # Water 1
    particle1 = pygcmc.DrudeParticle()
    particle1.drudeIndex = 1
    particle1.parentIndex = 0
    particle1.charge = charges[1]
    particle1.polarizability = 0.00097822  # SWM4-NDP
    particle1.aniso1Index = particle1.aniso2Index = particle1.aniso3Index = particle1.aniso4Index = -1
    particle1.aniso12 = particle1.aniso34 = 1.0
    particle1.computeSpringConstants()
    pygcmc.DrudeComplete.addParticle(particle1)
    
    # Water 2
    particle2 = pygcmc.DrudeParticle()
    particle2.drudeIndex = 3
    particle2.parentIndex = 2
    particle2.charge = charges[3]
    particle2.polarizability = 0.00097822
    particle2.aniso1Index = particle2.aniso2Index = particle2.aniso3Index = particle2.aniso4Index = -1
    particle2.aniso12 = particle2.aniso34 = 1.0
    particle2.computeSpringConstants()
    pygcmc.DrudeComplete.addParticle(particle2)
    
    # Add Thole screening
    screened = pygcmc.ScreenedPair()
    screened.dipole1 = 0  # First Drude particle
    screened.dipole2 = 1  # Second Drude particle
    screened.thole = 1.3  # SWM4-NDP value
    pygcmc.DrudeComplete.addScreenedPair(screened)
    
    # Set parameters
    params = pygcmc.DrudeSCFParams()
    params.tolerance = 0.01
    params.maxIterations = 100
    params.enableHardWall = False
    params.dampingFactor = 0.5
    pygcmc.DrudeComplete.setParameters(params)
    
    # Calculate
    energy = pygcmc.DrudeComplete.calculateEnergy(state)
    
    # Check that Drudes have moved
    disp1 = np.sqrt((state.atoms[1].x - state.atoms[0].x)**2 + 
                    (state.atoms[1].y - state.atoms[0].y)**2 + 
                    (state.atoms[1].z - state.atoms[0].z)**2)
    disp2 = np.sqrt((state.atoms[3].x - state.atoms[2].x)**2 + 
                    (state.atoms[3].y - state.atoms[2].y)**2 + 
                    (state.atoms[3].z - state.atoms[2].z)**2)
    
    print(f"\nDrude displacements:")
    print(f"  Water 1: {disp1:.6f} nm")
    print(f"  Water 2: {disp2:.6f} nm")
    print(f"  Energy: {energy:.3f} kJ/mol")
    
    # Both should have non-zero displacement due to mutual polarization
    assert disp1 > 1e-6, "First Drude should move"
    assert disp2 > 1e-6, "Second Drude should move"
    
    pygcmc.DrudeComplete.clear()


@pytest.mark.skipif(not HAS_OPENMM, reason="OpenMM not available")
def test_convergence_tolerance_effect():
    """Test effect of convergence tolerance on results"""
    print("\n=== Convergence Tolerance Effect ===")
    
    # Simple test system
    state = pygcmc.MCState()
    state.info.box = [3.0, 3.0, 3.0]
    
    parent = pygcmc.MCAtom()
    parent.x = parent.y = parent.z = 0.0
    parent.charge = 0.0
    parent.type = 0
    
    drude = pygcmc.MCAtom()
    drude.x = drude.y = drude.z = 0.0
    drude.charge = -1.0
    drude.type = 1
    
    external = pygcmc.MCAtom()
    external.x = 0.5
    external.y = external.z = 0.0
    external.charge = 2.0
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
    
    # Test different tolerances
    tolerances = [0.1, 0.01, 0.001, 0.0001]
    results = []
    
    for tol in tolerances:
        pygcmc.DrudeComplete.clear()
        
        particle = pygcmc.DrudeParticle()
        particle.drudeIndex = 1
        particle.parentIndex = 0
        particle.charge = -1.0
        particle.polarizability = 0.001
        particle.aniso1Index = particle.aniso2Index = particle.aniso3Index = particle.aniso4Index = -1
        particle.aniso12 = particle.aniso34 = 1.0
        particle.computeSpringConstants()
        pygcmc.DrudeComplete.addParticle(particle)
        
        params = pygcmc.DrudeSCFParams()
        params.tolerance = tol
        params.maxIterations = 1000
        params.enableHardWall = False
        params.dampingFactor = 0.5
        pygcmc.DrudeComplete.setParameters(params)
        
        # Reset Drude position
        state.atoms[1].x = state.atoms[1].y = state.atoms[1].z = 0.0
        
        energy = pygcmc.DrudeComplete.calculateEnergy(state)
        displacement = state.atoms[1].x
        
        results.append((tol, displacement, energy))
        
    print(f"\n{'Tolerance':>10} {'Displacement':>12} {'Energy':>12}")
    print("-" * 36)
    
    for tol, disp, energy in results:
        print(f"{tol:10.4f} {disp:12.8f} {energy:12.6f}")
    
    # Check convergence
    # Tighter tolerance should give more consistent results
    disp_tight = results[-1][1]  # Tightest tolerance
    for i in range(len(results)-1):
        tol, disp, _ = results[i]
        rel_diff = abs(disp - disp_tight) / disp_tight
        print(f"\nTol {tol}: {rel_diff*100:.2f}% difference from tightest")
        
        # Looser tolerances can have up to 1% difference
        if tol >= 0.01:
            assert rel_diff < 0.01, f"Tolerance {tol} gives too large error"
    
    pygcmc.DrudeComplete.clear()


if __name__ == "__main__":
    test_single_drude_external_field()
    test_water_dimer_thole_screening()
    test_convergence_tolerance_effect()
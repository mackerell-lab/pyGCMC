"""
SCF consistency tests inspired by OpenMM's TestDrudeSCFIntegrator
Note: This tests SCF repeat consistency, not true energy conservation in MD
"""

import pytest
import numpy as np
import pygcmc
import math


def test_scf_repeat_consistency():
    """Test SCF repeat consistency (inspired by OpenMM's energy conservation test)
    
    Note: OpenMM's actual test runs dynamics with DrudeSCFIntegrator for thousands
    of steps and checks total energy conservation. Here we only test that repeated
    SCF calculations give consistent results.
    """
    
    # Create a simple system with 2 Drude oscillators
    state = pygcmc.MCState()
    state.info.box = [3.0, 3.0, 3.0]
    state.info.cutoff = 1.5
    
    # Force field
    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 2
    ff.numMovementTypes = 2
    ff.ljEps = [0.0, 0.0]
    ff.ljSigma = [0.1, 0.1]
    state.forcefield = ff
    
    # Create two molecules, each with parent and Drude
    atoms = []
    
    # Molecule 1 - neutral parent with Drude
    parent1 = pygcmc.MCAtom()
    parent1.x, parent1.y, parent1.z = 0.0, 0.0, 0.0
    parent1.charge = 0.0  # Neutral
    parent1.type = 0
    atoms.append(parent1)
    
    drude1 = pygcmc.MCAtom()
    drude1.x, drude1.y, drude1.z = 0.0, 0.0, 0.0
    drude1.charge = -1.0
    drude1.type = 1
    atoms.append(drude1)
    
    # Molecule 2 - has net positive charge to create field
    parent2 = pygcmc.MCAtom()
    parent2.x, parent2.y, parent2.z = 1.0, 0.0, 0.0
    parent2.charge = 2.0  # Net +1 charge after Drude
    parent2.type = 0
    atoms.append(parent2)
    
    drude2 = pygcmc.MCAtom()
    drude2.x, drude2.y, drude2.z = 1.0, 0.0, 0.0
    drude2.charge = -1.0
    drude2.type = 1
    atoms.append(drude2)
    
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
    
    for i in range(2):
        particle = pygcmc.DrudeParticle()
        particle.drudeIndex = i * 2 + 1
        particle.parentIndex = i * 2
        particle.charge = -1.0
        particle.polarizability = 0.001  # nm^3
        particle.computeSpringConstants()
        pygcmc.DrudeComplete.addParticle(particle)
    
    # SCF parameters (matching OpenMM test)
    params = pygcmc.DrudeSCFParams()
    params.tolerance = 0.1  # 0.1 kJ/mol/nm like OpenMM
    params.maxIterations = 100
    params.enableHardWall = False  # Match OpenMM SCF
    pygcmc.DrudeComplete.setParameters(params)
    
    # Calculate initial energy
    initial_energy = pygcmc.DrudeComplete.calculateEnergy(state)
    
    # With two Drude oscillators and charge separation, they should induce dipoles
    print(f"Initial energy: {initial_energy:.6f} kJ/mol")
    
    # Check that Drude particles have moved from their parent positions
    # due to induced dipole interactions
    drude1_disp = math.sqrt((state.atoms[1].x - state.atoms[0].x)**2 + 
                           (state.atoms[1].y - state.atoms[0].y)**2 +
                           (state.atoms[1].z - state.atoms[0].z)**2)
    drude2_disp = math.sqrt((state.atoms[3].x - state.atoms[2].x)**2 + 
                           (state.atoms[3].y - state.atoms[2].y)**2 +
                           (state.atoms[3].z - state.atoms[2].z)**2)
    
    print(f"Drude 1 displacement: {drude1_disp:.6f} nm")
    print(f"Drude 2 displacement: {drude2_disp:.6f} nm")
    
    # Energy should be negative due to dipole-dipole attraction
    assert initial_energy < 0, f"Energy should be negative for induced dipoles, got {initial_energy}"
    
    # Test that energy is consistent after multiple calculations
    for _ in range(5):
        energy = pygcmc.DrudeComplete.calculateEnergy(state)
        assert abs(energy - initial_energy) < 1e-6, "Energy not consistent"
    
    pygcmc.DrudeComplete.clear()
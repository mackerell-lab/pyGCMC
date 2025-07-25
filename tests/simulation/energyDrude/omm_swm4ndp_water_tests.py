"""
SWM4-NDP water system test inspired by OpenMM's testWater
"""

import pytest
import numpy as np
import pygcmc
import math


def test_swm4_ndp_water_system():
    """Test SWM4-NDP water system (simplified version of OpenMM's testWater)"""
    
    # Create a single SWM4-NDP water molecule
    state = pygcmc.MCState()
    state.info.box = [3.0, 3.0, 3.0]
    state.info.cutoff = 1.0
    
    # Force field
    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 4  # O, D, H, M
    ff.numMovementTypes = 4
    ff.ljEps = [0.21094*4.184, 0.0, 0.0, 0.0]  # Only O has LJ
    ff.ljSigma = [0.318395, 0.1, 0.1, 0.1]
    state.forcefield = ff
    
    # Create water molecule atoms
    atoms = []
    
    # Oxygen
    oxygen = pygcmc.MCAtom()
    oxygen.x, oxygen.y, oxygen.z = 0.0, 0.0, 0.0
    oxygen.charge = 1.71636
    oxygen.type = 0
    atoms.append(oxygen)
    
    # Drude
    drude = pygcmc.MCAtom()
    drude.x, drude.y, drude.z = 0.0, 0.0, 0.0
    drude.charge = -1.71636
    drude.type = 1
    atoms.append(drude)
    
    # Hydrogen 1
    h1 = pygcmc.MCAtom()
    h1.x = 0.09572
    h1.y = 0.0
    h1.z = 0.0
    h1.charge = 0.55733
    h1.type = 2
    atoms.append(h1)
    
    # Hydrogen 2
    angle = 104.52 * math.pi / 180.0
    h2 = pygcmc.MCAtom()
    h2.x = 0.09572 * math.cos(angle)
    h2.y = 0.09572 * math.sin(angle)
    h2.z = 0.0
    h2.charge = 0.55733
    h2.type = 2
    atoms.append(h2)
    
    # Virtual site M
    bisector_angle = angle / 2.0
    m_site = pygcmc.MCAtom()
    m_site.x = 0.024034 * math.cos(bisector_angle)
    m_site.y = 0.024034 * math.sin(bisector_angle)
    m_site.z = 0.0
    m_site.charge = -1.11466
    m_site.type = 3
    atoms.append(m_site)
    
    state.atoms = atoms
    state.activeAtomCount = 5
    
    # Single residue for water
    res = pygcmc.MCResidue()
    res.atomStart = 0
    res.atomCount = 5
    res.active = True
    res.type = 0
    state.residues = [res]
    state.activeResidueCount = 1
    
    # Setup Drude particle
    pygcmc.DrudeComplete.clear()
    
    # From OpenMM test: polarizability = ONE_4PI_EPS0*1.71636*1.71636/(100000*4.184)
    polarizability = 138.935456 * 1.71636 * 1.71636 / (100000 * 4.184)
    
    particle = pygcmc.DrudeParticle()
    particle.drudeIndex = 1
    particle.parentIndex = 0
    particle.charge = -1.71636
    particle.polarizability = polarizability
    particle.computeSpringConstants()
    pygcmc.DrudeComplete.addParticle(particle)
    
    # SCF parameters
    params = pygcmc.DrudeSCFParams()
    params.tolerance = 0.1  # Match OpenMM test
    params.maxIterations = 100
    params.enableHardWall = False
    pygcmc.DrudeComplete.setParameters(params)
    
    # Calculate energy
    energy = pygcmc.DrudeComplete.calculateEnergy(state)
    
    print(f"SWM4-NDP single water energy: {energy:.6f} kJ/mol")
    
    # For a single water molecule with no external field,
    # the energy should be close to zero (only self-polarization)
    # OpenMM test uses 0.1 kJ/mol tolerance
    assert abs(energy) < 0.1, f"Single water energy too large: {energy} kJ/mol (threshold: 0.1)"
    
    pygcmc.DrudeComplete.clear()
#!/usr/bin/env python
"""Test SCF convergence with different parameters"""

import sys
import numpy as np
import math
sys.path.append('/home/zhaomt/gcmc/test107/pygcmc_dev/build/modules/bindings')

import pygcmc

def create_water_dimer():
    """Create two water molecules for testing"""
    state = pygcmc.MCState()
    state.info.box = [10.0, 10.0, 10.0]
    state.info.cutoff = 5.0
    
    # Initialize Drude force
    pygcmc.initializeDrudeForce()
    
    atoms = []
    
    # SWM4-NDP parameters
    qO = 1.71636
    qD = -1.71636
    qH = 0.55733
    qM = -1.11466
    rOH = 0.09572  # nm
    aHOH = 104.52 * math.pi / 180  # radians
    
    # Create two waters separated by different distances
    water_positions = [(5.0, 5.0, 5.0), (5.3, 5.0, 5.0)]  # 0.3 nm O-O distance
    
    for water_idx, (x, y, z) in enumerate(water_positions):
        # Oxygen
        o = pygcmc.MCAtom()
        o.x, o.y, o.z = x, y, z
        o.charge = qO
        o.type = 0
        atoms.append(o)
        
        # Drude - start with small displacement
        d = pygcmc.MCAtom()
        d.x, d.y, d.z = x + 0.001, y, z  # Small initial displacement
        d.charge = qD
        d.type = 1
        atoms.append(d)
        
        # Hydrogen 1
        h1 = pygcmc.MCAtom()
        h1.x = x + rOH
        h1.y = y
        h1.z = z
        h1.charge = qH
        h1.type = 2
        atoms.append(h1)
        
        # Hydrogen 2
        h2 = pygcmc.MCAtom()
        h2.x = x + rOH * math.cos(aHOH)
        h2.y = y + rOH * math.sin(aHOH)
        h2.z = z
        h2.charge = qH
        h2.type = 2
        atoms.append(h2)
        
        # M-site
        m = pygcmc.MCAtom()
        w_O = 0.786646558
        w_H = 0.106676721
        m.x = w_O * o.x + w_H * h1.x + w_H * h2.x
        m.y = w_O * o.y + w_H * h1.y + w_H * h2.y
        m.z = w_O * o.z + w_H * h1.z + w_H * h2.z
        m.charge = qM
        m.type = 3
        atoms.append(m)
    
    state.atoms = atoms
    state.activeAtomCount = len(atoms)
    
    # Create residues
    residues = []
    for i in range(2):
        res = pygcmc.MCResidue()
        res.atomStart = 5 * i
        res.atomCount = 5
        res.active = True
        res.type = 0
        residues.append(res)
    
    state.residues = residues
    state.activeResidueCount = len(residues)
    
    # Add Drude particles
    for i in range(2):
        pygcmc.addDrudeParticle(
            drudeIndex=5*i + 1,
            parentIndex=5*i,
            charge=qD,
            polarizability=0.000978253  # nm^3
        )
    
    # Add Thole screening between waters
    pygcmc.addDrudeScreenedPair(0, 1, 1.3)
    
    return state

def test_scf_parameters():
    """Test different SCF parameter combinations"""
    print("=== Testing SCF Convergence Parameters ===\n")
    
    # Test different parameter combinations
    test_cases = [
        {"tolerance": 1.0, "max_iter": 50, "damping": 0.5, "max_dist": 0.02},
        {"tolerance": 10.0, "max_iter": 100, "damping": 0.5, "max_dist": 0.02},
        {"tolerance": 50.0, "max_iter": 200, "damping": 0.3, "max_dist": 0.05},
        {"tolerance": 100.0, "max_iter": 500, "damping": 0.2, "max_dist": 0.1},
    ]
    
    for idx, params in enumerate(test_cases):
        print(f"Test {idx + 1}: tolerance={params['tolerance']}, max_iter={params['max_iter']}, " +
              f"damping={params['damping']}, max_dist={params['max_dist']}")
        
        # Create fresh water dimer
        state = create_water_dimer()
        
        # Set SCF parameters
        scf_params = pygcmc.DrudeSCFParams()
        scf_params.tolerance = params['tolerance']
        scf_params.maxIterations = params['max_iter']
        scf_params.dampingFactor = params['damping']
        scf_params.maxDrudeDistance = params['max_dist']
        pygcmc.setDrudeSCFParameters(scf_params)
        
        # Calculate energy
        result = pygcmc.computeSystemEnergyDrude(state)
        
        if isinstance(result, tuple):
            energy, components = result
            print(f"  Energy: {energy:.2f} kJ/mol")
            print(f"  Components: {components}")
        else:
            energy = result
            print(f"  Energy: {energy:.2f} kJ/mol")
        
        # Check Drude displacements
        for i in range(2):
            dx = state.atoms[5*i + 1].x - state.atoms[5*i].x
            dy = state.atoms[5*i + 1].y - state.atoms[5*i].y
            dz = state.atoms[5*i + 1].z - state.atoms[5*i].z
            dist = math.sqrt(dx*dx + dy*dy + dz*dz)
            print(f"  Water {i+1} Drude displacement: {dist*1000:.3f} pm")
        
        print()

def test_initial_displacement():
    """Test effect of initial Drude displacement"""
    print("\n=== Testing Initial Drude Displacement ===\n")
    
    initial_displacements = [0.0, 0.0001, 0.001, 0.01]  # nm
    
    # Set relaxed SCF parameters
    scf_params = pygcmc.DrudeSCFParams()
    scf_params.tolerance = 100.0
    scf_params.maxIterations = 500
    scf_params.dampingFactor = 0.2
    scf_params.maxDrudeDistance = 0.1
    pygcmc.setDrudeSCFParameters(scf_params)
    
    for disp in initial_displacements:
        print(f"Initial displacement: {disp} nm")
        
        # Create water dimer with specified initial displacement
        state = create_water_dimer()
        
        # Adjust initial Drude positions
        for i in range(2):
            state.atoms[5*i + 1].x = state.atoms[5*i].x + disp
        
        # Calculate energy
        result = pygcmc.computeSystemEnergyDrude(state)
        
        if isinstance(result, tuple):
            energy, components = result
            print(f"  Energy: {energy:.2f} kJ/mol")
            
            # Check final Drude positions
            for i in range(2):
                dx = state.atoms[5*i + 1].x - state.atoms[5*i].x
                dy = state.atoms[5*i + 1].y - state.atoms[5*i].y
                dz = state.atoms[5*i + 1].z - state.atoms[5*i].z
                dist = math.sqrt(dx*dx + dy*dy + dz*dz)
                print(f"  Water {i+1} final displacement: {dist*1000:.3f} pm")
        
        print()

def test_system_size():
    """Test convergence for different system sizes"""
    print("\n=== Testing System Size Effect ===\n")
    
    # Set relaxed parameters
    scf_params = pygcmc.DrudeSCFParams()
    scf_params.tolerance = 100.0
    scf_params.maxIterations = 500
    scf_params.dampingFactor = 0.2
    scf_params.maxDrudeDistance = 0.1
    pygcmc.setDrudeSCFParameters(scf_params)
    
    system_sizes = [1, 2, 4, 8]  # Number of waters
    
    for n_waters in system_sizes:
        print(f"System with {n_waters} water(s)")
        
        # Create state
        state = pygcmc.MCState()
        box_size = (n_waters / 10.0) ** (1/3) * 3.0  # Approximate density
        state.info.box = [box_size, box_size, box_size]
        state.info.cutoff = min(box_size/2 - 0.1, 1.2)
        
        # Initialize Drude force
        pygcmc.initializeDrudeForce()
        
        # Create waters in a grid
        atoms = []
        residues = []
        
        # Simplified: just place waters randomly
        for i in range(n_waters):
            x = np.random.random() * box_size
            y = np.random.random() * box_size
            z = np.random.random() * box_size
            
            # Create water atoms (simplified from earlier code)
            # ... (similar to create_water_dimer but with random positions)
            
        # For now, just report the test setup
        print(f"  Box size: {box_size:.2f} nm")
        print(f"  Would test {n_waters} waters")
        print()

if __name__ == "__main__":
    test_scf_parameters()
    test_initial_displacement()
    test_system_size()
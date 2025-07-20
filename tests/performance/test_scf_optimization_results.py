#!/usr/bin/env python
"""Test SCF optimization results"""

import sys
import numpy as np
import math
import time
sys.path.append('/home/zhaomt/gcmc/test107/pygcmc_dev/build/modules/bindings')

import pygcmc

def create_water_system(n_waters, box_size):
    """Create water system with proper residues"""
    state = pygcmc.MCState()
    state.info.box = [box_size, box_size, box_size]
    state.info.cutoff = min(box_size/2 - 0.1, 1.2)
    
    # Initialize Drude force
    pygcmc.initializeDrudeForce()
    
    atoms = []
    residues = []
    
    # SWM4-NDP parameters
    qO = 1.71636
    qD = -1.71636
    qH = 0.55733
    qM = -1.11466
    rOH = 0.09572  # nm
    aHOH = 104.52 * math.pi / 180
    
    # Create waters in a grid
    n_per_dim = int(np.ceil(n_waters**(1/3)))
    spacing = box_size / n_per_dim
    
    water_id = 0
    for i in range(n_per_dim):
        for j in range(n_per_dim):
            for k in range(n_per_dim):
                if water_id >= n_waters:
                    break
                
                x = (i + 0.5) * spacing
                y = (j + 0.5) * spacing
                z = (k + 0.5) * spacing
                
                # Oxygen
                o = pygcmc.MCAtom()
                o.x, o.y, o.z = x, y, z
                o.charge = qO
                o.type = 0
                atoms.append(o)
                
                # Drude - start at parent
                d = pygcmc.MCAtom()
                d.x, d.y, d.z = x, y, z
                d.charge = qD
                d.type = 1
                atoms.append(d)
                
                # Hydrogens
                h1 = pygcmc.MCAtom()
                h1.x = x + rOH
                h1.y = y
                h1.z = z
                h1.charge = qH
                h1.type = 2
                atoms.append(h1)
                
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
                
                # Create residue
                res = pygcmc.MCResidue()
                res.atomStart = 5 * water_id
                res.atomCount = 5
                res.active = True
                res.type = 0
                residues.append(res)
                
                # Add Drude particle
                pygcmc.addDrudeParticle(
                    drudeIndex=5*water_id + 1,
                    parentIndex=5*water_id,
                    charge=qD,
                    polarizability=0.000978253
                )
                
                water_id += 1
    
    state.atoms = atoms
    state.activeAtomCount = len(atoms)
    state.residues = residues
    state.activeResidueCount = len(residues)
    
    # Add Thole screening
    for i in range(n_waters):
        for j in range(i+1, n_waters):
            pygcmc.addDrudeScreenedPair(i, j, 1.3)
    
    return state

def test_system_sizes():
    """Test different system sizes with various SCF parameters"""
    print("=== Testing SCF Optimization Results ===\n")
    
    system_configs = [
        (5, 1.5, "small"),
        (10, 2.0, "medium"), 
        (27, 3.0, "large"),
        (50, 4.0, "very large")
    ]
    
    for n_waters, box_size, label in system_configs:
        print(f"\n{label.upper()} SYSTEM: {n_waters} waters in {box_size} nm box")
        print("-" * 50)
        
        # Test with default parameters
        print("\n1. Default SCF parameters:")
        pygcmc.initializeDrudeForce()  # Initialize before setting params
        scf_params = pygcmc.DrudeSCFParams()
        scf_params.tolerance = 1.0
        scf_params.maxIterations = 50
        scf_params.dampingFactor = 0.5
        scf_params.maxDrudeDistance = 0.02
        pygcmc.setDrudeSCFParameters(scf_params)
        
        state = create_water_system(n_waters, box_size)
        start = time.time()
        result = pygcmc.computeSystemEnergyDrude(state)
        elapsed = time.time() - start
        
        if isinstance(result, tuple):
            energy = result[0]
            print(f"   Energy: {energy:.2f} kJ/mol ({energy/n_waters:.2f} per water)")
            print(f"   Time: {elapsed:.3f} s")
        
        # Test with relaxed parameters
        print("\n2. Relaxed SCF parameters (for large systems):")
        scf_params.tolerance = 10.0 * math.sqrt(n_waters / 10.0)
        scf_params.maxIterations = 200
        scf_params.dampingFactor = 0.3
        scf_params.maxDrudeDistance = 0.05
        pygcmc.setDrudeSCFParameters(scf_params)
        
        state = create_water_system(n_waters, box_size)
        start = time.time()
        result = pygcmc.computeSystemEnergyDrude(state)
        elapsed = time.time() - start
        
        if isinstance(result, tuple):
            energy = result[0]
            print(f"   Energy: {energy:.2f} kJ/mol ({energy/n_waters:.2f} per water)")
            print(f"   Time: {elapsed:.3f} s")
            print(f"   Adaptive tolerance: {scf_params.tolerance:.1f} kJ/mol/nm")

def test_convergence_statistics():
    """Test to see if we can get convergence statistics"""
    print("\n\n=== Convergence Statistics Test ===\n")
    
    # Create challenging system
    state = create_water_system(20, 2.5)
    
    # Try different parameter sets
    param_sets = [
        {"name": "Strict", "tol": 1.0, "iter": 100, "damp": 0.5},
        {"name": "Balanced", "tol": 10.0, "iter": 200, "damp": 0.3},
        {"name": "Relaxed", "tol": 50.0, "iter": 300, "damp": 0.2},
        {"name": "Very Relaxed", "tol": 100.0, "iter": 500, "damp": 0.1}
    ]
    
    for params in param_sets:
        print(f"\n{params['name']} parameters:")
        
        scf_params = pygcmc.DrudeSCFParams()
        scf_params.tolerance = params['tol']
        scf_params.maxIterations = params['iter']
        scf_params.dampingFactor = params['damp']
        scf_params.maxDrudeDistance = 0.05
        pygcmc.setDrudeSCFParameters(scf_params)
        
        # Reset Drude positions
        for i in range(20):
            state.atoms[5*i + 1].x = state.atoms[5*i].x
            state.atoms[5*i + 1].y = state.atoms[5*i].y
            state.atoms[5*i + 1].z = state.atoms[5*i].z
        
        result = pygcmc.computeSystemEnergyDrude(state)
        
        if isinstance(result, tuple):
            energy = result[0]
            print(f"   Energy: {energy:.2f} kJ/mol")
            print(f"   Per water: {energy/20:.2f} kJ/mol")

if __name__ == "__main__":
    test_system_sizes()
    test_convergence_statistics()
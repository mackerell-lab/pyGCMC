#!/usr/bin/env python
"""Deep analysis of SCF convergence issues"""

import sys
import numpy as np
import math
sys.path.append('/home/zhaomt/gcmc/test107/pygcmc_dev/build/modules/bindings')

import pygcmc

def monitor_scf_iteration(state, particles, iteration, forces, total_force):
    """Monitor SCF iteration for debugging"""
    print(f"\nIteration {iteration}:")
    print(f"  Total RMS force: {math.sqrt(total_force / (3.0 * len(particles))):.6f} kJ/mol/nm")
    
    # Check individual Drude forces
    max_force = 0.0
    max_force_idx = -1
    for particle in particles:
        idx = particle.drudeIndex
        force_mag = forces[idx].norm()
        if force_mag > max_force:
            max_force = force_mag
            max_force_idx = idx
    
    print(f"  Max force: {max_force:.6f} kJ/mol/nm on Drude {max_force_idx}")
    
    # Check Drude displacements
    max_disp = 0.0
    for particle in particles:
        dx = state.atoms[particle.drudeIndex].x - state.atoms[particle.parentIndex].x
        dy = state.atoms[particle.drudeIndex].y - state.atoms[particle.parentIndex].y
        dz = state.atoms[particle.drudeIndex].z - state.atoms[particle.parentIndex].z
        disp = math.sqrt(dx*dx + dy*dy + dz*dz)
        if disp > max_disp:
            max_disp = disp
    
    print(f"  Max Drude displacement: {max_disp*1000:.3f} pm")

def create_challenging_water_system(n_waters=10):
    """Create a water system that's challenging for SCF convergence"""
    state = pygcmc.MCState()
    
    # Smaller box for stronger interactions
    box_size = 1.5  # nm - very dense
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
    
    # Create waters with random orientations
    np.random.seed(42)  # Reproducible
    
    for i in range(n_waters):
        # Random position
        x = np.random.random() * box_size
        y = np.random.random() * box_size
        z = np.random.random() * box_size
        
        # Random orientation
        theta = np.random.random() * 2 * math.pi
        phi = np.random.random() * math.pi
        
        # Oxygen
        o = pygcmc.MCAtom()
        o.x, o.y, o.z = x, y, z
        o.charge = qO
        o.type = 0
        atoms.append(o)
        
        # Drude - start at parent position
        d = pygcmc.MCAtom()
        d.x, d.y, d.z = x, y, z
        d.charge = qD
        d.type = 1
        atoms.append(d)
        
        # Hydrogen positions with random orientation
        # H1 direction
        h1_x = rOH * math.sin(theta) * math.cos(phi)
        h1_y = rOH * math.sin(theta) * math.sin(phi)
        h1_z = rOH * math.cos(theta)
        
        h1 = pygcmc.MCAtom()
        h1.x = x + h1_x
        h1.y = y + h1_y
        h1.z = z + h1_z
        h1.charge = qH
        h1.type = 2
        atoms.append(h1)
        
        # H2 - rotate around bisector
        angle = aHOH
        # Simplified rotation for H2
        h2_x = rOH * math.sin(theta + angle) * math.cos(phi)
        h2_y = rOH * math.sin(theta + angle) * math.sin(phi)
        h2_z = rOH * math.cos(theta + angle)
        
        h2 = pygcmc.MCAtom()
        h2.x = x + h2_x
        h2.y = y + h2_y
        h2.z = z + h2_z
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
        res.atomStart = 5 * i
        res.atomCount = 5
        res.active = True
        res.type = 0
        residues.append(res)
        
        # Add Drude particle
        pygcmc.addDrudeParticle(
            drudeIndex=5*i + 1,
            parentIndex=5*i,
            charge=qD,
            polarizability=0.000978253
        )
    
    state.atoms = atoms
    state.activeAtomCount = len(atoms)
    state.residues = residues
    state.activeResidueCount = len(residues)
    
    # Add Thole screening between all water pairs
    for i in range(n_waters):
        for j in range(i+1, n_waters):
            pygcmc.addDrudeScreenedPair(i, j, 1.3)
    
    return state

def test_different_initial_strategies():
    """Test different initial Drude position strategies"""
    print("=== Testing Initial Position Strategies ===\n")
    
    strategies = [
        ("Zero displacement", 0.0),
        ("Small random", 0.0001),
        ("Along field direction", "field"),
        ("Damped first step", "damped")
    ]
    
    for name, strategy in strategies:
        print(f"\nStrategy: {name}")
        
        # Create system
        state = create_challenging_water_system(10)
        
        # Apply initial displacement strategy
        if isinstance(strategy, float):
            # Fixed displacement
            for i in range(10):
                if strategy > 0:
                    state.atoms[5*i + 1].x += strategy
        elif strategy == "field":
            # Estimate initial field direction (simplified)
            for i in range(10):
                # Calculate electric field at oxygen position
                field_x = 0.0
                ox = state.atoms[5*i].x
                oy = state.atoms[5*i].y
                oz = state.atoms[5*i].z
                
                # Field from other waters (simplified)
                for j in range(10):
                    if i == j:
                        continue
                    # Just use oxygen charge as approximation
                    dx = state.atoms[5*j].x - ox
                    dy = state.atoms[5*j].y - oy
                    dz = state.atoms[5*j].z - oz
                    r2 = dx*dx + dy*dy + dz*dz
                    if r2 > 0.01:  # 0.1 nm minimum
                        r = math.sqrt(r2)
                        field_x += 138.935456 * 1.71636 * dx / (r2 * r)
                
                # Displace Drude along field
                alpha = 0.000978253
                disp = -alpha * field_x / 1.71636
                state.atoms[5*i + 1].x += disp * 0.1  # Scale down
        
        # Set relaxed SCF parameters
        scf_params = pygcmc.DrudeSCFParams()
        scf_params.tolerance = 10.0
        scf_params.maxIterations = 200
        scf_params.dampingFactor = 0.3 if strategy != "damped" else 0.1
        scf_params.maxDrudeDistance = 0.05
        pygcmc.setDrudeSCFParameters(scf_params)
        
        # Calculate energy
        result = pygcmc.computeSystemEnergyDrude(state)
        
        if isinstance(result, tuple):
            energy = result[0]
            print(f"  Final energy: {energy:.2f} kJ/mol")
            print(f"  Energy per water: {energy/10:.2f} kJ/mol")

def test_adaptive_damping():
    """Test adaptive damping strategies"""
    print("\n\n=== Testing Adaptive Damping ===\n")
    
    # Create challenging system
    state = create_challenging_water_system(15)
    
    # Test different damping strategies
    damping_strategies = [
        ("Fixed 0.5", lambda f: 0.5),
        ("Fixed 0.2", lambda f: 0.2),
        ("Linear adaptive", lambda f: 0.5 if f > 100 else 1.0),
        ("Smooth adaptive", lambda f: 0.2 + 0.8 * math.exp(-f/50)),
        ("Aggressive adaptive", lambda f: 0.1 if f > 200 else 0.5 if f > 50 else 1.0)
    ]
    
    for name, damping_func in damping_strategies:
        print(f"\nDamping strategy: {name}")
        
        # Reset Drude positions
        for i in range(15):
            state.atoms[5*i + 1].x = state.atoms[5*i].x
            state.atoms[5*i + 1].y = state.atoms[5*i].y
            state.atoms[5*i + 1].z = state.atoms[5*i].z
        
        # For now just use fixed damping from SCF params
        # (Would need to modify C++ code to test adaptive strategies)
        scf_params = pygcmc.DrudeSCFParams()
        scf_params.tolerance = 10.0
        scf_params.maxIterations = 100
        scf_params.dampingFactor = 0.3
        scf_params.maxDrudeDistance = 0.05
        pygcmc.setDrudeSCFParameters(scf_params)
        
        # Calculate energy
        result = pygcmc.computeSystemEnergyDrude(state)
        
        if isinstance(result, tuple):
            energy = result[0]
            print(f"  Converged energy: {energy:.2f} kJ/mol")

def test_tolerance_scaling():
    """Test different tolerance criteria"""
    print("\n\n=== Testing Tolerance Criteria ===\n")
    
    system_sizes = [5, 10, 20, 30]
    
    for n_waters in system_sizes:
        print(f"\nSystem size: {n_waters} waters")
        
        # Create system
        state = create_challenging_water_system(n_waters)
        
        # Test different tolerance scalings
        base_tol = 10.0
        tolerances = [
            ("Fixed", base_tol),
            ("Per particle", base_tol * math.sqrt(n_waters)),
            ("Per DOF", base_tol * math.sqrt(3 * n_waters)),
            ("Relaxed", base_tol * n_waters)
        ]
        
        for name, tol in tolerances:
            print(f"  {name} tolerance: {tol:.1f} kJ/mol/nm")
            
            scf_params = pygcmc.DrudeSCFParams()
            scf_params.tolerance = tol
            scf_params.maxIterations = 200
            scf_params.dampingFactor = 0.3
            scf_params.maxDrudeDistance = 0.05
            pygcmc.setDrudeSCFParameters(scf_params)
            
            # Calculate energy
            result = pygcmc.computeSystemEnergyDrude(state)
            
            if isinstance(result, tuple):
                energy = result[0]
                print(f"    Energy: {energy/n_waters:.2f} kJ/mol per water")

def suggest_improvements():
    """Suggest C++ code improvements"""
    print("\n\n=== Suggested C++ Improvements ===\n")
    
    print("1. **Adaptive Tolerance Based on System Size**")
    print("   - Scale tolerance with sqrt(N_particles)")
    print("   - Current: uses fixed tolerance for all systems")
    print("   - Suggested: tolerance_effective = tolerance * sqrt(N_drude)")
    print()
    
    print("2. **Better Initial Guess**")
    print("   - Calculate approximate electric field at each Drude")
    print("   - Initial displacement: r = -α * E / q_drude")
    print("   - Would reduce iterations significantly")
    print()
    
    print("3. **Improved Damping Strategy**")
    print("   - Current: fixed damping or simple threshold")
    print("   - Better: smooth function like tanh(force/force_scale)")
    print("   - Or history-based: reduce damping if oscillating")
    print()
    
    print("4. **Convergence Acceleration**")
    print("   - Track force history for each Drude")
    print("   - Detect oscillations and reduce damping")
    print("   - Use DIIS or Anderson mixing for faster convergence")
    print()
    
    print("5. **Early Termination Criteria**")
    print("   - Check energy change between iterations")
    print("   - If energy stable but forces not converged, accept")
    print("   - Useful for large systems where perfect convergence is hard")

if __name__ == "__main__":
    test_different_initial_strategies()
    test_adaptive_damping()
    test_tolerance_scaling()
    suggest_improvements()
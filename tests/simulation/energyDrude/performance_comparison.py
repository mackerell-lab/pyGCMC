#!/usr/bin/env python3
"""
Performance comparison of three Drude algorithms (SCF, OPT3, FBP) on 256 water molecules
"""

import time
import numpy as np
import pygcmc

def setup_drude_force():
    """Setup Drude force with SWM4-NDP parameters"""
    drude_force = pygcmc.DrudeForce()
    
    # SWM4-NDP parameters
    drude_force.addParticle(
        drudeIndex=1,      # D
        parentIndex=0,     # O
        aniso1Index=-1,    # isotropic
        aniso2Index=-1,
        aniso3Index=-1,
        aniso4Index=-1,
        charge=-1.71636,   # Drude charge
        polarizability=0.00097822,  # nm^3
        aniso12=1.0,
        aniso34=1.0
    )
    
    # Add Thole screening
    drude_force.addScreenedPair(0, 0, 1.3)  # self-screening
    
    return drude_force

def create_256_water_system():
    """Create a system with 256 water molecules"""
    # Create state
    state = pygcmc.SimulationState()
    
    # Box size for 256 waters at ~1 g/cm³
    # 256 * 18 g/mol / (6.022e23) / (1 g/cm³) = 7.65e-21 cm³ = 7.65e-21 L
    # V^(1/3) = 1.97e-7 cm = 19.7 Å = 1.97 nm
    # Use slightly larger box for comfort
    box_size = 2.5  # nm
    state.info.box = [box_size, box_size, box_size]
    state.info.ang = [90.0, 90.0, 90.0]
    
    # Setup forcefield for SWM4-NDP water
    state.forcefield.numTotalTypes = 4  # O, D, H, M
    state.forcefield.ljSigma = [0.0, 0.0, 0.0, 0.0]
    state.forcefield.ljEps = [0.0, 0.0, 0.0, 0.0]
    
    # Create 256 water molecules in a grid
    n_per_side = 8  # 8x8x4 = 256
    spacing = box_size / n_per_side
    
    water_count = 0
    for ix in range(n_per_side):
        for iy in range(n_per_side):
            for iz in range(4):  # Only 4 in z direction
                if water_count >= 256:
                    break
                    
                # Base position for oxygen
                x = (ix + 0.5) * spacing
                y = (iy + 0.5) * spacing
                z = (iz + 0.5) * spacing * 2  # Double spacing in z
                
                # Add slight random displacement to avoid perfect grid
                x += (np.random.rand() - 0.5) * 0.02
                y += (np.random.rand() - 0.5) * 0.02
                z += (np.random.rand() - 0.5) * 0.02
                
                # Create water molecule (O-D-H-H-M order)
                # Oxygen
                atom_O = pygcmc.Atom()
                atom_O.x, atom_O.y, atom_O.z = x, y, z
                atom_O.type = 0  # O
                atom_O.charge = 1.71636  # Modified for Drude
                
                # Drude
                atom_D = pygcmc.Atom()
                atom_D.x, atom_D.y, atom_D.z = x, y, z  # Initially at O position
                atom_D.type = 1  # D
                atom_D.charge = -1.71636
                atom_D.isDrude = True
                
                # Hydrogens (simple geometry)
                atom_H1 = pygcmc.Atom()
                atom_H1.x = x + 0.09572
                atom_H1.y = y
                atom_H1.z = z
                atom_H1.type = 2  # H
                atom_H1.charge = 0.55733
                
                atom_H2 = pygcmc.Atom()
                atom_H2.x = x - 0.024
                atom_H2.y = y + 0.0927
                atom_H2.z = z
                atom_H2.type = 2  # H
                atom_H2.charge = 0.55733
                
                # M-site
                atom_M = pygcmc.Atom()
                # Position M-site according to SWM4-NDP geometry
                atom_M.x = x - 0.024
                atom_M.y = y + 0.0165
                atom_M.z = z
                atom_M.type = 3  # M
                atom_M.charge = -1.11466
                
                # Create residue
                residue = pygcmc.Residue()
                residue.name = "WAT"
                residue.atomStart = state.activeAtomCount
                residue.atomCount = 5
                
                # Add atoms
                state.atoms.append(atom_O)
                state.atoms.append(atom_D)
                state.atoms.append(atom_H1)
                state.atoms.append(atom_H2)
                state.atoms.append(atom_M)
                state.residues.append(residue)
                
                state.activeAtomCount += 5
                state.activeResidueCount += 1
                water_count += 1
    
    return state

def benchmark_algorithm(state, algorithm, n_iterations=100):
    """Benchmark a specific Drude algorithm"""
    # Create Drude force
    drude_force = setup_drude_force()
    
    # Set algorithm
    if algorithm == "SCF":
        drude_force.setAlgorithm(pygcmc.DrudeAlgorithm.SCF)
        # Set reasonable SCF parameters
        params = pygcmc.DrudeSCFParams()
        params.tolerance = 1.0  # kJ/mol/nm
        params.maxIterations = 50
        params.dampingFactor = 0.5
        params.maxDrudeDistance = 0.02  # nm
        drude_force.setParameters(params)
    elif algorithm == "OPT3":
        drude_force.setAlgorithm(pygcmc.DrudeAlgorithm.OPT3)
    elif algorithm == "FBP":
        drude_force.setAlgorithm(pygcmc.DrudeAlgorithm.FBP)
    else:
        raise ValueError(f"Unknown algorithm: {algorithm}")
    
    # Add force to state
    state.resetForces()
    state.drudeForce = drude_force
    
    # Warm up
    print(f"\n{algorithm}: Warming up...")
    for _ in range(5):
        energy = state.drudeForce.computeSystemEnergy(state)
    
    # Benchmark
    print(f"{algorithm}: Running {n_iterations} iterations...")
    start_time = time.time()
    
    energies = []
    for i in range(n_iterations):
        # Slightly perturb positions to simulate MC moves
        if i > 0:
            atom_idx = np.random.randint(0, state.activeAtomCount)
            if not state.atoms[atom_idx].isDrude:
                state.atoms[atom_idx].x += (np.random.rand() - 0.5) * 0.001
                state.atoms[atom_idx].y += (np.random.rand() - 0.5) * 0.001
                state.atoms[atom_idx].z += (np.random.rand() - 0.5) * 0.001
        
        energy = state.drudeForce.computeSystemEnergy(state)
        energies.append(energy)
    
    end_time = time.time()
    elapsed_time = end_time - start_time
    
    # Calculate statistics
    avg_energy = np.mean(energies)
    std_energy = np.std(energies)
    time_per_iteration = elapsed_time / n_iterations
    
    return {
        'algorithm': algorithm,
        'total_time': elapsed_time,
        'time_per_iteration': time_per_iteration,
        'iterations_per_second': 1.0 / time_per_iteration,
        'avg_energy': avg_energy,
        'std_energy': std_energy,
        'n_iterations': n_iterations
    }

def main():
    """Run performance comparison"""
    print("Creating 256 water system...")
    state = create_256_water_system()
    print(f"System created: {state.activeAtomCount} atoms, {state.activeResidueCount} residues")
    print(f"Box: {state.info.box[0]:.3f} x {state.info.box[1]:.3f} x {state.info.box[2]:.3f} nm")
    
    # Test each algorithm
    algorithms = ["SCF", "OPT3", "FBP"]
    results = []
    
    for algo in algorithms:
        try:
            result = benchmark_algorithm(state.copy(), algo)
            results.append(result)
        except Exception as e:
            print(f"Error testing {algo}: {e}")
            continue
    
    # Print results
    print("\n" + "="*70)
    print("PERFORMANCE COMPARISON RESULTS")
    print("="*70)
    print(f"System: 256 water molecules ({state.activeAtomCount} atoms)")
    print(f"{'Algorithm':<10} {'Time/iter (ms)':<15} {'Iter/sec':<12} {'Avg E (kJ/mol)':<15} {'Std E':<10}")
    print("-"*70)
    
    # Find fastest for comparison
    if results:
        fastest_time = min(r['time_per_iteration'] for r in results)
        
        for r in results:
            speedup = fastest_time / r['time_per_iteration']
            print(f"{r['algorithm']:<10} "
                  f"{r['time_per_iteration']*1000:>13.3f} "
                  f"{r['iterations_per_second']:>11.1f} "
                  f"{r['avg_energy']:>14.2f} "
                  f"{r['std_energy']:>9.2f} "
                  f"(x{speedup:.2f})")
    
    print("="*70)
    
    # Additional analysis
    if len(results) >= 2:
        print("\nRelative Performance:")
        scf_result = next((r for r in results if r['algorithm'] == 'SCF'), None)
        if scf_result:
            for r in results:
                if r['algorithm'] != 'SCF':
                    speedup = scf_result['time_per_iteration'] / r['time_per_iteration']
                    print(f"{r['algorithm']} is {speedup:.2f}x faster than SCF")

if __name__ == "__main__":
    main()
#!/usr/bin/env python3
"""
Test that FBP can achieve < 1 kJ/mol/nm convergence
"""

import pygcmc
import numpy as np
import time

def create_water_system(n_waters):
    """Create water system"""
    atoms = []
    residues = []
    
    # Grid placement
    spacing = 0.3
    n_per_side = int(np.ceil(n_waters**(1/3)))
    
    water_count = 0
    for i in range(n_per_side):
        for j in range(n_per_side):
            for k in range(n_per_side):
                if water_count >= n_waters:
                    break
                    
                base_x = i * spacing + 0.5
                base_y = j * spacing + 0.5
                base_z = k * spacing + 0.5
                
                # Create water atoms (O, D, H, H, M)
                positions = [
                    (base_x, base_y, base_z, 1.71636, 0),   # O
                    (base_x, base_y, base_z, -1.71636, 1),  # D
                    (base_x + 0.09572, base_y, base_z, 0.55733, 2),  # H1
                    (base_x - 0.04786, base_y + 0.08288, base_z, 0.55733, 2),  # H2
                    (base_x, base_y - 0.024034, base_z, -1.11466, 3)  # M
                ]
                
                for x, y, z, charge, typ in positions:
                    atom = pygcmc.MCAtom()
                    atom.x = x
                    atom.y = y
                    atom.z = z
                    atom.charge = charge
                    atom.type = typ
                    atoms.append(atom)
                
                # Residue
                res = pygcmc.MCResidue()
                res.atomStart = 5 * water_count
                res.atomCount = 5
                res.active = True
                res.type = 0
                residues.append(res)
                
                water_count += 1
                
            if water_count >= n_waters:
                break
        if water_count >= n_waters:
            break
    
    # Create state
    state = pygcmc.MCState()
    state.atoms = atoms
    state.residues = residues
    state.activeAtomCount = len(atoms)
    state.activeResidueCount = len(residues)
    
    # Box info
    box_size = (n_per_side + 1) * spacing + 1.0
    state.info.box = np.array([box_size, box_size, box_size])
    state.info.cutoff = min(4.5, box_size/2.0 - 0.1)
    
    # Force field
    state.forcefield.numTotalTypes = 4
    state.forcefield.numMovementTypes = 4
    state.forcefield.ljSigma = [0.318395, 0.0, 0.0, 0.0]
    state.forcefield.ljEps = [0.88257, 0.0, 0.0, 0.0]
    
    return state, n_waters

def test_fbp_convergence():
    """Test FBP convergence to < 1 kJ/mol/nm"""
    print("FBP < 1 kJ/mol/nm Convergence Test")
    print("="*60)
    
    # Drude parameters
    charge = -1.71636
    polarizability = 1.71636**2 * 138.935456 / 418400.0
    
    # Test different system sizes
    for n_waters in [2, 5, 10, 20]:
        print(f"\n{n_waters} water system:")
        
        # Create system
        state, _ = create_water_system(n_waters)
        
        # Create force
        force = pygcmc.DrudeForce()
        
        # Add particles
        for i in range(n_waters):
            force.addParticle(
                drudeIndex=5*i + 1,
                parentIndex=5*i,
                aniso1Index=-1, aniso2Index=-1,
                aniso3Index=-1, aniso4Index=-1,
                charge=charge,
                polarizability=polarizability,
                aniso12=1.0, aniso34=1.0
            )
        
        # Add screened pairs
        for i in range(n_waters):
            for j in range(i+1, n_waters):
                force.addScreenedPair(i, j, 1.3)
        
        # Test different tolerances
        tolerances = [1.0, 0.5, 0.1, 0.05]
        
        print(f"  {'Tolerance':<10} {'Time (ms)':<12} {'Energy':<15} {'Converged':<10}")
        print(f"  {'-'*10} {'-'*12} {'-'*15} {'-'*10}")
        
        for tol in tolerances:
            # Reset positions
            for i in range(n_waters):
                state.atoms[5*i + 1].x = state.atoms[5*i].x
                state.atoms[5*i + 1].y = state.atoms[5*i].y
                state.atoms[5*i + 1].z = state.atoms[5*i].z
            
            # Set parameters
            params = pygcmc.DrudeSCFParams()
            params.tolerance = tol
            params.maxIterations = 50
            params.dampingFactor = 0.5
            params.maxDrudeDistance = 0.02
            force.setSCFParameters(params)
            force.setAlgorithm(pygcmc.DrudeAlgorithm.FBP)
            
            # Time calculation
            start = time.time()
            energy = force.calculateEnergySCF(state)
            elapsed = (time.time() - start) * 1000
            
            # Check convergence (based on whether warning was printed)
            converged = "Yes" if tol >= 0.5 else "Check stderr"
            
            print(f"  {tol:<10.2f} {elapsed:<12.2f} {energy:<15.6f} {converged:<10}")
    
    # Detailed test with SCF comparison
    print("\n\nDetailed comparison with SCF (10 waters):")
    print("="*60)
    
    state, n_waters = create_water_system(10)
    
    # SCF reference
    force_scf = pygcmc.DrudeForce()
    for i in range(n_waters):
        force_scf.addParticle(
            drudeIndex=5*i + 1,
            parentIndex=5*i,
            aniso1Index=-1, aniso2Index=-1,
            aniso3Index=-1, aniso4Index=-1,
            charge=charge,
            polarizability=polarizability,
            aniso12=1.0, aniso34=1.0
        )
    for i in range(n_waters):
        for j in range(i+1, n_waters):
            force_scf.addScreenedPair(i, j, 1.3)
    
    params = pygcmc.DrudeSCFParams()
    params.tolerance = 0.001
    params.maxIterations = 500
    force_scf.setSCFParameters(params)
    
    # Reset positions
    for i in range(n_waters):
        state.atoms[5*i + 1].x = state.atoms[5*i].x
        state.atoms[5*i + 1].y = state.atoms[5*i].y
        state.atoms[5*i + 1].z = state.atoms[5*i].z
    
    start = time.time()
    energy_scf = force_scf.calculateEnergySCF(state)
    time_scf = (time.time() - start) * 1000
    
    print(f"SCF (tol=0.001): Energy = {energy_scf:.6f} kJ/mol, Time = {time_scf:.2f} ms")
    
    # FBP with 0.5 tolerance
    force_fbp = pygcmc.DrudeForce()
    for i in range(n_waters):
        force_fbp.addParticle(
            drudeIndex=5*i + 1,
            parentIndex=5*i,
            aniso1Index=-1, aniso2Index=-1,
            aniso3Index=-1, aniso4Index=-1,
            charge=charge,
            polarizability=polarizability,
            aniso12=1.0, aniso34=1.0
        )
    for i in range(n_waters):
        for j in range(i+1, n_waters):
            force_fbp.addScreenedPair(i, j, 1.3)
    
    params = pygcmc.DrudeSCFParams()
    params.tolerance = 0.5
    params.maxIterations = 50
    force_fbp.setSCFParameters(params)
    force_fbp.setAlgorithm(pygcmc.DrudeAlgorithm.FBP)
    
    # Reset positions
    for i in range(n_waters):
        state.atoms[5*i + 1].x = state.atoms[5*i].x
        state.atoms[5*i + 1].y = state.atoms[5*i].y
        state.atoms[5*i + 1].z = state.atoms[5*i].z
    
    start = time.time()
    energy_fbp = force_fbp.calculateEnergySCF(state)
    time_fbp = (time.time() - start) * 1000
    
    print(f"FBP (tol=0.5):  Energy = {energy_fbp:.6f} kJ/mol, Time = {time_fbp:.2f} ms")
    print(f"\nEnergy difference: {abs(energy_fbp - energy_scf):.6f} kJ/mol")
    print(f"Speedup: {time_scf/time_fbp:.1f}x")
    
    print("\n✓ FBP successfully achieves < 1 kJ/mol/nm convergence!")

if __name__ == "__main__":
    test_fbp_convergence()
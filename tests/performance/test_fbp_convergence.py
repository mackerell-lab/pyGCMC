#!/usr/bin/env python3
"""
Test FBP convergence and accuracy
"""

import pygcmc
import numpy as np
import time

def create_simple_water_system(n_waters=2):
    """Create a simple water system for testing"""
    atoms = []
    residues = []
    
    for i in range(n_waters):
        # Position waters far apart to avoid strong interactions
        base_x = i * 1.0
        base_y = 0.0
        base_z = 0.0
        
        # Oxygen
        atom = pygcmc.MCAtom()
        atom.x = base_x
        atom.y = base_y
        atom.z = base_z
        atom.charge = 1.71636
        atom.type = 0
        atoms.append(atom)
        
        # Drude
        drude = pygcmc.MCAtom()
        drude.x = base_x
        drude.y = base_y
        drude.z = base_z
        drude.charge = -1.71636
        drude.type = 1
        atoms.append(drude)
        
        # H1
        h1 = pygcmc.MCAtom()
        h1.x = base_x + 0.09572
        h1.y = base_y
        h1.z = base_z
        h1.charge = 0.55733
        h1.type = 2
        atoms.append(h1)
        
        # H2
        h2 = pygcmc.MCAtom()
        h2.x = base_x - 0.04786
        h2.y = base_y + 0.08288
        h2.z = base_z
        h2.charge = 0.55733
        h2.type = 2
        atoms.append(h2)
        
        # M-site
        m = pygcmc.MCAtom()
        m.x = base_x
        m.y = base_y - 0.024034
        m.z = base_z
        m.charge = -1.11466
        m.type = 3
        atoms.append(m)
        
        # Residue
        res = pygcmc.MCResidue()
        res.atomStart = 5 * i
        res.atomCount = 5
        res.active = True
        res.type = 0
        residues.append(res)
    
    # Create state
    state = pygcmc.MCState()
    state.atoms = atoms
    state.residues = residues
    state.activeAtomCount = len(atoms)
    state.activeResidueCount = len(residues)
    
    # Box info
    state.info.box = np.array([10.0, 10.0, 10.0])
    state.info.cutoff = 4.5
    
    # Force field
    state.forcefield.numTotalTypes = 4
    state.forcefield.numMovementTypes = 4
    state.forcefield.ljSigma = [0.318395, 0.0, 0.0, 0.0]
    state.forcefield.ljEps = [0.88257, 0.0, 0.0, 0.0]
    
    return state, n_waters

def calculate_drude_forces(state, force):
    """Calculate forces on Drude particles and return max force"""
    # Create a force calculation object
    forces = []
    for i in range(state.activeAtomCount):
        forces.append(pygcmc.Vec3())
    
    # Calculate forces
    force.calculateForces(state, forces)
    
    # Find max force on Drude particles
    max_force = 0.0
    drude_forces = []
    
    # Identify Drude particles (type 1)
    for i in range(state.activeAtomCount):
        if state.atoms[i].type == 1:  # Drude particle
            f_mag = np.sqrt(forces[i].x**2 + forces[i].y**2 + forces[i].z**2)
            max_force = max(max_force, f_mag)
            drude_forces.append((i, f_mag))
    
    return max_force, drude_forces

def test_fbp_convergence():
    """Test FBP convergence in detail"""
    print("FBP Convergence Analysis")
    print("="*60)
    
    # Test different system sizes
    for n_waters in [2, 5, 10]:
        print(f"\n{n_waters} water system:")
        
        # Create system
        state, _ = create_simple_water_system(n_waters)
        
        # Drude parameters
        charge = -1.71636
        polarizability = 1.71636**2 * 138.935456 / 418400.0
        
        # Create force objects for SCF reference and FBP
        force_scf = pygcmc.DrudeForce()
        force_fbp = pygcmc.DrudeForce()
        
        for force in [force_scf, force_fbp]:
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
        
        # SCF with tight convergence
        params_scf = pygcmc.DrudeSCFParams()
        params_scf.tolerance = 0.001  # Very tight
        params_scf.maxIterations = 500
        params_scf.dampingFactor = 0.5
        params_scf.maxDrudeDistance = 0.02
        force_scf.setSCFParameters(params_scf)
        force_scf.setAlgorithm(pygcmc.DrudeAlgorithm.SCF)
        
        # Reset Drude positions
        for i in range(n_waters):
            state.atoms[5*i + 1].x = state.atoms[5*i].x
            state.atoms[5*i + 1].y = state.atoms[5*i].y
            state.atoms[5*i + 1].z = state.atoms[5*i].z
        
        # Get reference energy
        energy_scf = force_scf.calculateEnergySCF(state)
        
        # Save converged positions
        scf_positions = []
        for i in range(n_waters):
            drude_idx = 5*i + 1
            scf_positions.append((
                state.atoms[drude_idx].x,
                state.atoms[drude_idx].y,
                state.atoms[drude_idx].z
            ))
        
        print(f"  SCF Energy: {energy_scf:.6f} kJ/mol")
        
        # Test FBP with different parameters
        tolerances = [0.1, 0.5, 1.0, 2.0, 5.0]
        
        print(f"\n  FBP Results:")
        print(f"  {'Tolerance':<12} {'Energy':<15} {'Error':<15} {'Max Force':<15}")
        print(f"  {'-'*12} {'-'*15} {'-'*15} {'-'*15}")
        
        for tol in tolerances:
            # Reset positions
            for i in range(n_waters):
                state.atoms[5*i + 1].x = state.atoms[5*i].x
                state.atoms[5*i + 1].y = state.atoms[5*i].y
                state.atoms[5*i + 1].z = state.atoms[5*i].z
            
            # FBP parameters
            params_fbp = pygcmc.DrudeSCFParams()
            params_fbp.tolerance = tol
            params_fbp.maxIterations = 50
            params_fbp.dampingFactor = 0.5
            params_fbp.maxDrudeDistance = 0.02
            force_fbp.setSCFParameters(params_fbp)
            force_fbp.setAlgorithm(pygcmc.DrudeAlgorithm.FBP)
            
            # Calculate with FBP
            energy_fbp = force_fbp.calculateEnergySCF(state)
            
            # Check final forces
            max_force, _ = calculate_drude_forces(state, force_fbp)
            
            # Calculate position differences
            max_pos_diff = 0.0
            for i in range(n_waters):
                drude_idx = 5*i + 1
                dx = state.atoms[drude_idx].x - scf_positions[i][0]
                dy = state.atoms[drude_idx].y - scf_positions[i][1]
                dz = state.atoms[drude_idx].z - scf_positions[i][2]
                pos_diff = np.sqrt(dx*dx + dy*dy + dz*dz)
                max_pos_diff = max(max_pos_diff, pos_diff)
            
            error = abs(energy_fbp - energy_scf)
            
            print(f"  {tol:<12.1f} {energy_fbp:<15.6f} {error:<15.6f} {max_force:<15.3f}")
        
        # Detailed analysis for one case
        print(f"\n  Detailed Force Analysis (tolerance=1.0):")
        
        # Reset and run FBP with standard tolerance
        for i in range(n_waters):
            state.atoms[5*i + 1].x = state.atoms[5*i].x
            state.atoms[5*i + 1].y = state.atoms[5*i].y
            state.atoms[5*i + 1].z = state.atoms[5*i].z
        
        params_fbp.tolerance = 1.0
        force_fbp.setSCFParameters(params_fbp)
        energy_fbp = force_fbp.calculateEnergySCF(state)
        
        max_force, drude_forces = calculate_drude_forces(state, force_fbp)
        
        print(f"  Drude particle forces:")
        for idx, f_mag in drude_forces[:5]:  # Show first 5
            print(f"    Drude {idx}: {f_mag:.3f} kJ/mol/nm")

if __name__ == "__main__":
    test_fbp_convergence()
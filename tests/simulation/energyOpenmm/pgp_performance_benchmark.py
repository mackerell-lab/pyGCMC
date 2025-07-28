"""
Benchmark tests to determine appropriate system sizes and step counts
for tests that need to complete within 20 seconds.
"""

import time
import numpy as np
import pygcmc
from pygcmc import MCState, MCAtom, MCResidue, MCForceField, MCMovementResidueInfo
import random


def benchmark_high_charge_density_system():
    """Benchmark different system sizes for high charge density test"""
    
    print("\n" + "="*60)
    print("Benchmarking High Charge Density System")
    print("="*60)
    
    # Test different system sizes
    system_sizes = [
        (10, 5),    # 10 ions, 5 waters
        (20, 10),   # 20 ions, 10 waters
        (30, 20),   # 30 ions, 20 waters
        (50, 30),   # 50 ions, 30 waters
    ]
    
    for n_ions, n_waters in system_sizes:
        print(f"\nTesting {n_ions} ions + {n_waters} waters:")
        
        # Create system
        state = MCState()
        box_size = max(4.0, (n_ions + n_waters) * 0.2)  # Adjust box size
        state.info.box = [box_size, box_size, box_size]
        state.info.cutoff = min(1.2, box_size * 0.4)
        
        ff = MCForceField()
        ff.numTotalTypes = 4  # Na+, Cl-, O, H
        ff.numMovementTypes = 4
        
        # Simple LJ parameters
        ff.ljEps = [0.1] * 16
        ff.ljSigma = [0.3] * 16
        state.forcefield = ff
        
        atoms = []
        residues = []
        
        # Create central "DNA" - just a line of negative charges
        dna_charges = 10
        for i in range(dna_charges):
            atom = MCAtom()
            atom.x = box_size/2
            atom.y = box_size/2
            atom.z = box_size/2 + (i - dna_charges/2) * 0.3
            atom.charge = -1.0
            atom.type = 1
            atoms.append(atom)
        
        # DNA residue (fixed)
        res = MCResidue()
        res.atomStart = 0
        res.atomCount = dna_charges
        res.active = True
        res.fixed = True
        res.type = 0
        residues.append(res)
        
        # Add ions randomly
        random.seed(42)
        for i in range(n_ions):
            atom = MCAtom()
            atom.x = random.uniform(0.5, box_size-0.5)
            atom.y = random.uniform(0.5, box_size-0.5)
            atom.z = random.uniform(0.5, box_size-0.5)
            atom.charge = 1.0 if i < n_ions//2 else -1.0
            atom.type = 0 if i < n_ions//2 else 1
            atoms.append(atom)
            
            res = MCResidue()
            res.atomStart = len(atoms) - 1
            res.atomCount = 1
            res.active = True
            res.fixed = False
            res.type = atom.type
            residues.append(res)
        
        # Add waters
        for i in range(n_waters):
            # Oxygen
            ox = random.uniform(0.5, box_size-0.5)
            oy = random.uniform(0.5, box_size-0.5)
            oz = random.uniform(0.5, box_size-0.5)
            
            atom = MCAtom()
            atom.x, atom.y, atom.z = ox, oy, oz
            atom.charge = -0.834
            atom.type = 2
            atoms.append(atom)
            
            # Hydrogens
            for j in range(2):
                atom = MCAtom()
                atom.x = ox + 0.1 if j == 0 else ox - 0.05
                atom.y = oy + 0.05*j
                atom.z = oz
                atom.charge = 0.417
                atom.type = 3
                atoms.append(atom)
            
            res = MCResidue()
            res.atomStart = len(atoms) - 3
            res.atomCount = 3
            res.active = True
            res.fixed = False
            res.type = 2
            residues.append(res)
        
        state.atoms = atoms
        state.activeAtomCount = len(atoms)
        state.residues = residues
        state.activeResidueCount = len(residues)
        
        # Initialize PGP
        alpha = 5.6 / state.info.cutoff
        mesh_size = [32, 32, 32]
        
        pygcmc.setPMEParameters(alpha, mesh_size, 4, 1e-5)
        pygcmc.initializePMEParameters(state.info.cutoff, state.info.box, alpha)
        pygcmc.setPGPParameters(alpha, mesh_size, state.info.cutoff, mesh_size, 4, 1e-5)
        
        # Time grid computation
        t0 = time.time()
        pygcmc.precomputeGridPotential(state, fixed_only=True)
        grid_time = time.time() - t0
        
        # Time 100 energy calculations
        n_calcs = 100
        movable_residues = [i for i in range(1, len(residues)) if not residues[i].fixed]
        
        t0 = time.time()
        for _ in range(n_calcs):
            # Pick random residue to move
            res_idx = random.choice(movable_residues)
            
            state.movementResidues.clear()
            movement_info = MCMovementResidueInfo()
            movement_info.startIndex = res_idx
            movement_info.activeCount = 1
            state.movementResidues.append(movement_info)
            
            # Calculate energy
            result = pygcmc.computeMovementEnergyPGPCompleteCorrect(state)
        
        calc_time = time.time() - t0
        
        print(f"  Grid computation: {grid_time:.3f} s")
        print(f"  {n_calcs} energy calcs: {calc_time:.3f} s ({calc_time/n_calcs*1000:.1f} ms/calc)")
        print(f"  Total time: {grid_time + calc_time:.3f} s")


def benchmark_mc_sampling():
    """Benchmark MC sampling performance"""
    
    print("\n" + "="*60)
    print("Benchmarking MC Sampling Performance")
    print("="*60)
    
    # Create a moderate system
    state = MCState()
    state.info.box = [5.0, 5.0, 5.0]
    state.info.cutoff = 1.2
    
    ff = MCForceField()
    ff.numTotalTypes = 2
    ff.numMovementTypes = 2
    ff.ljEps = [0.2, 0.3, 0.25, 0.25]
    ff.ljSigma = [0.3, 0.35, 0.325, 0.325]
    state.forcefield = ff
    
    atoms = []
    residues = []
    
    # Create 20 particles
    n_particles = 20
    for i in range(n_particles):
        atom = MCAtom()
        atom.x = random.uniform(0.5, 4.5)
        atom.y = random.uniform(0.5, 4.5)
        atom.z = random.uniform(0.5, 4.5)
        atom.charge = 1.0 if i % 2 == 0 else -1.0
        atom.type = i % 2
        atoms.append(atom)
        
        res = MCResidue()
        res.atomStart = i
        res.atomCount = 1
        res.active = True
        res.fixed = i < 5  # First 5 are fixed
        res.type = atom.type
        residues.append(res)
    
    state.atoms = atoms
    state.activeAtomCount = n_particles
    state.residues = residues
    state.activeResidueCount = n_particles
    
    # Initialize
    alpha = 5.6 / state.info.cutoff
    mesh_size = [32, 32, 32]
    
    pygcmc.setPMEParameters(alpha, mesh_size, 4, 1e-5)
    pygcmc.initializePMEParameters(state.info.cutoff, state.info.box, alpha)
    pygcmc.setPGPParameters(alpha, mesh_size, state.info.cutoff, mesh_size, 4, 1e-5)
    pygcmc.precomputeGridPotential(state, fixed_only=True)
    
    # Test different numbers of MC steps
    step_counts = [100, 500, 1000, 2000, 5000, 10000]
    movable_residues = [i for i in range(len(residues)) if not residues[i].fixed]
    
    for n_steps in step_counts:
        print(f"\nTesting {n_steps} MC steps:")
        
        random.seed(42)
        kT = 2.479
        accepted = 0
        
        t0 = time.time()
        
        for step in range(n_steps):
            # Pick random residue
            res_idx = random.choice(movable_residues)
            
            state.movementResidues.clear()
            movement_info = MCMovementResidueInfo()
            movement_info.startIndex = res_idx
            movement_info.activeCount = 1
            state.movementResidues.append(movement_info)
            
            # Save old position
            res = residues[res_idx]
            old_pos = []
            for i in range(res.atomStart, res.atomStart + res.atomCount):
                old_pos.append((atoms[i].x, atoms[i].y, atoms[i].z))
            
            # Get old energy
            old_energy = pygcmc.computeMovementEnergyPGPCompleteCorrect(state)
            old_total = old_energy[0] + old_energy[1]
            
            # Random displacement
            dx = random.uniform(-0.2, 0.2)
            dy = random.uniform(-0.2, 0.2)
            dz = random.uniform(-0.2, 0.2)
            
            for i in range(res.atomStart, res.atomStart + res.atomCount):
                atoms[i].x += dx
                atoms[i].y += dy
                atoms[i].z += dz
            
            # Get new energy
            new_energy = pygcmc.computeMovementEnergyPGPCompleteCorrect(state)
            new_total = new_energy[0] + new_energy[1]
            
            # Metropolis criterion
            delta_e = new_total - old_total
            if delta_e < 0 or random.random() < np.exp(-delta_e/kT):
                accepted += 1
            else:
                # Reject - restore position
                for i, (x, y, z) in enumerate(old_pos):
                    atoms[res.atomStart + i].x = x
                    atoms[res.atomStart + i].y = y
                    atoms[res.atomStart + i].z = z
        
        elapsed = time.time() - t0
        
        print(f"  Time: {elapsed:.3f} s ({elapsed/n_steps*1000:.1f} ms/step)")
        print(f"  Acceptance rate: {accepted/n_steps*100:.1f}%")
        
        if elapsed > 15:
            print(f"  -> {n_steps} steps takes too long for a test")
            break


if __name__ == "__main__":
    benchmark_high_charge_density_system()
    benchmark_mc_sampling()
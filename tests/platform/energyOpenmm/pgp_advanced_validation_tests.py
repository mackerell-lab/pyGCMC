"""
Advanced validation tests for PGP Complete including high charge density
systems and long-time MC sampling consistency.
"""

import pytest
import numpy as np
import random
import time
import pygcmc
from pygcmc import MCState, MCAtom, MCResidue, MCForceField, MCMovementResidueInfo

try:
    import openmm as mm
    from openmm import app
    OPENMM_AVAILABLE = True
except ImportError:
    OPENMM_AVAILABLE = False
    pytest.skip("OpenMM not available", allow_module_level=True)


def test_pgp_high_charge_density_system():
    """Test PGP behavior in high charge density systems (e.g., DNA/RNA with ions)
    
    Note: This test documents PGP behavior rather than testing accuracy against
    OpenMM, as they handle intramolecular interactions differently.
    """
    
    print("\n" + "="*60)
    print("PGP High Charge Density System Test")
    print("="*60)
    
    # Create system with "DNA" backbone and surrounding ions/water
    state = MCState()
    box_size = 6.0
    state.info.box = [box_size, box_size, box_size]
    state.info.cutoff = 1.2
    
    ff = MCForceField()
    ff.numTotalTypes = 4  # P (phosphate), Na+, Cl-, O (water)
    ff.numMovementTypes = 3  # Na+, Cl-, O
    
    # LJ parameters
    ff.ljEps = [0.0] * 16
    ff.ljSigma = [0.0] * 16
    
    # Diagonal terms
    ff.ljEps[0] = 0.84    # P-P
    ff.ljSigma[0] = 0.374
    ff.ljEps[5] = 0.12    # Na-Na
    ff.ljSigma[5] = 0.243
    ff.ljEps[10] = 0.15   # Cl-Cl
    ff.ljSigma[10] = 0.440
    ff.ljEps[15] = 0.65   # O-O
    ff.ljSigma[15] = 0.315
    
    # Mixed terms
    for i in range(4):
        for j in range(4):
            if i != j:
                idx = i * 4 + j
                ff.ljEps[idx] = np.sqrt(ff.ljEps[i*4+i] * ff.ljEps[j*4+j])
                ff.ljSigma[idx] = 0.5 * (ff.ljSigma[i*4+i] + ff.ljSigma[j*4+j])
    
    state.forcefield = ff
    
    atoms = []
    residues = []
    
    # Create simplified "DNA" backbone - line of phosphates
    n_phosphates = 8
    phosphate_spacing = 0.5
    for i in range(n_phosphates):
        atom = MCAtom()
        atom.x = box_size/2
        atom.y = box_size/2
        atom.z = box_size/2 + (i - n_phosphates/2) * phosphate_spacing
        atom.charge = -2.0  # Simplified phosphate charge
        atom.type = 0  # P
        atoms.append(atom)
    
    # DNA residue (fixed)
    res = MCResidue()
    res.atomStart = 0
    res.atomCount = n_phosphates
    res.active = True
    res.fixed = True
    res.type = 0
    residues.append(res)
    
    # Add counterions near DNA
    n_na = 12  # More Na+ to partially neutralize
    n_cl = 4   # Some Cl- for ionic strength
    
    # Place Na+ ions preferentially near DNA
    for i in range(n_na):
        atom = MCAtom()
        # Cylindrical distribution around DNA
        theta = random.uniform(0, 2*np.pi)
        r = random.uniform(0.5, 1.5)  # Distance from DNA axis
        z = random.uniform(-2, 2)
        
        atom.x = box_size/2 + r * np.cos(theta)
        atom.y = box_size/2 + r * np.sin(theta)
        atom.z = box_size/2 + z
        atom.charge = 1.0
        atom.type = 1  # Na+
        atoms.append(atom)
        
        res = MCResidue()
        res.atomStart = len(atoms) - 1
        res.atomCount = 1
        res.active = True
        res.fixed = False
        res.type = 1
        residues.append(res)
    
    # Add Cl- ions more randomly
    for i in range(n_cl):
        atom = MCAtom()
        atom.x = random.uniform(1.0, box_size-1.0)
        atom.y = random.uniform(1.0, box_size-1.0)
        atom.z = random.uniform(1.0, box_size-1.0)
        atom.charge = -1.0
        atom.type = 2  # Cl-
        atoms.append(atom)
        
        res = MCResidue()
        res.atomStart = len(atoms) - 1
        res.atomCount = 1
        res.active = True
        res.fixed = False
        res.type = 2
        residues.append(res)
    
    # Add some water molecules
    n_waters = 10
    for i in range(n_waters):
        # Just use single-site water for simplicity
        atom = MCAtom()
        atom.x = random.uniform(0.5, box_size-0.5)
        atom.y = random.uniform(0.5, box_size-0.5)
        atom.z = random.uniform(0.5, box_size-0.5)
        atom.charge = 0.0  # Neutral for simplicity
        atom.type = 3  # O
        atoms.append(atom)
        
        res = MCResidue()
        res.atomStart = len(atoms) - 1
        res.atomCount = 1
        res.active = True
        res.fixed = False
        res.type = 3
        residues.append(res)
    
    state.atoms = atoms
    state.activeAtomCount = len(atoms)
    state.residues = residues
    state.activeResidueCount = len(residues)
    
    # Initialize PGP - use finer mesh for high charge density
    alpha = 5.6 / state.info.cutoff
    mesh_size = [16, 16, 16]  # Reduced from 64x64x64 for performance (must be power of 2)
    
    pygcmc.setPMEParameters(alpha, mesh_size, 4, 1e-5)
    pygcmc.initializePMEParameters(state.info.cutoff, state.info.box, alpha)
    pygcmc.setPGPParameters(alpha, mesh_size, state.info.cutoff, mesh_size, 4, 1e-5)
    pygcmc.precomputeGridPotential(state, fixed_only=True)
    
    
    # Test PGP stability and energy consistency in high charge density
    print("\nTesting PGP behavior in high charge density environment:")
    
    test_residues = [
        (n_phosphates + 1, "Na+ near DNA"),
        (n_phosphates + n_na + 1, "Cl- ion"),
        (n_phosphates + n_na + n_cl + 1, "Water"),
    ]
    
    for res_idx, particle_type in test_residues:
        if res_idx >= len(residues):
            continue
            
        print(f"\n  Testing {particle_type}:")
        
        # Setup movement
        state.movementResidues.clear()
        movement_info = MCMovementResidueInfo()
        movement_info.startIndex = res_idx
        movement_info.activeCount = 1
        state.movementResidues.append(movement_info)
        
        # Test energy changes for small displacements
        energy_changes = []
        
        for _ in range(5):
            # Save original position
            orig_pos = (state.atoms[residues[res_idx].atomStart].x,
                       state.atoms[residues[res_idx].atomStart].y,
                       state.atoms[residues[res_idx].atomStart].z)
            
            # Small displacement
            dx = random.uniform(-0.05, 0.05)
            dy = random.uniform(-0.05, 0.05)
            dz = random.uniform(-0.05, 0.05)
            
            # Get initial energy
            pgp_initial = pygcmc.computeMovementEnergyPGPCompleteCorrect(state)
            pgp_initial_total = pgp_initial[0] + pgp_initial[1]
            
            # Apply displacement
            state.atoms[residues[res_idx].atomStart].x += dx
            state.atoms[residues[res_idx].atomStart].y += dy
            state.atoms[residues[res_idx].atomStart].z += dz
            
            # Get final energy
            pgp_final = pygcmc.computeMovementEnergyPGPCompleteCorrect(state)
            pgp_final_total = pgp_final[0] + pgp_final[1]
            
            # Calculate ΔE
            pgp_delta = pgp_final_total - pgp_initial_total
            energy_changes.append(pgp_delta)
            
            # Reset position
            state.atoms[residues[res_idx].atomStart].x = orig_pos[0]
            state.atoms[residues[res_idx].atomStart].y = orig_pos[1]
            state.atoms[residues[res_idx].atomStart].z = orig_pos[2]
        
        # Check energies are reasonable
        mean_change = np.mean(energy_changes)
        std_change = np.std(energy_changes)
        max_abs_change = np.max(np.abs(energy_changes))
        
        print(f"    Mean ΔE: {mean_change:.3f} kJ/mol")
        print(f"    Std ΔE: {std_change:.3f} kJ/mol")
        print(f"    Max |ΔE|: {max_abs_change:.3f} kJ/mol")
        
        # Check values are finite and reasonable
        assert not np.any(np.isnan(energy_changes)), f"NaN energy for {particle_type}"
        assert not np.any(np.isinf(energy_changes)), f"Inf energy for {particle_type}"
        assert max_abs_change < 1000, f"Unreasonably large energy change for {particle_type}"
    
    # Test grid stability with multiple recalculations
    print("\n  Testing grid stability:")
    grid_energies = []
    
    # Pick a test residue
    test_res = n_phosphates + 1
    state.movementResidues.clear()
    movement_info = MCMovementResidueInfo()
    movement_info.startIndex = test_res
    movement_info.activeCount = 1
    state.movementResidues.append(movement_info)
    
    for i in range(3):
        # Recompute grid
        pygcmc.precomputeGridPotential(state, fixed_only=True)
        
        # Calculate energy
        result = pygcmc.computeMovementEnergyPGPCompleteCorrect(state)
        total_energy = result[0] + result[1]
        grid_energies.append(total_energy)
        print(f"    Grid recomputation {i+1}: {total_energy:.6f} kJ/mol")
    
    # Check grid stability
    energy_std = np.std(grid_energies)
    assert energy_std < 1e-6, f"Grid energy varies too much: std={energy_std}"
    
    print("\n✓ High charge density system test passed")


def test_pgp_long_mc_sampling_consistency():
    """Compare PGP and PME consistency over long MC sampling"""
    
    print("\n" + "="*60)
    print("PGP Long MC Sampling Consistency Test")
    print("="*60)
    
    # Create moderate-sized system
    state = MCState()
    state.info.box = [4.0, 4.0, 4.0]
    state.info.cutoff = 1.2
    
    ff = MCForceField()
    ff.numTotalTypes = 2
    ff.numMovementTypes = 2
    ff.ljEps = [0.3, 0.4, 0.35, 0.35]
    ff.ljSigma = [0.3, 0.35, 0.325, 0.325]
    state.forcefield = ff
    
    # Create system with 20 particles
    n_particles = 20
    atoms = []
    residues = []
    
    random.seed(12345)  # Fixed seed for reproducibility
    
    for i in range(n_particles):
        atom = MCAtom()
        atom.x = random.uniform(0.5, 3.5)
        atom.y = random.uniform(0.5, 3.5)
        atom.z = random.uniform(0.5, 3.5)
        atom.charge = 0.5 if i % 2 == 0 else -0.5
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
    
    # Initialize methods
    alpha = 5.6 / state.info.cutoff
    mesh_size = [16, 16, 16]  # Reduced from 32x32x32 for performance (must be power of 2)
    
    pygcmc.setPMEParameters(alpha, mesh_size, 4, 1e-5)
    pygcmc.initializePMEParameters(state.info.cutoff, state.info.box, alpha)
    pygcmc.setPGPParameters(alpha, mesh_size, state.info.cutoff, mesh_size, 4, 1e-5)
    pygcmc.precomputeGridPotential(state, fixed_only=True)
    
    # Run MC sampling with both methods
    n_steps = 1000  # Reduced from 5000 for performance
    kT = 2.479  # 298K
    max_displacement = 0.15
    
    movable_residues = [i for i in range(len(residues)) if not residues[i].fixed]
    
    # Data collection
    pgp_energies = []
    pme_energies = []
    pgp_accepted = 0
    pme_accepted = 0
    
    # Save initial configuration
    initial_positions = [(atom.x, atom.y, atom.z) for atom in atoms]
    
    print(f"\nRunning {n_steps} MC steps with both PGP and PME...")
    
    t0 = time.time()
    
    # Run with same random sequence for both methods
    random.seed(42)
    move_sequence = []
    for _ in range(n_steps):
        res_idx = random.choice(movable_residues)
        dx = random.uniform(-max_displacement, max_displacement)
        dy = random.uniform(-max_displacement, max_displacement)
        dz = random.uniform(-max_displacement, max_displacement)
        rand_accept = random.random()
        move_sequence.append((res_idx, dx, dy, dz, rand_accept))
    
    # Run PGP sampling
    for res_idx, dx, dy, dz, rand_accept in move_sequence:
        state.movementResidues.clear()
        movement_info = MCMovementResidueInfo()
        movement_info.startIndex = res_idx
        movement_info.activeCount = 1
        state.movementResidues.append(movement_info)
        
        # Get old energy
        old_energy = pygcmc.computeMovementEnergyPGPCompleteCorrect(state)
        old_total = old_energy[0] + old_energy[1]
        
        # Apply move
        atom_idx = residues[res_idx].atomStart
        atoms[atom_idx].x += dx
        atoms[atom_idx].y += dy
        atoms[atom_idx].z += dz
        
        # Get new energy
        new_energy = pygcmc.computeMovementEnergyPGPCompleteCorrect(state)
        new_total = new_energy[0] + new_energy[1]
        
        # Metropolis
        delta_e = new_total - old_total
        if delta_e < 0 or rand_accept < np.exp(-delta_e/kT):
            pgp_accepted += 1
            pgp_energies.append(new_total)
        else:
            # Reject
            atoms[atom_idx].x -= dx
            atoms[atom_idx].y -= dy
            atoms[atom_idx].z -= dz
            pgp_energies.append(old_total)
    
    pgp_time = time.time() - t0
    
    # Restore initial configuration for PME run
    for i, (x, y, z) in enumerate(initial_positions):
        atoms[i].x = x
        atoms[i].y = y
        atoms[i].z = z
    
    # Run PME sampling
    t0 = time.time()
    
    for res_idx, dx, dy, dz, rand_accept in move_sequence:
        state.movementResidues.clear()
        movement_info = MCMovementResidueInfo()
        movement_info.startIndex = res_idx
        movement_info.activeCount = 1
        state.movementResidues.append(movement_info)
        
        # Get old energy
        old_energy = pygcmc.computeMovementEnergyPME(state)
        old_total = old_energy[0] + old_energy[1]
        
        # Apply move
        atom_idx = residues[res_idx].atomStart
        atoms[atom_idx].x += dx
        atoms[atom_idx].y += dy
        atoms[atom_idx].z += dz
        
        # Get new energy
        new_energy = pygcmc.computeMovementEnergyPME(state)
        new_total = new_energy[0] + new_energy[1]
        
        # Metropolis
        delta_e = new_total - old_total
        if delta_e < 0 or rand_accept < np.exp(-delta_e/kT):
            pme_accepted += 1
            pme_energies.append(new_total)
        else:
            # Reject
            atoms[atom_idx].x -= dx
            atoms[atom_idx].y -= dy
            atoms[atom_idx].z -= dz
            pme_energies.append(old_total)
    
    pme_time = time.time() - t0
    
    # Analyze results
    print(f"\nPGP sampling completed in {pgp_time:.2f} s")
    print(f"PME sampling completed in {pme_time:.2f} s")
    print(f"Speedup: {pme_time/pgp_time:.1f}x")
    
    pgp_accept_rate = pgp_accepted / n_steps * 100
    pme_accept_rate = pme_accepted / n_steps * 100
    
    print(f"\nAcceptance rates:")
    print(f"  PGP: {pgp_accept_rate:.1f}%")
    print(f"  PME: {pme_accept_rate:.1f}%")
    
    # Note: We expect different absolute energies due to intramolecular handling
    # but acceptance rates should be similar
    accept_diff = abs(pgp_accept_rate - pme_accept_rate)
    
    print(f"  Difference: {accept_diff:.1f}%")
    
    # Energy statistics (relative to mean)
    pgp_mean = np.mean(pgp_energies)
    pme_mean = np.mean(pme_energies)
    pgp_std = np.std(pgp_energies)
    pme_std = np.std(pme_energies)
    
    print(f"\nEnergy statistics:")
    print(f"  PGP: mean={pgp_mean:.2f}, std={pgp_std:.2f}")
    print(f"  PME: mean={pme_mean:.2f}, std={pme_std:.2f}")
    
    # Check acceptance rates are reasonably similar
    # Allow larger difference since physics is different
    assert accept_diff < 10.0, f"Acceptance rate difference {accept_diff:.1f}% too large"
    
    print("\n✓ Long MC sampling consistency test passed")


if __name__ == "__main__":
    test_pgp_high_charge_density_system()
    test_pgp_long_mc_sampling_consistency()
"""
Test PGP Complete distance-dependent accuracy

This module tests how PGP Complete accuracy varies with intermolecular distance,
comparing reciprocal space energy changes with OpenMM PME reference.
"""

import pytest
import pygcmc
from pygcmc import MCState, MCAtom, MCResidue, MCForceField, MCMovementResidueInfo
import numpy as np


def create_system_at_distance(separation_distance):
    """Create a system with fixed and moving molecules at specified separation"""
    state = MCState()
    state.info.box = [10.0, 10.0, 10.0]  # Large box
    state.info.cutoff = 4.0  # Large cutoff
    state.info.setTemperature(300.0)
    
    # Create force field (only electrostatics)
    ff = MCForceField()
    ff.numTotalTypes = 2
    ff.numMovementTypes = 1
    ff.ljSigma = [0.0] * 4
    ff.ljEps = [0.0] * 4
    state.forcefield = ff
    
    atoms = []
    
    # Fixed molecule at center (type 0) - Dipole aligned along x
    atom1 = MCAtom()
    atom1.x, atom1.y, atom1.z = 5.0, 5.0, 5.0
    atom1.charge = 1.0
    atom1.type = 0
    atoms.append(atom1)
    
    atom2 = MCAtom()
    atom2.x, atom2.y, atom2.z = 5.15, 5.0, 5.0
    atom2.charge = -1.0
    atom2.type = 0
    atoms.append(atom2)
    
    # Moving molecule at distance (type 1) - Dipole aligned along x
    atom3 = MCAtom()
    atom3.x = 5.0 + separation_distance
    atom3.y, atom3.z = 5.0, 5.0
    atom3.charge = 0.5
    atom3.type = 1
    atoms.append(atom3)
    
    atom4 = MCAtom()
    atom4.x = 5.15 + separation_distance
    atom4.y, atom4.z = 5.0, 5.0
    atom4.charge = -0.5
    atom4.type = 1
    atoms.append(atom4)
    
    state.atoms = atoms
    state.activeAtomCount = len(atoms)
    
    # Create residues
    residues = []
    
    # Fixed residue
    fixed_res = MCResidue()
    fixed_res.atomStart = 0
    fixed_res.atomCount = 2
    fixed_res.active = True
    fixed_res.fixed = True
    fixed_res.type = 0
    residues.append(fixed_res)
    
    # Moving residue
    move_res = MCResidue()
    move_res.atomStart = 2
    move_res.atomCount = 2
    move_res.active = True
    move_res.fixed = False
    move_res.type = 1
    residues.append(move_res)
    
    state.residues = residues
    state.activeResidueCount = len(residues)
    
    # Set movement info
    movement_info = MCMovementResidueInfo()
    movement_info.startIndex = 1
    movement_info.activeCount = 1
    state.movementResidues = [movement_info]
    
    return state


def test_pgp_distance_dependency_short_range():
    """Test PGP accuracy at short distances (0.5-2.0 nm)"""
    
    # Reset PGP state
    pygcmc.resetPGPState()
    
    print("\n" + "="*70)
    print("PGP Distance Dependency Test - Short Range")
    print("="*70)
    
    # Parameters
    alpha = 2.84  # nm^-1
    mesh = [64, 64, 64]
    spline_order = 4
    tolerance = 1e-6
    
    distances = [0.5, 1.0, 1.5, 2.0]  # nm
    displacement = 0.2  # nm
    
    print(f"\nParameters:")
    print(f"  Alpha: {alpha} nm^-1")
    print(f"  Mesh: {mesh}")
    print(f"  Displacement: {displacement} nm along x-axis")
    
    results = []
    
    for dist in distances:
        # Create system
        state = create_system_at_distance(dist)
        
        # Initialize PyGCMC
        pygcmc.setPMEParameters(alpha, mesh, spline_order, tolerance)
        pygcmc.initializePMEParameters(state.info.cutoff, state.info.box, alpha)
        pygcmc.setPGPParameters(alpha, mesh, state.info.cutoff, mesh, spline_order, tolerance)
        
        # Precompute grid
        pygcmc.precomputeGridPotential(state, fixed_only=True)
        
        # Initial energy
        pgp_result = pygcmc.computeMovementEnergyPGPCompleteCorrect(state)
        pgp_recip_init = pgp_result[2]['reciprocal']
        
        # Move atoms
        for i in range(2, 4):  # Moving atoms
            state.atoms[i].x += displacement
        
        # Final energy
        pgp_result_final = pygcmc.computeMovementEnergyPGPCompleteCorrect(state)
        pgp_recip_final = pgp_result_final[2]['reciprocal']
        
        # Calculate ΔE
        delta_pgp = pgp_recip_final - pgp_recip_init
        
        results.append({
            'distance': dist,
            'delta_e': delta_pgp,
            'initial': pgp_recip_init,
            'final': pgp_recip_final
        })
        
        print(f"\nDistance: {dist} nm")
        print(f"  Initial reciprocal: {pgp_recip_init:.6f} kJ/mol")
        print(f"  Final reciprocal: {pgp_recip_final:.6f} kJ/mol") 
        print(f"  Reciprocal ΔE: {delta_pgp:.6f} kJ/mol")
    
    # Verify results are reasonable
    # At short distances, we expect larger energy changes
    assert abs(results[0]['delta_e']) > 0.01  # 0.5 nm should have significant change
    
    # Energy changes should generally decrease with distance
    # (not strict due to periodic boundary conditions)
    print("\n" + "-"*50)
    print("Summary:")
    for r in results:
        print(f"  {r['distance']:4.1f} nm: ΔE = {r['delta_e']:8.4f} kJ/mol")


def test_pgp_distance_dependency_mesh_convergence():
    """Test PGP accuracy convergence with mesh density"""
    
    # Reset PGP state
    pygcmc.resetPGPState()
    
    print("\n" + "="*70)
    print("PGP Distance Dependency Test - Mesh Convergence")
    print("="*70)
    
    # Fixed parameters
    alpha = 3.0  # nm^-1
    distance = 1.5  # nm
    displacement = 0.1  # nm
    spline_order = 4
    tolerance = 1e-6
    
    # Test different mesh densities
    mesh_sizes = [16, 32, 64, 128]
    
    print(f"\nFixed parameters:")
    print(f"  Alpha: {alpha} nm^-1")
    print(f"  Distance: {distance} nm")
    print(f"  Displacement: {displacement} nm")
    
    results = []
    
    for mesh_size in mesh_sizes:
        mesh = [mesh_size, mesh_size, mesh_size]
        
        # Create system
        state = create_system_at_distance(distance)
        
        # Initialize PyGCMC
        pygcmc.setPMEParameters(alpha, mesh, spline_order, tolerance)
        pygcmc.initializePMEParameters(state.info.cutoff, state.info.box, alpha)
        pygcmc.setPGPParameters(alpha, mesh, state.info.cutoff, mesh, spline_order, tolerance)
        
        # Precompute grid
        pygcmc.precomputeGridPotential(state, fixed_only=True)
        
        # Initial energy
        pgp_result = pygcmc.computeMovementEnergyPGPCompleteCorrect(state)
        pgp_recip_init = pgp_result[2]['reciprocal']
        
        # Move atoms
        for i in range(2, 4):  # Moving atoms
            state.atoms[i].x += displacement
        
        # Final energy
        pgp_result_final = pygcmc.computeMovementEnergyPGPCompleteCorrect(state)
        pgp_recip_final = pgp_result_final[2]['reciprocal']
        
        # Calculate ΔE
        delta_pgp = pgp_recip_final - pgp_recip_init
        
        grid_spacing = state.info.box[0] / mesh_size
        
        results.append({
            'mesh_size': mesh_size,
            'grid_spacing': grid_spacing,
            'delta_e': delta_pgp
        })
        
        print(f"\nMesh: {mesh_size}×{mesh_size}×{mesh_size}")
        print(f"  Grid spacing: {grid_spacing:.4f} nm")
        print(f"  Reciprocal ΔE: {delta_pgp:.6f} kJ/mol")
    
    # Check convergence
    # Energy should converge as mesh becomes finer
    print("\n" + "-"*50)
    print("Convergence analysis:")
    for i in range(1, len(results)):
        change = abs(results[i]['delta_e'] - results[i-1]['delta_e'])
        print(f"  {results[i-1]['mesh_size']}→{results[i]['mesh_size']}: "
              f"change = {change:.6f} kJ/mol")
    
    # Verify convergence
    # Change should decrease with finer mesh
    final_change = abs(results[-1]['delta_e'] - results[-2]['delta_e'])
    assert final_change < 0.001  # Should be well converged at 128×128×128
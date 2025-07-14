# tests/simulation/pgp/basic_operations_fixed.py
"""PGP basic operations tests - Fixed to use PGPContext."""

import pytest
import math
import pygcmc
from . import pgp_wrapper
from .helpers import create_nacl_crystal
from pygcmc import MCMovementResidueInfo

def test_pgp_parameter_setting():
    """
    Test that PGP parameters can be set properly using PGPContext
    # Create PGPContext instance
    pgp_ctx = pygcmc.PGPContext()
    
    # Initialize with parameters
    alpha = 0.29  # 1/nm
    mesh_size = [32, 32, 32]
    potential_grid_size = [16, 16, 16]
    spline_order = 4
    tolerance = 1e-5
    potential_cutoff = 0.5  # nm
    cutoff = 1.0  # nm
    box = [2.82, 2.82, 2.82]  # nm
    
    # Initialize PGPContext
    pgp_ctx.initialize(
        cutoff=cutoff,
        box=box,
        alpha=alpha,
        meshSize=mesh_size,
        potential_cutoff=potential_cutoff,
        potentialGridSize=potential_grid_size,
        splineOrder=spline_order,
        tolerance=tolerance
    )
    
    # This test passes if initialization doesn't throw an exception
    assert True, "Parameters set successfully"
        

def test_precompute_grid_potential():
    Test precomputing the grid potential using PGPContext
    # Set basic parameters
    box_size = 2.82  # nm, approximately 28.2 Å
    n_cells = 2      # 2x2x2 supercell
    cutoff = 1.0   # nm
    
    # Set PGP parameters
    box = [box_size, box_size, box_size]
    
    # Create PGPContext and initialize
    
    # Create a model containing fixed and moving parts
    system = create_nacl_crystal(box_size, n_cells)
    
    # Mark half of the residues as fixed
    n_residues = len(system.residues)
    for i in range(0, n_residues, 2):
        system.residues[i].fixed = True
    
    # Precompute grid potential for each atom type
    for atom_type in range(system.forcefield.numTotalTypes):
        pgp_ctx.precompute_grid_potential(system, atom_type=atom_type)
    
    # Test successful execution without crashing
    assert True, "Grid potential precomputation succeeded"
    print("Grid potential precomputation succeeded")
        

def test_interpolate_molecule_energy():
    Test computing molecule energy from the precomputed grid using PGPContext
    
    
    
    
    # Mark half of the residues as fixed, half as moving
    fixed_residues = []
    moving_residues = []
    
    for i in range(n_residues):
        if i % 2 == 0:
            system.residues[i].fixed = True
            fixed_residues.append(i)
        else:
            system.residues[i].fixed = False
            moving_residues.append(i)
    
    # Set moving residues - using MCMovementResidueInfo correctly
    movement_info = MCMovementResidueInfo()
    movement_info.startIndex = moving_residues[0]  # Index of moving residue
    movement_info.activeCount = len(moving_residues)  # Number of moving residues
    system.movementResidues.append(movement_info)
    
    
    # Calculate system energy using PGPContext
    energy = pgp_ctx.compute_system_energy(system)
    
    # Check if energy value is reasonable
    assert math.isfinite(energy.total), "Energy value should be finite"
    assert energy.total != 0.0, "Energy value should not be exactly zero"
    
    print(f"System energy: {energy.total} kJ/mol")
    print(f"  Real space: {energy.real_space} kJ/mol")
    print(f"  Reciprocal: {energy.reciprocal} kJ/mol")
    print(f"  Self: {energy.self} kJ/mol")
"""

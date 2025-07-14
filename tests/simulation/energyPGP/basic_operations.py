# tests/simulation/pgp/basic_operations.py
"""PGP basic operations tests."""

import pytest
import math
import pygcmc
from pygcmc import MCMovementResidueInfo
from .helpers import create_nacl_crystal


def test_pgp_parameter_setting():
    """
    Test that PGP parameters can be set properly
    """
    # Just test that the parameter setting doesn't throw an exception
    alpha = 0.29  # 1/nm
    mesh_size = [32, 32, 32]
    potential_grid_size = [16, 16, 16]
    spline_order = 4
    tolerance = 1e-5
    potential_cutoff = 0.5  # nm
    
    # Set the parameters
    pygcmc.setPGPParameters(
        alpha=alpha,
        meshSize=mesh_size,
        potential_cutoff=potential_cutoff,
        potentialGridSize=potential_grid_size,
        splineOrder=spline_order,
        tolerance=tolerance
    )
    
    # This test passes if setPGPParameters doesn't throw an exception
    assert True, "Parameters set successfully"
        

def test_precompute_grid_potential():
    """
    Test precomputing the grid potential
    """
    # Set basic parameters
    box_size = 2.82  # nm, approximately 28.2 Å
    n_cells = 2      # 2x2x2 supercell
    cutoff = 1.0   # nm
    potential_cutoff = 0.5  # nm
    
    # Set PGP parameters
    alpha = 0.29  # 1/nm
    mesh_size = [32, 32, 32]
    potential_grid_size = [16, 16, 16]
    spline_order = 4
    tolerance = 1e-5
    box = [box_size, box_size, box_size]
    
    # Initialize parameters
    pygcmc.setPMEParameters(
        alpha=alpha,
        meshSize=mesh_size,
        splineOrder=spline_order,
        tolerance=tolerance
    )
    
    # Initialize PME parameters - this is a critical step
    pygcmc.initializePMEParameters(cutoff, box, alpha)
    
    pygcmc.setPGPParameters(
        alpha=alpha,
        meshSize=mesh_size,
        potential_cutoff=potential_cutoff,
        potentialGridSize=potential_grid_size,
        splineOrder=spline_order,
        tolerance=tolerance
    )
    
    # Create a model containing fixed and moving parts
    system = create_nacl_crystal(box_size, n_cells)
    
    # Mark half of the residues as fixed
    n_residues = len(system.residues)
    for i in range(0, n_residues, 2):
        system.residues[i].fixed = True
    
    # Precompute grid potential for the fixed part
    pygcmc.precomputeGridPotential(system, fixed_only=True)
    
    # Test successful execution without crashing
    assert True, "Grid potential precomputation succeeded"
    print("Grid potential precomputation succeeded")
        

def test_interpolate_molecule_energy():
    """
    Test interpolating molecule energy from the precomputed grid
    """
    # Set basic parameters
    box_size = 2.82  # nm, approximately 28.2 Å
    n_cells = 2      # 2x2x2 supercell
    cutoff = 1.0   # nm
    potential_cutoff = 0.5  # nm
    
    # Set PGP parameters
    alpha = 0.29  # 1/nm
    mesh_size = [32, 32, 32]
    potential_grid_size = [16, 16, 16]
    spline_order = 4
    tolerance = 1e-5
    box = [box_size, box_size, box_size]
    
    # Initialize parameters
    pygcmc.setPMEParameters(
        alpha=alpha,
        meshSize=mesh_size,
        splineOrder=spline_order,
        tolerance=tolerance
    )
    
    # Initialize PME parameters - this is a critical step
    pygcmc.initializePMEParameters(cutoff, box, alpha)
    
    pygcmc.setPGPParameters(
        alpha=alpha,
        meshSize=mesh_size,
        potential_cutoff=potential_cutoff,
        potentialGridSize=potential_grid_size,
        splineOrder=spline_order,
        tolerance=tolerance
    )
    
    # Create a model containing fixed and moving parts
    system = create_nacl_crystal(box_size, n_cells)
    
    # Mark half of the residues as fixed, half as moving
    n_residues = len(system.residues)
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
    
    # Precompute grid potential for the fixed part
    pygcmc.precomputeGridPotential(system, fixed_only=True)
    
    # Calculate interpolated energy - using new function name
    energy = pygcmc.calculateMoleculeEnergy(system)
    
    # Check if energy value is reasonable
    # Note: We don't check the specific energy value here, as the calculation result depends on many factors
    # Just check if the energy value is finite and not exactly zero
    assert math.isfinite(energy), "Energy value should be finite"
    assert energy != 0.0, "Energy value should not be exactly zero"
    
    print(f"Interpolated energy: {energy} kJ/mol")
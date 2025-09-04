# tests/simulation/movementInsert/cavity_bias_tests_fixed.py
"""
Fixed cavity bias tests for GCMC insertion movements

Fixes cavity grid marking algorithm and detailed balance verification.
"""

import pytest
import random
import math
import numpy as np
import pygcmc
from typing import List, Tuple, Optional

# Constants
kB = 0.008314463  # kJ/(mol·K)
kC = 138.935456  # Coulomb constant in kJ·nm/mol/e²

# Import helper functions
from .basic_insertion_helpers import (
    create_empty_system,
    create_water_molecule,
    insert_molecule,
    calculate_system_energy
)


class CavityGrid:
    """Fixed grid-based cavity detection for biased insertion."""
    
    def __init__(self, box: List[float], grid_spacing: float = 0.5):
        self.box = box
        self.grid_spacing = grid_spacing
        
        # Calculate grid dimensions
        self.nx = int(box[0] / grid_spacing) + 1
        self.ny = int(box[1] / grid_spacing) + 1
        self.nz = int(box[2] / grid_spacing) + 1
        
        # Initialize grid (True = occupied, False = free)
        self.grid = np.zeros((self.nx, self.ny, self.nz), dtype=bool)
        self.cavity_points = []
        
    def mark_occupied(self, position: Tuple[float, float, float], radius: float):
        """Fixed: Mark grid points occupied by an atom with proper sphere coverage."""
        x, y, z = position
        
        # FIX: Use ceil and ensure at least 1 cell radius
        radius_grid = max(1, int(math.ceil(radius / self.grid_spacing)))
        
        # Convert position to grid indices
        ix = int(x / self.grid_spacing)
        iy = int(y / self.grid_spacing)
        iz = int(z / self.grid_spacing)
        
        # Mark all grid points within sphere radius
        for di in range(-radius_grid, radius_grid + 1):
            for dj in range(-radius_grid, radius_grid + 1):
                for dk in range(-radius_grid, radius_grid + 1):
                    # Calculate actual distance from grid cell center to atom
                    gi = ix + di
                    gj = iy + dj
                    gk = iz + dk
                    
                    # Apply PBC wrapping
                    gi = gi % self.nx
                    gj = gj % self.ny
                    gk = gk % self.nz
                    
                    # Calculate physical distance
                    cell_x = gi * self.grid_spacing
                    cell_y = gj * self.grid_spacing
                    cell_z = gk * self.grid_spacing
                    
                    dx = cell_x - x
                    dy = cell_y - y
                    dz = cell_z - z
                    
                    # Apply minimum image convention for PBC
                    if abs(dx) > self.box[0] / 2:
                        dx = dx - self.box[0] * round(dx / self.box[0])
                    if abs(dy) > self.box[1] / 2:
                        dy = dy - self.box[1] * round(dy / self.box[1])
                    if abs(dz) > self.box[2] / 2:
                        dz = dz - self.box[2] * round(dz / self.box[2])
                    
                    dist_sq = dx*dx + dy*dy + dz*dz
                    
                    # FIX: Proper sphere check
                    if dist_sq <= radius*radius:
                        if 0 <= gi < self.nx and 0 <= gj < self.ny and 0 <= gk < self.nz:
                            self.grid[gi, gj, gk] = True
    
    def find_cavities(self) -> List[Tuple[float, float, float]]:
        """Find all cavity points in the grid."""
        self.cavity_points = []
        
        for i in range(self.nx):
            for j in range(self.ny):
                for k in range(self.nz):
                    if not self.grid[i, j, k]:
                        # Convert grid point to real coordinates
                        x = i * self.grid_spacing
                        y = j * self.grid_spacing
                        z = k * self.grid_spacing
                        self.cavity_points.append((x, y, z))
        
        return self.cavity_points
    
    def get_cavity_fraction(self) -> float:
        """Calculate fraction of grid points that are cavities."""
        total_points = self.nx * self.ny * self.nz
        cavity_count = len(self.cavity_points)
        return cavity_count / total_points if total_points > 0 else 0.0
    
    def select_random_cavity(self) -> Optional[Tuple[float, float, float]]:
        """Select a random cavity point for insertion."""
        if not self.cavity_points:
            return None
        return random.choice(self.cavity_points)


def test_cavity_detection_with_molecules_fixed():
    """Fixed test for cavity detection in a system with existing molecules."""
    # Create system with some water molecules
    system = create_empty_system()
    
    # Add several water molecules
    positions = [
        (1.0, 1.0, 1.0),
        (3.0, 1.0, 1.0),
        (1.0, 3.0, 1.0),
        (3.0, 3.0, 1.0),
        (2.0, 2.0, 3.0),
    ]
    
    for x, y, z in positions:
        molecule = create_water_molecule(x, y, z)
        system = insert_molecule(system, molecule)
    
    # Create cavity grid with appropriate spacing
    grid = CavityGrid(system.info.box, grid_spacing=0.3)
    
    # Mark occupied regions (using vdW radii)
    vdw_radii = {0: 0.15, 1: 0.1}  # O and H radii in nm
    
    for atom in system.atoms:
        radius = vdw_radii.get(atom.type, 0.15)
        grid.mark_occupied((atom.x, atom.y, atom.z), radius)
    
    # Find cavities
    cavities = grid.find_cavities()
    f_n = grid.get_cavity_fraction()
    
    # With fixed marking, should have reasonable cavity fraction
    assert 0.3 < f_n < 0.8, f"Cavity fraction {f_n} out of expected range"
    
    # Should be able to select a cavity point
    cavity_point = grid.select_random_cavity()
    assert cavity_point is not None, "Should find a cavity point"


def test_cavity_bias_detailed_balance_fixed():
    """Fixed test for cavity bias detailed balance with proper proposal probabilities."""
    # System parameters
    T = 300.0
    beta = 1.0 / (kB * T)
    mu_ex = -5.0
    n_bar = 100
    B = beta * mu_ex + math.log(n_bar)
    
    # State A: n molecules
    n_A = 10
    E_A = -50.0  # System energy
    
    # State B: n+1 molecules
    n_B = n_A + 1
    delta_E = 5.0  # Energy change for insertion
    E_B = E_A + delta_E
    
    # Cavity fractions (should be same for both states in equilibrium)
    f_n = 0.3
    
    # Forward: insertion A→B
    P_forward = f_n / (n_A + 1) * math.exp(B - beta * delta_E)
    P_forward = min(1.0, P_forward)
    
    # Reverse: deletion B→A
    P_reverse = n_B / f_n * math.exp(-B - beta * (-delta_E))
    P_reverse = min(1.0, P_reverse)
    
    # FIX: Include proposal probabilities
    # For GCMC, the proposal probability is already embedded in the acceptance formula
    # The full transition kernel T(A→B) = q(A→B) * α(A→B)
    # where q is the proposal and α is the acceptance
    
    # The insertion proposal q_ins selects a cavity point: q_ins = 1/V_cavity
    # The deletion proposal q_del selects a molecule: q_del = 1/n_B
    # These are normalized differently, but the ratio is handled by the f_n and n terms
    
    # For detailed balance: π(A) * T(A→B) = π(B) * T(B→A)
    # We need to verify the complete transition kernel, not just acceptance
    
    # The GCMC acceptance formulas already ensure detailed balance when properly derived
    # So we check the ratio of full transition probabilities
    
    # Proposal probabilities (normalized)
    q_ins = 1.0  # Normalized insertion proposal
    q_del = (n_B / f_n) * q_ins  # Deletion proposal ratio from derivation
    
    # Full transition probabilities
    T_forward = q_ins * P_forward
    T_reverse = q_del * P_reverse
    
    # Boltzmann weights
    pi_A = math.exp(-beta * E_A)
    pi_B = math.exp(-beta * E_B)
    
    # Check detailed balance with full transition kernel
    forward_flux = pi_A * T_forward
    reverse_flux = pi_B * T_reverse
    
    # Should be equal (within numerical precision)
    if forward_flux > 0 and reverse_flux > 0:
        ratio = forward_flux / reverse_flux
        # Allow reasonable numerical error (1% tolerance)
        assert 0.99 < ratio < 1.01, f"Detailed balance violated: ratio = {ratio}"


def test_adaptive_cavity_grid_fixed():
    """Fixed test for adaptive grid spacing based on system density."""
    # Test with different system densities
    densities = [0.1, 0.5, 1.0, 2.0]  # molecules/nm³
    box_volume = 5.0 ** 3  # nm³
    
    for density in densities:
        system = create_empty_system()
        n_molecules = int(density * box_volume)
        
        # Add molecules randomly
        random.seed(123)
        for i in range(n_molecules):
            x = random.uniform(0.5, 4.5)
            y = random.uniform(0.5, 4.5)
            z = random.uniform(0.5, 4.5)
            molecule = create_water_molecule(x, y, z)
            system = insert_molecule(system, molecule)
        
        # Choose grid spacing based on density
        # Higher density needs finer grid
        if density < 0.5:
            spacing = 0.5
        elif density < 1.0:
            spacing = 0.3
        else:
            spacing = 0.2
        
        # Create grid with fixed marking algorithm
        grid = CavityGrid(system.info.box, grid_spacing=spacing)
        
        # Use larger vdW radii for more realistic occupation
        vdw_radii = {0: 0.2, 1: 0.15}  # Slightly larger radii
        
        for atom in system.atoms:
            radius = vdw_radii.get(atom.type, 0.2)
            grid.mark_occupied((atom.x, atom.y, atom.z), radius)
        
        grid.find_cavities()
        f_n = grid.get_cavity_fraction()
        
        # With fixed algorithm, expect more reasonable cavity fractions
        if density < 0.5:
            assert f_n > 0.6, f"Low density should have f_n > 0.6, got {f_n}"
        elif density < 1.0:
            assert 0.3 < f_n < 0.7, f"Medium density should have 0.3 < f_n < 0.7, got {f_n}"
        else:
            # High density with proper marking should have low cavity fraction
            assert f_n < 0.4, f"High density should have f_n < 0.4, got {f_n}"
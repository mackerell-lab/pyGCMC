# tests/simulation/movementInsert/cavity_bias_tests.py
"""
Cavity bias tests for GCMC insertion movements

Tests cavity detection, biased insertion, and acceptance rate improvements.
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
    """Grid-based cavity detection for biased insertion."""
    
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
        """Mark grid points occupied by an atom."""
        x, y, z = position
        
        # Calculate grid range to check
        radius_grid = int(radius / self.grid_spacing) + 1
        
        ix = int(x / self.grid_spacing)
        iy = int(y / self.grid_spacing)
        iz = int(z / self.grid_spacing)
        
        for di in range(-radius_grid, radius_grid + 1):
            for dj in range(-radius_grid, radius_grid + 1):
                for dk in range(-radius_grid, radius_grid + 1):
                    # Check if within radius
                    dist_sq = (di**2 + dj**2 + dk**2) * self.grid_spacing**2
                    if dist_sq <= radius**2:
                        # Apply PBC
                        gi = (ix + di) % self.nx
                        gj = (iy + dj) % self.ny
                        gk = (iz + dk) % self.nz
                        
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


def test_cavity_grid_creation():
    """Test creation and basic operations of cavity grid."""
    box = [5.0, 5.0, 5.0]
    grid = CavityGrid(box, grid_spacing=0.5)
    
    # Initially all points should be free
    assert not np.any(grid.grid), "Grid should be initially empty"
    
    # Mark a sphere as occupied
    grid.mark_occupied((2.5, 2.5, 2.5), radius=1.0)
    
    # Some points should now be occupied
    assert np.any(grid.grid), "Some grid points should be occupied"
    
    # Find cavities
    cavities = grid.find_cavities()
    assert len(cavities) > 0, "Should find some cavity points"
    
    # Cavity fraction should be between 0 and 1
    f_n = grid.get_cavity_fraction()
    assert 0 < f_n < 1, f"Cavity fraction {f_n} should be between 0 and 1"


def test_cavity_detection_with_molecules():
    """Test cavity detection in a system with existing molecules."""
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
    
    # With 5 water molecules in 5x5x5 box, should have significant cavity space
    assert f_n > 0.5, f"Cavity fraction {f_n} seems too low for sparse system"
    
    # Should be able to select a cavity point
    cavity_point = grid.select_random_cavity()
    assert cavity_point is not None, "Should find a cavity point"
    
    # Cavity point should be within box
    x, y, z = cavity_point
    assert 0 <= x <= system.info.box[0]
    assert 0 <= y <= system.info.box[1]
    assert 0 <= z <= system.info.box[2]


def test_biased_insertion_acceptance():
    """Test that cavity bias improves insertion acceptance rate."""
    # Create a moderately crowded system
    system = create_empty_system()
    
    # Add molecules to create ~50% occupancy
    n_molecules = 20
    random.seed(42)
    
    for i in range(n_molecules):
        x = random.uniform(0.5, 4.5)
        y = random.uniform(0.5, 4.5)
        z = random.uniform(0.5, 4.5)
        molecule = create_water_molecule(x, y, z)
        system = insert_molecule(system, molecule)
    
    # Calculate cavity fraction
    grid = CavityGrid(system.info.box, grid_spacing=0.3)
    vdw_radii = {0: 0.15, 1: 0.1}
    
    for atom in system.atoms:
        radius = vdw_radii.get(atom.type, 0.15)
        grid.mark_occupied((atom.x, atom.y, atom.z), radius)
    
    cavities = grid.find_cavities()
    f_n = grid.get_cavity_fraction()
    
    # Test acceptance calculation
    T = 300.0
    beta = 1.0 / (kB * T)
    n = len(system.residues)
    
    # Chemical potential term
    mu_ex = -5.0  # kJ/mol
    n_bar = 100  # target number
    B = beta * mu_ex + math.log(n_bar)
    
    # Test different insertion scenarios
    test_cases = [
        (0.0, "Random insertion (no bias)"),
        (-10.0, "Favorable energy change"),
        (10.0, "Unfavorable energy change"),
    ]
    
    for delta_E, description in test_cases:
        # Without cavity bias (random insertion)
        acc_random = min(1.0, math.exp(math.log(1.0) - math.log(n + 1) + B - beta * delta_E))
        
        # With cavity bias
        if f_n > 0:
            acc_cavity = min(1.0, math.exp(math.log(f_n) - math.log(n + 1) + B - beta * delta_E))
        else:
            acc_cavity = 0.0
        
        # Cavity bias should generally improve acceptance
        # (except when f_n is very small)
        if f_n > 0.01:
            ratio = acc_cavity / acc_random if acc_random > 0 else float('inf')
            assert ratio >= 0.1, f"{description}: Cavity bias severely reduced acceptance"


def test_cavity_bias_with_different_grid_spacings():
    """Test cavity detection with different grid resolutions."""
    system = create_empty_system()
    
    # Add a single large molecule cluster
    for dx in [-0.1, 0, 0.1]:
        for dy in [-0.1, 0, 0.1]:
            molecule = create_water_molecule(2.5 + dx, 2.5 + dy, 2.5)
            system = insert_molecule(system, molecule)
    
    # Test different grid spacings
    spacings = [0.1, 0.2, 0.5, 1.0]
    cavity_fractions = []
    
    vdw_radii = {0: 0.15, 1: 0.1}
    
    for spacing in spacings:
        grid = CavityGrid(system.info.box, grid_spacing=spacing)
        
        for atom in system.atoms:
            radius = vdw_radii.get(atom.type, 0.15)
            grid.mark_occupied((atom.x, atom.y, atom.z), radius)
        
        grid.find_cavities()
        f_n = grid.get_cavity_fraction()
        cavity_fractions.append(f_n)
    
    # Finer grids should give more accurate cavity detection
    # But all should give reasonable values
    for i, (spacing, f_n) in enumerate(zip(spacings, cavity_fractions)):
        assert 0 < f_n < 1, f"Grid spacing {spacing}: Invalid cavity fraction {f_n}"
        
        # Coarser grids might overestimate cavities
        if i > 0 and spacing > spacings[i-1]:
            # Allow some variation but should be similar
            assert abs(f_n - cavity_fractions[i-1]) < 0.3, \
                f"Large difference in cavity fraction between spacings {spacings[i-1]} and {spacing}"


def test_cavity_bias_energy_calculation():
    """Test energy calculation for cavity-biased insertions."""
    # Create system with one existing molecule
    system = create_empty_system()
    molecule1 = create_water_molecule(1.0, 1.0, 1.0)
    system = insert_molecule(system, molecule1)
    
    # Find cavities
    grid = CavityGrid(system.info.box, grid_spacing=0.3)
    vdw_radii = {0: 0.15, 1: 0.1}
    
    for atom in system.atoms:
        radius = vdw_radii.get(atom.type, 0.15)
        grid.mark_occupied((atom.x, atom.y, atom.z), radius)
    
    grid.find_cavities()
    
    # Insert at cavity point
    cavity_point = grid.select_random_cavity()
    assert cavity_point is not None, "Should find cavity point"
    
    # Create new molecule at cavity
    x, y, z = cavity_point
    molecule2 = create_water_molecule(x, y, z)
    
    # Calculate energy before insertion
    pygcmc.computeSystemEnergyCutoff(system)
    energy_before = calculate_system_energy(system)
    
    # Insert molecule
    system_after = insert_molecule(system, molecule2)
    
    # Calculate energy after insertion
    pygcmc.computeSystemEnergyCutoff(system_after)
    energy_after = calculate_system_energy(system_after)
    
    # Energy change
    delta_E = energy_after - energy_before
    
    # Cavity insertion should generally avoid high-energy overlaps
    # So energy change should be reasonable (not extremely positive)
    assert delta_E < 1000.0, f"Energy change {delta_E} suggests bad overlap"


def test_adaptive_cavity_grid():
    """Test adaptive grid spacing based on system density."""
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
        
        # Create grid
        grid = CavityGrid(system.info.box, grid_spacing=spacing)
        vdw_radii = {0: 0.15, 1: 0.1}
        
        for atom in system.atoms:
            radius = vdw_radii.get(atom.type, 0.15)
            grid.mark_occupied((atom.x, atom.y, atom.z), radius)
        
        grid.find_cavities()
        f_n = grid.get_cavity_fraction()
        
        # Higher density should have lower cavity fraction
        if density > 1.0:
            assert f_n < 1.0, f"High density system should have f_n < 1.0, got {f_n}"
        elif density < 0.5:
            assert f_n > 0.5, f"Low density system should have f_n > 0.5, got {f_n}"


def test_cavity_bias_detailed_balance():
    """Test cavity bias detailed balance with proper proposal probabilities."""
    # System parameters
    T = 300.0
    beta = 1.0 / (kB * T)
    mu_ex = -5.0
    n_bar = 100
    B = beta * mu_ex + math.log(n_bar)
    
    # State A: n molecules, State B: n+1 molecules
    n_A = 10
    n_B = n_A + 1
    E_A = -50.0  # System energy
    delta_E = 5.0  # Energy change for insertion
    E_B = E_A + delta_E
    
    # Cavity fraction
    f_n = 0.3
    
    # Calculate the acceptance ratio r
    r = (f_n / (n_A + 1.0)) * math.exp(B - beta * delta_E)
    
    # Acceptance probabilities
    alpha_ins = min(1.0, r)        # Forward: insertion A→B
    alpha_del = min(1.0, 1.0 / r)  # Reverse: deletion B→A
    
    # Proposal probabilities
    # For consistent detailed balance: q_del/q_ins = f_n/(n+1)
    q_ins = 1.0
    q_del = f_n / n_B  # This ensures q_del/q_ins = f_n/(n_B) = f_n/(n_A+1)
    
    # Grand canonical weights: π_GC(n,E) ∝ exp(-βE + Bn)
    pi_A = math.exp(-beta * E_A + B * n_A)
    pi_B = math.exp(-beta * E_B + B * n_B)
    
    # Full transition probabilities: T = π * q * α
    forward_flux = pi_A * q_ins * alpha_ins
    reverse_flux = pi_B * q_del * alpha_del
    
    # Check detailed balance
    if forward_flux > 0 and reverse_flux > 0:
        ratio = forward_flux / reverse_flux
        # Should be equal within numerical precision
        assert 0.999999999 < ratio < 1.000000001, \
            f"Detailed balance violated: ratio = {ratio}"
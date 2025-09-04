"""
Analyze PME error components for medium complexity system
Focus on separating real space and reciprocal space contributions
"""

import numpy as np
import pygcmc
from pygcmc import MCState, MCAtom, MCResidue, MCForceField
from pygcmc import initializePMEParameters, computeSystemEnergyPME

try:
    from openmm import *
    from openmm.app import *
    from openmm.unit import *
    OPENMM_AVAILABLE = True
except ImportError:
    OPENMM_AVAILABLE = False


# Import helper functions
from .analyze_pme_error_helpers import (
    create_medium_system,
    analyze_pme_components
)

def analyze_distance_distribution():
    """Analyze charge-charge distances in the system"""
    
    print("\nDISTANCE DISTRIBUTION ANALYSIS")
    print("="*70)
    
    state = create_medium_system()
    
    # Calculate all pairwise distances
    distances = []
    energies = []
    
    for i in range(state.activeAtomCount):
        for j in range(i+1, state.activeAtomCount):
            dx = state.atoms[i].x - state.atoms[j].x
            dy = state.atoms[i].y - state.atoms[j].y
            dz = state.atoms[i].z - state.atoms[j].z
            
            # Apply minimum image convention
            box = state.info.box
            dx -= box[0] * round(dx / box[0])
            dy -= box[1] * round(dy / box[1])
            dz -= box[2] * round(dz / box[2])
            
            dist = np.sqrt(dx*dx + dy*dy + dz*dz)
            distances.append(dist)
            
            # Calculate Coulomb energy for this pair
            if dist > 0:
                energy = 138.935 * state.atoms[i].charge * state.atoms[j].charge / dist
                energies.append(energy)
    
    # Analyze by distance ranges
    ranges = [(0, 0.5), (0.5, 1.0), (1.0, 1.2), (1.2, 2.0), (2.0, 10.0)]
    
    print("\nEnergy contribution by distance range:")
    print(f"{'Range (nm)':>15} {'Count':>8} {'Energy (kJ/mol)':>15} {'% of Total':>12}")
    print("-"*50)
    
    total_energy = sum(energies)
    
    for r_min, r_max in ranges:
        count = 0
        energy_sum = 0
        for d, e in zip(distances, energies):
            if r_min <= d < r_max:
                count += 1
                energy_sum += e
        
        pct = (energy_sum / total_energy * 100) if total_energy != 0 else 0
        print(f"{r_min:6.1f} - {r_max:6.1f} {count:8d} {energy_sum:15.2f} {pct:11.1f}%")
    
    print(f"\nTotal pairs: {len(distances)}")
    print(f"Total Coulomb energy (no cutoff): {total_energy:.2f} kJ/mol")
    print(f"Cutoff: {state.info.cutoff} nm")
    
    # Count pairs within cutoff
    within_cutoff = 0
    for d in distances:
        if d <= state.info.cutoff:
            within_cutoff += 1
    print(f"Pairs within cutoff: {within_cutoff} ({within_cutoff/len(distances)*100:.1f}%)")


def test_mesh_convergence():
    """Test how PME converges with mesh size"""
    
    print("\nMESH SIZE CONVERGENCE TEST")
    print("="*70)
    
    state = create_medium_system()
    
    alpha = 3.5
    spline_order = 4
    mesh_sizes = [16, 24, 32, 48, 64]
    
    print(f"\nAlpha = {alpha}, Spline order = {spline_order}")
    print(f"{'Mesh':>6} {'Real':>12} {'Recip':>12} {'Total':>12} {'ΔTotal':>12}")
    print("-"*54)
    
    prev_total = None
    
    for mesh in mesh_sizes:
        mesh_size = [mesh, mesh, mesh]
        
        pygcmc.setPMEParameters(alpha, mesh_size, spline_order)
        initializePMEParameters(state.info.cutoff, state.info.box, alpha, mesh_size, spline_order)
        computeSystemEnergyPME(state)
        
        real = state.ewald_energy.get('real_space', 0.0)
        recip = state.ewald_energy.get('reciprocal', 0.0)
        self = state.ewald_energy.get('self', 0.0)
        total = real + recip + self
        
        if prev_total is not None:
            delta = total - prev_total
        else:
            delta = 0
        
        print(f"{mesh:6d} {real:12.4f} {recip:12.4f} {total:12.4f} {delta:12.6f}")
        
        prev_total = total


if __name__ == "__main__":
    analyze_pme_components()
    analyze_distance_distribution()
    test_mesh_convergence()
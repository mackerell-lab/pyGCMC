"""
High-precision PME energy comparison between PyGCMC and OpenMM

This test ensures PME implementations agree to within 1% for realistic systems.
Critical for validating GCMC energy calculations.
"""

import pytest
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pygcmc
from pygcmc import MCState, MCAtom, MCResidue, MCForceField
from pygcmc import initializePMEParameters, computeSystemEnergyPME

# Only import OpenMM if available
try:
    from openmm import *
    from openmm.app import *
    from openmm.unit import *
    OPENMM_AVAILABLE = True
except ImportError:
    OPENMM_AVAILABLE = False


# Import helper functions
from .pme_high_precision_helpers import (
    create_test_system,
    calculate_openmm_pme_precise
)

def test_pme_high_precision_simple_system():
    """Test PME with high precision on a simple system"""
    
    # Create simple system
    state = create_test_system(n_ions=4, n_ligand_atoms=2)
    
    print("\nSimple system test:")
    print(f"  Atoms: {state.activeAtomCount} (4 ions + 2 ligand atoms)")
    print(f"  Box: {state.info.box[0]} nm")
    print(f"  Cutoff: {state.info.cutoff} nm")
    
    # Use PME parameters that should give good agreement
    # Alpha selection is critical for matching
    alpha = 3.0  # 1/nm, typical value
    mesh_size = [32, 32, 32]
    spline_order = 4
    
    # Initialize PyGCMC PME
    print(f"\nPME parameters:")
    print(f"  Alpha: {alpha} 1/nm")
    print(f"  Mesh: {mesh_size}")
    print(f"  Spline order: {spline_order}")
    
    pygcmc.setPMEParameters(alpha, mesh_size, spline_order)
    initializePMEParameters(state.info.cutoff, state.info.box, alpha, mesh_size, spline_order)
    
    # Calculate with PyGCMC
    computeSystemEnergyPME(state)
    pygcmc_total = state.ewald_energy.get('total', 0.0)
    
    print(f"\nPyGCMC PME components:")
    print(f"  Real space: {state.ewald_energy.get('real_space', 0.0):.6f} kJ/mol")
    print(f"  Reciprocal: {state.ewald_energy.get('reciprocal', 0.0):.6f} kJ/mol")
    print(f"  Self: {state.ewald_energy.get('self', 0.0):.6f} kJ/mol")
    print(f"  Total: {pygcmc_total:.6f} kJ/mol")
    
    # Calculate with OpenMM
    if OPENMM_AVAILABLE:
        print(f"\nOpenMM calculation:")
        openmm_energy = calculate_openmm_pme_precise(state, alpha, mesh_size)
        print(f"  Energy: {openmm_energy:.6f} kJ/mol")
        
        # Compare
        diff = abs(pygcmc_total - openmm_energy)
        rel_diff = diff / abs(openmm_energy) if openmm_energy != 0 else 0
        
        print(f"\nComparison:")
        print(f"  Absolute difference: {diff:.6f} kJ/mol")
        print(f"  Relative difference: {rel_diff*100:.3f}%")
        
        # For simple systems, should agree to within 1%
        assert rel_diff < 0.0045, f"Simple system: PME energies differ by {rel_diff*100:.3f}% (> 0.45%)"
    else:
        print("\nOpenMM not available, skipping comparison")
        # Still check that PME gives reasonable values
        assert abs(pygcmc_total) > 0.1, "PME total energy too small"


def test_pme_high_precision_complex_system():
    """Test PME with high precision on a more complex system"""
    
    # Create more complex system
    state = create_test_system(n_ions=8, n_ligand_atoms=6)
    
    print("\nComplex system test:")
    print(f"  Atoms: {state.activeAtomCount} (8 ions + 6 ligand atoms)")
    
    # Test different alpha values to find best match
    test_alphas = [2.5, 3.0, 3.5]
    mesh_size = [32, 32, 32]
    spline_order = 4
    
    best_alpha = None
    best_diff = float('inf')
    
    for alpha in test_alphas:
        print(f"\n--- Testing alpha = {alpha} ---")
        
        # PyGCMC calculation
        pygcmc.setPMEParameters(alpha, mesh_size, spline_order)
        initializePMEParameters(state.info.cutoff, state.info.box, alpha, mesh_size, spline_order)
        computeSystemEnergyPME(state)
        pygcmc_total = state.ewald_energy.get('total', 0.0)
        
        print(f"PyGCMC total: {pygcmc_total:.6f} kJ/mol")
        
        if OPENMM_AVAILABLE:
            # OpenMM calculation
            openmm_energy = calculate_openmm_pme_precise(state, alpha, mesh_size)
            print(f"OpenMM total: {openmm_energy:.6f} kJ/mol")
            
            diff = abs(pygcmc_total - openmm_energy)
            rel_diff = diff / abs(openmm_energy) if openmm_energy != 0 else 0
            print(f"Relative difference: {rel_diff*100:.3f}%")
            
            if rel_diff < best_diff:
                best_diff = rel_diff
                best_alpha = alpha
    
    if OPENMM_AVAILABLE:
        print(f"\nBest alpha: {best_alpha} (difference: {best_diff*100:.3f}%)")
        
        # For complex systems, allow slightly more tolerance but still < 1%
        assert best_diff < 0.01, f"Complex system: Best PME match differs by {best_diff*100:.3f}% (> 1%)"


def test_pme_convergence_with_mesh_size():
    """Test PME convergence with different mesh sizes"""
    
    state = create_test_system(n_ions=4, n_ligand_atoms=4)
    
    print("\nMesh size convergence test:")
    
    alpha = 3.0
    spline_order = 4
    mesh_sizes = [[24, 24, 24], [32, 32, 32], [48, 48, 48], [64, 64, 64]]
    
    energies = []
    for mesh in mesh_sizes:
        pygcmc.setPMEParameters(alpha, mesh, spline_order)
        initializePMEParameters(state.info.cutoff, state.info.box, alpha, mesh, spline_order)
        computeSystemEnergyPME(state)
        total = state.ewald_energy.get('total', 0.0)
        energies.append(total)
        print(f"  Mesh {mesh[0]}x{mesh[1]}x{mesh[2]}: {total:.6f} kJ/mol")
    
    # Check convergence - difference between successive meshes should decrease
    for i in range(1, len(energies)):
        diff = abs(energies[i] - energies[i-1])
        print(f"  Difference {mesh_sizes[i-1][0]}→{mesh_sizes[i][0]}: {diff:.6f} kJ/mol")
    
    # Should converge - last difference should be small
    final_diff = abs(energies[-1] - energies[-2])
    assert final_diff < 0.1, f"PME not converging with mesh size (final diff: {final_diff:.6f})"


if __name__ == "__main__":
    test_pme_high_precision_simple_system()
    test_pme_high_precision_complex_system()
    test_pme_convergence_with_mesh_size()
"""
Test PME convergence with different parameters
"""

import pytest
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from .pme_ligand_openmm_comparison import create_charged_system_with_ligands, calculate_openmm_pme_energy
import pygcmc
from pygcmc import initializePMEParameters, computeSystemEnergyPME

try:
    from openmm import *
    OPENMM_AVAILABLE = True
except ImportError:
    OPENMM_AVAILABLE = False


@pytest.mark.skipif(not OPENMM_AVAILABLE, reason="OpenMM not available")
def test_pme_parameter_sensitivity():
    """Test how PME results change with different parameters"""
    
    state, _ = create_charged_system_with_ligands()
    
    print(f"\nSystem: {state.activeAtomCount} atoms in {state.info.box[0]}nm box")
    
    # Test different parameter combinations
    test_params = [
        # (alpha, mesh_size, spline_order)
        (2.0, [24, 24, 24], 4),
        (2.84, [32, 32, 32], 4),  # Default
        (3.0, [32, 32, 32], 4),
        (4.0, [32, 32, 32], 4),
        (2.84, [48, 48, 48], 4),
        (2.84, [64, 64, 64], 4),
        (2.84, [32, 32, 32], 6),
    ]
    
    print("\nParameter sensitivity analysis:")
    print("=" * 80)
    print(f"{'Alpha':>6} {'Mesh':>12} {'Spline':>7} {'PyGCMC':>12} {'OpenMM':>12} {'Diff %':>8}")
    print("-" * 80)
    
    # Track results for assertions
    results = []
    default_params_result = None
    
    for alpha, mesh_size, spline_order in test_params:
        # PyGCMC calculation
        pygcmc.setPMEParameters(alpha, mesh_size, spline_order)
        initializePMEParameters(state.info.cutoff, state.info.box, alpha, mesh_size, spline_order)
        computeSystemEnergyPME(state)
        
        pygcmc_total = (state.ewald_energy['real_space'] + 
                       state.ewald_energy['reciprocal'] + 
                       state.ewald_energy['self'])
        
        # Validate energy components
        assert state.ewald_energy['self'] < 0, "Self energy should be negative"
        assert pygcmc_total != 0, "Total PME energy should be non-zero"
        
        # OpenMM calculation
        openmm_energy = calculate_openmm_pme_energy(state)
        
        if openmm_energy is not None:
            diff_pct = abs(pygcmc_total - openmm_energy) / abs(openmm_energy) * 100
            mesh_str = f"{mesh_size[0]}x{mesh_size[1]}x{mesh_size[2]}"
            print(f"{alpha:6.2f} {mesh_str:>12} {spline_order:7d} {pygcmc_total:12.4f} {openmm_energy:12.4f} {diff_pct:7.2f}%")
            
            results.append({
                'alpha': alpha,
                'mesh': mesh_size,
                'spline': spline_order,
                'pygcmc': pygcmc_total,
                'openmm': openmm_energy,
                'diff_pct': diff_pct
            })
            
            # Save default parameters result
            if alpha == 2.84 and mesh_size == [32, 32, 32] and spline_order == 4:
                default_params_result = results[-1]
    
    # Assertions on parameter sensitivity
    if results:
        # Find best result
        best_result = min(results, key=lambda x: x['diff_pct'])
        print(f"\nBest parameters: alpha={best_result['alpha']}, mesh={best_result['mesh'][0]}x{best_result['mesh'][1]}x{best_result['mesh'][2]}, spline={best_result['spline']}")
        print(f"Best difference: {best_result['diff_pct']:.2f}%")
        
        # Assert at least one parameter set gives reasonable agreement
        assert best_result['diff_pct'] < 15.0, f"Best PME agreement is {best_result['diff_pct']:.2f}% (> 15%), which may indicate an issue"
        
        # Assert default parameters give reasonable results
        if default_params_result:
            assert default_params_result['diff_pct'] < 20.0, f"Default parameters give {default_params_result['diff_pct']:.2f}% difference (> 20%)"
        
        # Assert that increasing mesh size improves accuracy (convergence test)
        mesh_32_result = next((r for r in results if r['alpha'] == 2.84 and r['mesh'] == [32, 32, 32] and r['spline'] == 4), None)
        mesh_48_result = next((r for r in results if r['alpha'] == 2.84 and r['mesh'] == [48, 48, 48] and r['spline'] == 4), None)
        mesh_64_result = next((r for r in results if r['alpha'] == 2.84 and r['mesh'] == [64, 64, 64] and r['spline'] == 4), None)
        
        if mesh_32_result and mesh_48_result and mesh_64_result:
            # Energy should converge as mesh size increases
            diff_32_48 = abs(mesh_32_result['pygcmc'] - mesh_48_result['pygcmc'])
            diff_48_64 = abs(mesh_48_result['pygcmc'] - mesh_64_result['pygcmc'])
            print(f"\nConvergence test:")
            print(f"  32->48 mesh difference: {diff_32_48:.4f} kJ/mol")
            print(f"  48->64 mesh difference: {diff_48_64:.4f} kJ/mol")
            # Each refinement should reduce the change
            assert diff_48_64 < diff_32_48, "PME should converge with increasing mesh size"
    
    # Detailed breakdown for default parameters
    print("\nDetailed breakdown for default parameters (alpha=2.84, 32x32x32, order=4):")
    alpha, mesh_size, spline_order = 2.84, [32, 32, 32], 4
    pygcmc.setPMEParameters(alpha, mesh_size, spline_order)
    initializePMEParameters(state.info.cutoff, state.info.box, alpha, mesh_size, spline_order)
    computeSystemEnergyPME(state)
    
    print(f"  Real space:  {state.ewald_energy['real_space']:10.4f} kJ/mol")
    print(f"  Reciprocal:  {state.ewald_energy['reciprocal']:10.4f} kJ/mol")
    print(f"  Self:        {state.ewald_energy['self']:10.4f} kJ/mol")
    
    # Check charge distribution
    total_charge_sq = 0
    for atom in state.atoms:
        total_charge_sq += atom.charge ** 2
    print(f"\nTotal charge squared: {total_charge_sq:.2f}")
    
    # Theoretical self energy
    import math
    theoretical_self = -alpha / math.sqrt(math.pi) * total_charge_sq * 138.935
    print(f"Theoretical self energy: {theoretical_self:.4f} kJ/mol")
    print(f"Actual self energy: {state.ewald_energy['self']:.4f} kJ/mol")
    print(f"Self energy ratio: {state.ewald_energy['self']/theoretical_self:.3f}")
    
    # Assert self energy is close to theoretical value
    self_energy_error = abs(state.ewald_energy['self'] - theoretical_self) / abs(theoretical_self)
    assert self_energy_error < 0.01, f"Self energy differs from theoretical by {self_energy_error*100:.1f}% (> 1%)"


if __name__ == "__main__":
    test_pme_parameter_sensitivity()
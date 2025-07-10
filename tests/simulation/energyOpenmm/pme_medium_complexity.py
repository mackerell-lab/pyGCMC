"""
Test PME accuracy for medium complexity system (multiple ions + ligand, no water)

This test bridges the gap between simple ion systems (1-2% error) and 
complex water-containing systems (15% error) to understand PME accuracy scaling.
"""

import pytest
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pygcmc
from pygcmc import initializePMEParameters, computeSystemEnergyPME

from simulation.energyOpenmm.pme_medium_complexity_helpers import (
    create_medium_complexity_system,
    calculate_openmm_energy_medium,
    OPENMM_AVAILABLE
)


@pytest.mark.skipif(not OPENMM_AVAILABLE, reason="OpenMM not available")
def test_pme_medium_complexity():
    """Test PME accuracy for medium complexity system"""
    
    state = create_medium_complexity_system()
    
    print("\nMedium Complexity System:")
    print("=" * 70)
    print(f"  {state.activeAtomCount} atoms: 4 Na+, 4 Cl-, 12 ligand atoms")
    print(f"  Box: {state.info.box[0]} x {state.info.box[1]} x {state.info.box[2]} nm")
    print(f"  Cutoff: {state.info.cutoff} nm")
    
    # Test range of alpha values
    alphas = [4.0, 5.0, 5.6, 6.0, 7.0, 8.0]
    mesh_size = [32, 32, 32]
    spline_order = 4
    
    print(f"\n{'Alpha':>6} {'PyGCMC Energy':>15} {'OpenMM Energy':>15} {'Difference':>12} {'Rel Error %':>12}")
    print("-" * 70)
    
    errors = []
    
    for alpha in alphas:
        # Initialize PME
        pygcmc.setPMEParameters(alpha, mesh_size, spline_order)
        initializePMEParameters(state.info.cutoff, state.info.box, alpha, mesh_size, spline_order)
        
        # Calculate with PyGCMC
        computeSystemEnergyPME(state)
        # Get PME energy components
        pygcmc_energy = (state.ewald_energy['real_space'] + 
                        state.ewald_energy['reciprocal'] + 
                        state.ewald_energy['self'])
        
        # Calculate with OpenMM
        openmm_energy = calculate_openmm_energy_medium(state, alpha)
        
        if openmm_energy is not None:
            diff = pygcmc_energy - openmm_energy
            rel_error = abs(diff) / abs(openmm_energy) * 100 if openmm_energy != 0 else 0
            errors.append(rel_error)
            
            print(f"{alpha:6.1f} {pygcmc_energy:15.4f} {openmm_energy:15.4f} {diff:12.4f} {rel_error:12.2f}")
        else:
            print(f"{alpha:6.1f} {pygcmc_energy:15.4f} {'N/A':>15} {'N/A':>12} {'N/A':>12}")
    
    if errors:
        avg_error = sum(errors) / len(errors)
        print(f"\nAverage relative error: {avg_error:.2f}%")
        
        # Check that error is reasonable for medium complexity
        assert avg_error < 5.0, f"Average PME error {avg_error:.2f}% exceeds 5% threshold"
        
        # Best case should be under 2%
        min_error = min(errors)
        print(f"Best relative error: {min_error:.2f}%")
        assert min_error < 2.0, f"Best PME error {min_error:.2f}% exceeds 2% threshold"


if __name__ == "__main__":
    test_pme_medium_complexity()
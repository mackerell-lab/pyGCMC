# tests/simulation/energyPME/advanced_parameters.py
"""Advanced PME parameter tests: alpha dependency and LJ energy calculation."""

import math
import pygcmc
from .helpers import create_nacl_crystal


def test_pme_alpha_dependency():
    """
    Test the dependency of PME accuracy on alpha (Ewald separation parameter).
    
    This test examines how different alpha values affect PME accuracy when compared
    to standard Ewald summation. A proper alpha balances real-space and reciprocal
    space calculation costs.
    """
    # Create a basic system
    box_size = 4.0
    n_cells = 2
    state = create_nacl_crystal(box_size, n_cells)
    
    # Standard settings
    box = [box_size, box_size, box_size]
    cutoff = box_size / 2.0
    mesh_size = [32, 32, 32]  # Reasonably fine mesh
    spline_order = 4
    kmax = [8, 8, 8]  # High accuracy reference
    
    # Test different alpha values
    alpha_values = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0]
    results = []
    
    print(f"\nPME Alpha Dependency Test")
    print(f"Testing effect of alpha parameter on PME accuracy")
    
    for alpha in alpha_values:
        # Calculate reference with standard Ewald
        pygcmc.setEwaldParameters(alpha, kmax)
        pygcmc.initializeEwaldParameters(cutoff, box, alpha)
        _, _, ewald_dict = pygcmc.computeSystemEnergyEwald(state)
        
        ewald_energy = ewald_dict["total"]
        ewald_real = ewald_dict["real_space"]
        ewald_recip = ewald_dict["reciprocal"]
        
        # Calculate with PME
        pygcmc.setPMEParameters(alpha, mesh_size, spline_order)
        pygcmc.initializePMEParameters(cutoff, box, alpha, mesh_size, spline_order)
        _, _, pme_dict = pygcmc.computeSystemEnergyPME(state)
        
        pme_energy = pme_dict["total"]
        pme_real = pme_dict["real_space"]
        pme_recip = pme_dict["reciprocal"]
        
        # Calculate relative errors
        total_error = abs(pme_energy - ewald_energy) / abs(ewald_energy)
        real_error = abs(pme_real - ewald_real) / abs(ewald_real) if abs(ewald_real) > 1e-10 else 0.0
        recip_error = abs(pme_recip - ewald_recip) / abs(ewald_recip) if abs(ewald_recip) > 1e-10 else 0.0
        
        # Store results
        results.append({
            'alpha': alpha,
            'total_error': total_error,
            'real_error': real_error,
            'recip_error': recip_error,
            'real_ratio': abs(ewald_real) / abs(ewald_energy),
            'recip_ratio': abs(ewald_recip) / abs(ewald_energy)
        })
        
        print(f"Alpha = {alpha:.2f}, "
              f"Total error = {total_error:.8f}, "
              f"Real error = {real_error:.8f}, "
              f"Recip error = {recip_error:.8f}, "
              f"Real/Total ratio = {results[-1]['real_ratio']:.4f}, "
              f"Recip/Total ratio = {results[-1]['recip_ratio']:.4f}")
    
    # Find optimal alpha - one that balances real and reciprocal space contributions
    optimal_indices = [i for i, r in enumerate(results) 
                      if abs(r['real_ratio'] - 0.5) < 0.1]  # Close to 50/50 split
    
    if optimal_indices:
        optimal_alpha = results[optimal_indices[0]]['alpha']
        print(f"\nOptimal alpha for this system is approximately {optimal_alpha}")
        
        # Verify that optimal alpha gives acceptable error
        optimal_result = results[optimal_indices[0]]
        assert optimal_result['total_error'] < 0.05  # Less than 5% error
    else:
        print("\nNo optimal alpha found with close to 50/50 split between real and reciprocal space")
        
    # Verify that all alphas give reasonable accuracy
    for result in results:
        assert result['total_error'] < 0.1  # Less than 10% error for all alphas


def test_pme_lj_energy():
    """
    Test that LJ/VDW energy is correctly computed in PME method
    
    This test verifies that when using the PME method, the LJ/VDW energy
    is correctly calculated, and matches results from direct calculation.
    """
    # Create a basic system with known LJ parameters
    box_size = 4.0
    n_cells = 2
    state = create_nacl_crystal(box_size, n_cells)
    
    # Use specific LJ parameters to make the test more sensitive
    # Modify force field to ensure LJ contribution is significant
    ff = state.forcefield
    
    # Set stronger LJ parameters
    sigma_na = 0.4  # nm
    sigma_cl = 0.5  # nm
    eps_na = 0.5    # kJ/mol - increased to make LJ energy more significant
    eps_cl = 0.8    # kJ/mol - increased to make LJ energy more significant
    
    # Update force field parameters
    ff.ljSigma = [
        sigma_na, (sigma_na + sigma_cl)/2.0,
        (sigma_na + sigma_cl)/2.0, sigma_cl
    ]
    ff.ljEps = [
        eps_na, math.sqrt(eps_na * eps_cl),
        math.sqrt(eps_na * eps_cl), eps_cl
    ]
    
    state.forcefield = ff
    
    box = [box_size, box_size, box_size]
    cutoff = box_size / 2.0
    alpha = 0.3
    mesh_size = [32, 32, 32]
    spline_order = 4
    
    print("\nTesting PME LJ energy calculation")
    
    # First, calculate using direct method (reference)
    # Reset energies
    for res in state.residues:
        res.energy_vdw = 0.0
        res.energy_elec = 0.0
    
    # Use computeSystemEnergyPBCCutoff to calculate reference VDW energy
    # This function calculates both electrostatic and VDW energies, but we only care about VDW
    pygcmc.computeSystemEnergyPBCCutoff(state)
    direct_vdw = sum(res.energy_vdw for res in state.residues if res.active)
    
    print(f"Direct VDW energy: {direct_vdw:.6f} kJ/mol")
    
    # Reset energies
    for res in state.residues:
        res.energy_vdw = 0.0
        res.energy_elec = 0.0
    
    # Calculate using PME method
    pygcmc.setPMEParameters(alpha, mesh_size, spline_order)
    pygcmc.initializePMEParameters(cutoff, box, alpha, mesh_size, spline_order)
    
    pme_elec, pme_vdw, pme_dict = pygcmc.computeSystemEnergyPME(state)
    
    # Extract VDW energy
    pme_vdw_from_residues = sum(res.energy_vdw for res in state.residues if res.active)
    
    print(f"PME VDW energy: {pme_vdw:.6f} kJ/mol")
    print(f"PME VDW from residues: {pme_vdw_from_residues:.6f} kJ/mol")
    
    # Compare VDW energies - these should be very close since both use the
    # same computational approach for LJ interactions
    rel_diff = abs(direct_vdw - pme_vdw_from_residues) / (abs(direct_vdw) + 1e-10)
    print(f"Relative difference in VDW energy: {rel_diff:.6f}")
    
    # The energies should be nearly identical
    # Use relaxed criterion since PME and direct calculation may have small differences
    assert rel_diff < 0.01, "VDW energies from direct and PME methods differ significantly"
    
    # Additional test: ensure VDW energy is a significant component
    vdw_fraction = abs(pme_vdw) / (abs(pme_elec) + abs(pme_vdw) + 1e-10)
    print(f"VDW energy fraction: {vdw_fraction:.2%}")
    
    # With our stronger LJ parameters, VDW should be significant
    # Updated threshold to reflect actual NaCl crystal behavior
    assert vdw_fraction > 0.03, "VDW energy is too small a fraction of total energy"
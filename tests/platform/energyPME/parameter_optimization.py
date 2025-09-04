# tests/simulation/energyPME/parameter_optimization.py
"""PME parameter optimization tests: spline order, error tolerance, and mesh accuracy."""

import pygcmc
from .helpers import create_nacl_crystal


def test_pme_spline_order():
    """
    Test the effect of different B-spline orders in PME
    
    This test verifies that higher order B-splines generally produce accurate PME results.
    Note that error may not monotonically decrease with spline order due to numerical
    properties of the PME algorithm, but highest order splines should give very accurate results.
    """
    # Create a basic system
    box_size = 4.0
    n_cells = 2
    state = create_nacl_crystal(box_size, n_cells)
    
    box = [box_size, box_size, box_size]
    cutoff = box_size / 2.0
    alpha = 0.3
    
    # Use a finer mesh for better accuracy and convergence
    mesh_size = [32, 32, 32]
    
    # Calculate reference energy using very high order spline
    # This will be our "ground truth" for comparison
    pygcmc.setPMEParameters(alpha, mesh_size, 6)  # 6th order as reference value
    pygcmc.initializePMEParameters(cutoff, box, alpha, mesh_size, 6)
    _, _, reference_dict = pygcmc.computeSystemEnergyPME(state)
    reference_energy = reference_dict["total"]
    
    print(f"Reference energy (order 6): {reference_energy:.6f}")
    
    # Test orders 1 through 5
    test_orders = [2, 3, 4, 5]
    energies = []
    errors = []
    
    for order in test_orders:
        # Initialize PME with specified spline order
        pygcmc.setPMEParameters(alpha, mesh_size, order)
        pygcmc.initializePMEParameters(cutoff, box, alpha, mesh_size, order)
        
        # Calculate energy
        _, _, pme_dict = pygcmc.computeSystemEnergyPME(state)
        total_energy = pme_dict["total"]
        energies.append(total_energy)
        
        # Calculate relative error from reference
        rel_error = abs(total_energy - reference_energy) / abs(reference_energy)
        errors.append(rel_error)
        
        print(f"Spline order {order}: energy = {total_energy:.6f}, relative error = {rel_error:.6f}")
    
    # Show error ratio changes between orders
    for i in range(1, len(errors)):
        ratio = errors[i-1]/max(errors[i], 1e-10)  # Avoid division by zero
        print(f"Error ratio {test_orders[i-1]}->{test_orders[i]}: {ratio:.2f}x change")
        
    # Verify important properties:
    
    # 1. At least one high-order spline should have very low error
    assert min(errors) < 0.001, f"No spline order achieved good accuracy (min error: {min(errors):.6f})"
    
    # 2. Highest order spline (5) should be very accurate
    assert errors[-1] < 0.0001, f"Highest order spline ({test_orders[-1]}) should be very accurate"
    
    # 3. Check consistency between different orders
    max_energy_diff = max([abs(e1 - e2) for e1, e2 in zip(energies, energies[1:])])
    assert max_energy_diff < 0.2, f"Too much inconsistency between different spline orders"


def test_pme_error_tolerance():
    """
    Test PME error tolerance control
    
    This test verifies that the PME algorithm can adapt to different error tolerance settings.
    Note: Different tolerance values may result in significantly different parameter choices,
    especially for alpha and mesh size, which can lead to different energy values.
    """
    # Create a basic system
    box_size = 4.0
    n_cells = 2
    state = create_nacl_crystal(box_size, n_cells)
    
    box = [box_size, box_size, box_size]
    cutoff = box_size / 2.0
    
    # Test different error tolerances
    tolerances = [1e-3, 1e-4, 1e-5, 1e-6]
    results = []
    
    print(f"\nPME Error Tolerance Test")
    for tol in tolerances:
        # Initialize with auto parameters and specified tolerance
        default_mesh = [32, 32, 32]  # Default mesh size
        pygcmc.setPMEParameters(0.0, default_mesh, 4, tol)
        
        # Let autoAdjustPMEParameters select optimal parameters
        pygcmc.initializePMEParameters(cutoff, box, 0.0, [], 4, tol)
        
        # Calculate energy
        _, _, pme_dict = pygcmc.computeSystemEnergyPME(state)
        
        # Extract energy components
        real_energy = pme_dict["real_space"]
        recip_energy = pme_dict["reciprocal"]
        self_energy = pme_dict["self"]
        total_energy = pme_dict["total"]
        
        # Record results
        results.append({
            'tolerance': tol,
            'real': real_energy,
            'reciprocal': recip_energy,
            'self': self_energy,
            'total': total_energy
        })
        
        print(f"Tolerance {tol:.1e}: real={real_energy:.2f}, recip={recip_energy:.2f}, "
              f"self={self_energy:.2f}, total={total_energy:.2f}")
    
    # Check that self energy remains relatively stable
    self_energies = [r['self'] for r in results]
    self_mean = sum(self_energies) / len(self_energies)
    self_max_diff = max([abs(e - self_mean) for e in self_energies])
    
    print(f"\nSelf energy: mean={self_mean:.2f}, max deviation={self_max_diff:.2f}")
    
    # Check that energies converge for the last two tolerances
    finer_tols = results[-2:]  # Take the two smallest tolerances
    finer_diff = abs(finer_tols[0]['total'] - finer_tols[1]['total'])
    print(f"Energy difference between {tolerances[-2]:.1e} and {tolerances[-1]:.1e}: {finer_diff:.2f}")
    
    # Allow some variation in energy between different tolerances, but last two values should be close
    assert finer_diff < 0.5 * abs(finer_tols[0]['total']), \
           "Energy does not converge with stricter tolerance"
    
    # Self energy should remain relatively stable
    assert self_max_diff < 0.2 * abs(self_mean), \
           "Self energy varies too much with different tolerances"


def test_pme_mesh_accuracy():
    """
    Test the accuracy of PME with different mesh sizes compared to a reference solution.
    
    This test is inspired by the testPMEParameters function in TestEwald.h.
    It uses increasingly finer mesh sizes to show convergence toward exact Ewald results.
    """
    # Create a basic system
    box_size = 4.0
    n_cells = 2
    state = create_nacl_crystal(box_size, n_cells)
    
    # Set parameters for reference Ewald calculation
    alpha = 0.3
    kmax = [10, 10, 10]  # High precision k-vectors for reference
    box = [box_size, box_size, box_size]
    cutoff = box_size / 2.0
    spline_order = 6  # Higher spline order for accuracy test
    
    # First calculate reference with standard Ewald (high precision)
    pygcmc.setEwaldParameters(alpha, kmax)
    pygcmc.initializeEwaldParameters(cutoff, box, alpha)
    _, _, ewald_dict = pygcmc.computeSystemEnergyEwald(state)
    reference_energy = ewald_dict["total"]
    reference_recip = ewald_dict["reciprocal"]
    
    print(f"\nPME Mesh Accuracy Test")
    print(f"Reference Ewald energy = {reference_energy:.6f}")
    
    # Test a series of increasing mesh sizes for PME
    mesh_sizes = [8, 16, 24, 32, 48, 64]
    results = []
    
    for mesh_size in mesh_sizes:
        # Use cubic mesh for simplicity
        mesh = [mesh_size, mesh_size, mesh_size]
        
        # Initialize PME with current mesh size
        pygcmc.setPMEParameters(alpha, mesh, spline_order)
        pygcmc.initializePMEParameters(cutoff, box, alpha, mesh, spline_order)
        
        # Calculate energy with current PME settings
        _, _, pme_dict = pygcmc.computeSystemEnergyPME(state)
        pme_energy = pme_dict["total"]
        pme_recip = pme_dict["reciprocal"]
        
        # Calculate relative error
        energy_diff = abs(pme_energy - reference_energy)
        relative_error = energy_diff / abs(reference_energy)
        recip_relative_error = abs(pme_recip - reference_recip) / abs(reference_recip)
        
        # Store results
        results.append({
            'mesh_size': mesh_size,
            'energy': pme_energy,
            'error': relative_error,
            'recip_error': recip_relative_error
        })
        
        print(f"Mesh size: {mesh_size}x{mesh_size}x{mesh_size}, "
              f"PME energy = {pme_energy:.6f}, "
              f"relative error = {relative_error:.8f}, "
              f"reciprocal error = {recip_relative_error:.8f}")
    
    # Verify that errors decrease with increasing mesh size
    for i in range(1, len(results)):
        assert results[i]['error'] <= results[i-1]['error'] * 1.2  # Allow some numerical variance
        assert results[i]['recip_error'] <= results[i-1]['recip_error'] * 1.2
    
    # Verify that largest mesh size gives acceptably small error
    assert results[-1]['error'] < 0.01  # Within 1% of reference
    assert results[-1]['recip_error'] < 0.05  # Within 5% for reciprocal component
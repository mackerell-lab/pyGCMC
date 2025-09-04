# tests/simulation/energyPME/basic_tests.py
"""Basic PME tests: initialization, PME vs Ewald comparison, and movement energy."""

import pygcmc
from pygcmc import MCMovementResidueInfo
from .helpers import create_nacl_crystal


def test_pme_initialization():
    """
    Test PME parameter initialization
    """
    # Create a basic system
    box_size = 4.0
    n_cells = 2
    state = create_nacl_crystal(box_size, n_cells)
    
    # Set system box size
    box = [box_size, box_size, box_size]
    cutoff = box_size / 2.0
    
    # Explicitly set PME parameters, without using auto-adjust
    alpha = 0.3
    mesh_size = [32, 32, 32]
    spline_order = 4
    
    # Ensure using explicit call with all parameters
    print(f"Initializing PME with parameters: alpha={alpha}, mesh_size={mesh_size}, cutoff={cutoff}")
    pygcmc.setPMEParameters(alpha, mesh_size, spline_order)
    pygcmc.initializePMEParameters(cutoff, box, alpha, mesh_size, spline_order)
    
    # Output parameters for debugging
    print(f"PME parameters set. Ready to compute energy.")
    
    # Calculate energy
    elec_energy, vdw_energy, ewald_dict = pygcmc.computeSystemEnergyPME(state)
    
    # Check that energy components are reasonable
    assert ewald_dict["real_space"] != 0.0
    assert ewald_dict["reciprocal"] != 0.0
    assert ewald_dict["self"] != 0.0
    assert ewald_dict["total"] != 0.0
    
    print(f"PME energy components: real_space={ewald_dict['real_space']:.6f}, "
          f"reciprocal={ewald_dict['reciprocal']:.6f}, "
          f"self={ewald_dict['self']:.6f}, "
          f"total={ewald_dict['total']:.6f}")
    
    # Check that total energy is the sum of components
    assert abs(ewald_dict["total"] - (ewald_dict["real_space"] + 
                                    ewald_dict["reciprocal"] + 
                                    ewald_dict["self"] + vdw_energy)) < 1e-6


def test_pme_vs_ewald():
    """
    Compare PME and Ewald methods for a simple system
    
    This test verifies that our OpenMM-based PME implementation produces
    results consistent with the standard Ewald method.
    """
    # Create a basic system
    box_size = 4.0
    n_cells = 2
    state = create_nacl_crystal(box_size, n_cells)
    
    # Set parameters for both methods
    alpha = 0.3
    kmax = [5, 5, 5]
    # Use a reasonable mesh size for accuracy
    mesh_size = [32, 32, 32]
    spline_order = 4  # 4th order B-splines typically have good precision/performance balance
    
    box = [box_size, box_size, box_size]
    cutoff = box_size / 2.0
    
    # First calculate with Ewald
    pygcmc.setEwaldParameters(alpha, kmax)
    pygcmc.initializeEwaldParameters(cutoff, box, alpha)
    ewald_elec, ewald_vdw, ewald_dict = pygcmc.computeSystemEnergyEwald(state)
    
    ewald_real = ewald_dict["real_space"]
    ewald_recip = ewald_dict["reciprocal"]
    ewald_self = ewald_dict["self"]
    ewald_total = ewald_dict["total"]
    
    # Then calculate with PME
    pygcmc.setPMEParameters(alpha, mesh_size, spline_order)
    pygcmc.initializePMEParameters(cutoff, box, alpha, mesh_size, spline_order)
    
    # Add debug output
    print(f"\nPME vs Ewald Comparison:")
    print(f"Alpha: {alpha}, Mesh Size: {mesh_size}, Spline Order: {spline_order}")
    print(f"Box: {box}, Cutoff: {cutoff}")
    
    pme_elec, pme_vdw, pme_dict = pygcmc.computeSystemEnergyPME(state)
    
    pme_real = pme_dict["real_space"]
    pme_recip = pme_dict["reciprocal"]
    pme_self = pme_dict["self"]
    pme_total = pme_dict["total"]
    
    # Calculate relative errors
    real_rel_error = abs(ewald_real - pme_real) / abs(ewald_real)
    recip_rel_error = abs(ewald_recip - pme_recip) / abs(ewald_recip)
    self_rel_error = abs(ewald_self - pme_self) / abs(ewald_self)
    total_rel_error = abs(ewald_total - pme_total) / abs(ewald_total)
    
    # Print comparison
    print("Energy comparison:")
    print(f"              Real Space      Reciprocal     Self           Total")
    print(f"Ewald:        {ewald_real:.6f}    {ewald_recip:.6f}    {ewald_self:.6f}    {ewald_total:.6f}")
    print(f"PME:          {pme_real:.6f}    {pme_recip:.6f}    {pme_self:.6f}    {pme_total:.6f}")
    print(f"Rel. Error:   {real_rel_error:.6f}    {recip_rel_error:.6f}    {self_rel_error:.6f}    {total_rel_error:.6f}")
    
    # Verify results meet accuracy requirements
    # Real-space and self energies should be nearly identical
    assert real_rel_error < 1e-6, "Real-space energies don't match"
    assert self_rel_error < 1e-6, "Self energies don't match"
    
    # Reciprocal space energy should be within 0.5% for a mesh of 32³
    assert recip_rel_error < 0.005, f"Reciprocal space error ({recip_rel_error:.2%}) exceeds threshold"
    
    # Total energy should be within 0.2%
    assert total_rel_error < 0.002, f"Total energy error ({total_rel_error:.2%}) exceeds threshold"


def test_pme_movement_energy():
    """
    Test movement energy calculation with PME
    """
    # Create a basic system
    box_size = 4.0
    n_cells = 2
    state = create_nacl_crystal(box_size, n_cells)
    
    # Set PME parameters
    box = [box_size, box_size, box_size]
    cutoff = box_size / 2.0
    alpha = 0.3
    mesh_size = [32, 32, 32]
    
    pygcmc.setPMEParameters(alpha, mesh_size)
    pygcmc.initializePMEParameters(cutoff, box, alpha)
    
    # Set up movement residues
    movement_info = MCMovementResidueInfo()
    movement_info.startIndex = 0
    movement_info.activeCount = 2  # Move 2 residues
    state.movementResidues.append(movement_info)
    
    # First calculate full system energy
    _, _, system_pme_dict = pygcmc.computeSystemEnergyPME(state)
    full_energy = system_pme_dict["total"]
    
    # Then calculate movement energy
    _, _, movement_pme_dict = pygcmc.computeMovementEnergyPME(state)
    movement_energy = movement_pme_dict["total"]
    
    print(f"Full system energy: {full_energy:.6f}")
    print(f"Movement energy: {movement_energy:.6f}")
    
    # Movement energy should be a subset of the full energy
    assert abs(movement_energy) <= abs(full_energy) * 1.1  # Allow some numerical variance
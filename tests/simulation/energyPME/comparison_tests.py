# tests/simulation/energyPME/comparison_tests.py
"""PME comparison tests: Ewald vs PME with different systems."""

import random
import pygcmc
from pygcmc import MCState, MCAtom, MCResidue, MCForceField


def test_ewald_vs_pme_comparison():
    """
    Test consistency between Ewald and PME methods, replicating testEwaldVsPME from pme.cpp
    
    This test compares the energy calculations between a standard PME implementation
    and a high precision PME implementation that mimics the Ewald method.
    """
    print("\nRunning test_ewald_vs_pme_comparison (replicating testEwaldVsPME from pme.cpp)...")
    
    # Create an amorphous salt system of random particles
    num_particles = 100
    box_size = 3.0
    cutoff = 1.0
    
    # Create the state
    state = MCState()
    state.info.box = [box_size, box_size, box_size]
    state.info.cutoff = cutoff
    state.info.setTemperature(300.0)
    
    # Set a simple force field
    ff = MCForceField()
    ff.numTotalTypes = 1  # Single atom type for simplicity
    ff.ljSigma = [0.3]  # Dummy LJ parameters - size matches numTotalTypes^2
    ff.ljEps = [0.0]  # Set LJ epsilon to zero to eliminate LJ interactions
    
    state.forcefield = ff
    
    # Create atoms with random positions with fixed seed
    random.seed(98765)  # Same seed as in the C++ version
    
    atoms = []
    residues = []
    
    print(f"Creating a system with {num_particles} randomly positioned particles...")
    
    # Generate random positions and set alternating charges
    for i in range(num_particles):
        atom = MCAtom()
        # Set atom position
        atom.x = random.uniform(0, box_size)
        atom.y = random.uniform(0, box_size)
        atom.z = random.uniform(0, box_size)
        
        # Set charge - alternating positive and negative charges
        if i < num_particles // 2:
            atom.charge = 1.0  # Na+
            atom.type = 0
        else:
            atom.charge = -1.0  # Cl-
            atom.type = 1
        
        atoms.append(atom)
        
        # Create residues - one residue for every two atoms
        if i % 2 == 0:
            res = MCResidue()
            res.atomStart = i
            res.atomCount = 2 if i < num_particles-1 else 1  # Handle last atom
            res.active = True
            res.fixed = False
            residues.append(res)
    
    # Set active atom and residue counts
    state.atoms = atoms
    state.residues = residues  # Set residues
    state.activeAtomCount = len(atoms)  # Add active atom count
    state.activeResidueCount = len(residues)  # Add active residue count
    
    # Update force field parameters to support two atom types
    ff.numTotalTypes = 2  # Change to 2 atom types (Na+ and Cl-)
    ff.ljSigma = [0.3, 0.3, 0.3, 0.3]  # Expand to 2x2 matrix
    ff.ljEps = [0.0, 0.0, 0.0, 0.0]    # Expand to 2x2 matrix
    
    state.forcefield = ff
    
    # Calculate total charge to verify system is neutral
    total_charge = sum(atom.charge for atom in atoms)
    print(f"Total system charge: {total_charge}")
    assert abs(total_charge) < 1e-10, "System must be charge neutral"
    
    # Method 1: PME with standard parameters
    alpha_pme = 2.5 / cutoff
    grid_size_pme = 32  # Power of 2, as in C++ version
    spline_order_pme = 5  # As in C++ version
    
    # Method 2: PME with high precision parameters to mimic Ewald
    alpha_ewald = alpha_pme  # Same alpha
    grid_size_ewald = 64  # Higher precision grid, still power of 2
    spline_order_ewald = 6  # Higher order interpolation, as in C++ version
    
    # Setup box for both methods
    box = state.info.box
    
    # Calculate PME energy with standard parameters
    pme_mesh_standard = [grid_size_pme, grid_size_pme, grid_size_pme]
    pygcmc.setPMEParameters(alpha_pme, pme_mesh_standard, spline_order_pme)
    pygcmc.initializePMEParameters(cutoff, box, alpha_pme, pme_mesh_standard, spline_order_pme)
    
    _, _, pme_dict = pygcmc.computeSystemEnergyPME(state)
    
    # Extract energy components for standard PME
    energy_pme_real = pme_dict["real_space"]
    energy_pme_recip = pme_dict["reciprocal"]
    energy_pme_self = pme_dict["self"]
    energy_pme_total = pme_dict["total"]
    
    # Calculate with high precision parameters (Ewald-like)
    pme_mesh_ewald = [grid_size_ewald, grid_size_ewald, grid_size_ewald]
    pygcmc.setPMEParameters(alpha_ewald, pme_mesh_ewald, spline_order_ewald)
    pygcmc.initializePMEParameters(cutoff, box, alpha_ewald, pme_mesh_ewald, spline_order_ewald)
    
    _, _, ewald_dict = pygcmc.computeSystemEnergyPME(state)
    
    # Extract energy components for Ewald-like PME
    energy_ewald_real = ewald_dict["real_space"]
    energy_ewald_recip = ewald_dict["reciprocal"]
    energy_ewald_self = ewald_dict["self"]
    energy_ewald_total = ewald_dict["total"]
    
    # Compare results as in the C++ version
    print("PME Energy (standard parameters):")
    print(f"  Total:      {energy_pme_total:.6f}")
    print(f"  Real space: {energy_pme_real:.6f}")
    print(f"  Reciprocal: {energy_pme_recip:.6f}")
    print(f"  Self:       {energy_pme_self:.6f}")
    
    print("Ewald-like Energy (high precision parameters):")
    print(f"  Total:      {energy_ewald_total:.6f}")
    print(f"  Real space: {energy_ewald_real:.6f}")
    print(f"  Reciprocal: {energy_ewald_recip:.6f}")
    print(f"  Self:       {energy_ewald_self:.6f}")
    
    # Calculate relative differences
    abs_ewald_total = abs(energy_ewald_total)
    scale_factor = abs_ewald_total if abs_ewald_total > 1.0 else 1.0
    
    rel_diff_total = abs(energy_pme_total - energy_ewald_total) / scale_factor
    rel_diff_real = abs(energy_pme_real - energy_ewald_real) / (abs(energy_ewald_real) if abs(energy_ewald_real) > 1.0 else 1.0)
    rel_diff_recip = abs(energy_pme_recip - energy_ewald_recip) / (abs(energy_ewald_recip) if abs(energy_ewald_recip) > 1.0 else 1.0)
    
    print(f"Relative energy difference (total): {rel_diff_total:.6f}")
    print(f"Real space relative difference:     {rel_diff_real:.6f}")
    print(f"Reciprocal space relative difference: {rel_diff_recip:.6f}")
    
    # Assert that the difference is small - tolerate up to 0.2% difference
    assert rel_diff_total < 0.002, f"Total energy differs by more than 0.2%: {rel_diff_total:.6f}"
    
    print("test_ewald_vs_pme_comparison completed successfully")


def test_ewald_vs_pme_random():
    """
    Compare Ewald and PME methods using a random particle system.
    
    This test is modeled after the testEwaldVsPME function in pme.cpp.
    It creates a random distribution of charged particles and compares energy 
    calculations between standard Ewald and PME methods.
    """
    # Create a random particle system similar to pme.cpp's testEwaldVsPME
    num_particles = 100
    box_size = 3.0
    cutoff = 1.0
    
    # Create a new state with random positions
    state = MCState()
    state.info.box = [box_size, box_size, box_size]
    state.info.setTemperature(300.0)
    state.info.cutoff = cutoff
    
    # Set force field parameters for simple ions
    ff = MCForceField()
    ff.numTotalTypes = 2  # Positive and negative ions
    
    # Simple LJ parameters (identical for simplicity)
    sigma = 0.3  # nm
    epsilon = 0.1  # kJ/mol
    ff.ljSigma = [sigma, sigma, sigma, sigma]
    ff.ljEps = [epsilon, epsilon, epsilon, epsilon]
    
    state.forcefield = ff
    
    # Generate random particles with alternating charges
    random.seed(98765)  # Use same seed as in pme.cpp
    
    atoms = []
    residues = []
    
    print(f"\nCreating random system with {num_particles} particles...")
    
    for i in range(num_particles):
        # Create a new atom
        atom = MCAtom()
        
        # Random position within box
        atom.x = random.random() * box_size
        atom.y = random.random() * box_size
        atom.z = random.random() * box_size
        
        # Alternating charges
        if i < num_particles/2:
            atom.charge = 1.0  # Positive
            atom.type = 0
        else:
            atom.charge = -1.0  # Negative
            atom.type = 1
        
        atoms.append(atom)
        
        # Create one residue per atom pair
        if i % 2 == 0:
            res = MCResidue()
            res.atomStart = i
            res.atomCount = 2 if i < num_particles-1 else 1
            res.active = True
            res.fixed = False
            residues.append(res)
    
    state.atoms = atoms
    state.residues = residues
    state.activeAtomCount = len(atoms)
    state.activeResidueCount = len(residues)
    
    # System parameters
    box = [box_size, box_size, box_size]
    
    # Set parameters for standard Ewald calculation (reference)
    alpha_ewald = 2.5 / cutoff  # Same as in pme.cpp
    kmax_ewald = [8, 8, 8]  # Higher precision for reference
    
    # Calculate with standard Ewald
    pygcmc.setEwaldParameters(alpha_ewald, kmax_ewald)
    pygcmc.initializeEwaldParameters(cutoff, box, alpha_ewald)
    ewald_elec, ewald_vdw, ewald_dict = pygcmc.computeSystemEnergyEwald(state)
    
    ewald_real = ewald_dict["real_space"]
    ewald_recip = ewald_dict["reciprocal"]
    ewald_self = ewald_dict["self"]
    ewald_total = ewald_dict["total"]
    
    # Set parameters for standard PME calculation
    alpha_pme = alpha_ewald  # Use same alpha for comparison
    mesh_size_pme = [32, 32, 32]  # Standard mesh size
    spline_order_pme = 5  # Same as in pme.cpp
    
    # Calculate with PME
    pygcmc.setPMEParameters(alpha_pme, mesh_size_pme, spline_order_pme)
    pygcmc.initializePMEParameters(cutoff, box, alpha_pme, mesh_size_pme, spline_order_pme)
    pme_elec, pme_vdw, pme_dict = pygcmc.computeSystemEnergyPME(state)
    
    pme_real = pme_dict["real_space"]
    pme_recip = pme_dict["reciprocal"]
    pme_self = pme_dict["self"]
    pme_total = pme_dict["total"]
    
    # Calculate relative errors for PME vs Ewald
    real_rel_error = abs(ewald_real - pme_real) / max(abs(ewald_real), 1.0)
    recip_rel_error = abs(ewald_recip - pme_recip) / max(abs(ewald_recip), 1.0)
    self_rel_error = abs(ewald_self - pme_self) / max(abs(ewald_self), 1.0)
    total_rel_error = abs(ewald_total - pme_total) / max(abs(ewald_total), 1.0)
    
    # Print comparison
    print("\nEwald vs PME Comparison for Random System:")
    print(f"Parameters: Alpha={alpha_pme:.4f}, Cutoff={cutoff:.2f}nm, Box={box_size:.2f}nm")
    print(f"PME: Mesh={mesh_size_pme}, Spline Order={spline_order_pme}")
    print(f"Ewald: kmax={kmax_ewald}")
    
    print("\nEnergy Comparison:")
    print(f"                  Real Space      Reciprocal     Self           Total")
    print(f"Ewald:            {ewald_real:.6f}    {ewald_recip:.6f}    {ewald_self:.6f}    {ewald_total:.6f}")
    print(f"PME:              {pme_real:.6f}    {pme_recip:.6f}    {pme_self:.6f}    {pme_total:.6f}")
    
    print("\nRelative Errors (vs Ewald):")
    print(f"                  Real Space      Reciprocal     Self           Total")
    print(f"PME:              {real_rel_error:.6f}    {recip_rel_error:.6f}    {self_rel_error:.6f}    {total_rel_error:.6f}")
    
    # Verify results meet accuracy requirements
    assert real_rel_error < 0.01, "Real-space energies don't match"
    assert self_rel_error < 0.01, "Self energies don't match"
    assert recip_rel_error < 0.05, f"Reciprocal space error ({recip_rel_error:.2%}) exceeds threshold"
    assert total_rel_error < 0.05, f"Total energy error ({total_rel_error:.2%}) exceeds threshold"
    
    print(f"\nTest passed: PME accuracy is within expected thresholds")
    print(f"PME total energy error: {total_rel_error:.2%}")
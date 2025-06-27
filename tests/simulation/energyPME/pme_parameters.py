# tests/simulation/energyPME/pme_parameters.py
"""PME parameters influence test."""

import math
import random
import pygcmc
from pygcmc import MCState, MCAtom, MCResidue, MCForceField


def test_pme_parameters():
    """
    Test PME parameters influence on accuracy, replicating testPMEParameters from pme.cpp
    
    This test examines how different PME parameters (alpha, grid size, interpolation order,
    and dielectric constant) affect the calculation accuracy when compared to a high precision
    reference calculation.
    """
    print("\nRunning test_pme_parameters (replicating testPMEParameters from pme.cpp)...")
    
    # Create a simple system of random charges similar to the C++ version
    num_particles = 51
    box_size = 4.7
    
    state = MCState()
    
    # Set box size and cutoff
    state.info.box = [box_size, box_size, box_size]
    state.info.setTemperature(300.0)
    state.info.cutoff = 2.0  # Same as in C++ version
    
    # Set a simple force field
    ff = MCForceField()
    ff.numTotalTypes = 1  # Single atom type for simplicity
    ff.ljSigma = [0.3]  # Dummy LJ parameters - size matches numTotalTypes^2
    ff.ljEps = [0.0]  # Set LJ epsilon to zero to eliminate LJ interactions
    
    state.forcefield = ff
    
    # Create atoms with random positions
    # Use a fixed seed for reproducibility, matching the C++ version
    import random
    random.seed(12345)
    
    atoms = []
    residues = []
    
    print(f"Creating a system with {num_particles} randomly positioned particles...")
    
    # Generate random positions and distribute charges uniformly between -1 and +1
    for i in range(num_particles):
        atom = MCAtom()
        
        # Random position within box
        atom.x = random.random() * box_size
        atom.y = random.random() * box_size
        atom.z = random.random() * box_size
        
        # Distribute charges uniformly between -1 and +1
        atom.charge = -1.0 + i * 2.0 / (num_particles - 1)
        atom.type = 0
        
        atoms.append(atom)
        
        # Create one residue per atom for simplicity
        if i % 1 == 0:  # Every atom gets its own residue
            res = MCResidue()
            res.atomStart = i
            res.atomCount = 1
            res.active = True
            res.fixed = False
            residues.append(res)
    
    state.atoms = atoms
    state.residues = residues
    state.activeAtomCount = len(atoms)
    state.activeResidueCount = len(residues)
    
    # System parameters
    box = [box_size, box_size, box_size]
    cutoff = 2.0  # Same as in C++ version
    
    # First compute energy with a high precision setting
    # Use high alpha (more real space) and finer grid
    print("\nCalculating high precision reference energy...")
    alpha_high = 2.5
    mesh_size_high = [64, 64, 64]  # High precision mesh
    spline_order_high = 5
    
    # Set high precision PME parameters
    pygcmc.setPMEParameters(alpha_high, mesh_size_high, spline_order_high)
    pygcmc.initializePMEParameters(cutoff, box, alpha_high, mesh_size_high, spline_order_high)
    
    # Calculate high precision reference energy
    high_elec, high_vdw, high_dict = pygcmc.computeSystemEnergyPME(state)
    energy_high = high_dict["total"]
    
    print(f"Reference energy (high precision): {energy_high:.6f}")
    
    # Test 1: Low alpha (more reciprocal space)
    print("\nTest 1: Testing with lower alpha (more reciprocal space work)...")
    alpha_low = 1.5
    mesh_size_low_alpha = [64, 64, 64]  # Keep high mesh
    
    pygcmc.setPMEParameters(alpha_low, mesh_size_low_alpha, spline_order_high)
    pygcmc.initializePMEParameters(cutoff, box, alpha_low, mesh_size_low_alpha, spline_order_high)
    
    low_alpha_elec, low_alpha_vdw, low_alpha_dict = pygcmc.computeSystemEnergyPME(state)
    energy_low_alpha = low_alpha_dict["total"]
    
    rel_error_low_alpha = abs(energy_high - energy_low_alpha) / max(abs(energy_high), 1.0)
    
    print(f"Energy with low alpha: {energy_low_alpha:.6f}")
    print(f"Relative error with low alpha: {rel_error_low_alpha:.8f}")
    
    # Test 2: Coarser grid
    print("\nTest 2: Testing with coarser grid...")
    mesh_size_coarse = [32, 32, 32]  # Coarser grid
    
    pygcmc.setPMEParameters(alpha_high, mesh_size_coarse, spline_order_high)
    pygcmc.initializePMEParameters(cutoff, box, alpha_high, mesh_size_coarse, spline_order_high)
    
    coarse_elec, coarse_vdw, coarse_dict = pygcmc.computeSystemEnergyPME(state)
    energy_coarse = coarse_dict["total"]
    
    rel_error_coarse = abs(energy_high - energy_coarse) / max(abs(energy_high), 1.0)
    
    print(f"Energy with coarse grid: {energy_coarse:.6f}")
    print(f"Relative error with coarse grid: {rel_error_coarse:.8f}")
    
    # Test 3: Lower interpolation order
    print("\nTest 3: Testing with lower interpolation order...")
    spline_order_low = 3  # Lower order
    
    pygcmc.setPMEParameters(alpha_high, mesh_size_high, spline_order_low)
    pygcmc.initializePMEParameters(cutoff, box, alpha_high, mesh_size_high, spline_order_low)
    
    low_order_elec, low_order_vdw, low_order_dict = pygcmc.computeSystemEnergyPME(state)
    energy_low_order = low_order_dict["total"]
    
    rel_error_low_order = abs(energy_high - energy_low_order) / max(abs(energy_high), 1.0)
    
    print(f"Energy with low order: {energy_low_order:.6f}")
    print(f"Relative error with low order: {rel_error_low_order:.8f}")
    
    # Test 4: Different dielectric constant
    # Note: In Python/PyGCMC we may need to handle dielectric differently than the C++ version
    # Here we're simulating it by scaling the charges
    print("\nTest 4: Testing with different dielectric constant (simulated by scaling charges)...")
    
    # Create a copy of the state with scaled charges to simulate dielectric = 2.0
    dielectric_state = MCState()
    dielectric_state.info = state.info
    dielectric_state.forcefield = state.forcefield
    
    dielectric_atoms = []
    for atom in state.atoms:
        dielectric_atom = MCAtom()
        dielectric_atom.x = atom.x
        dielectric_atom.y = atom.y
        dielectric_atom.z = atom.z
        dielectric_atom.charge = atom.charge / math.sqrt(2.0)  # Scale charge to simulate ε = 2.0
        dielectric_atom.type = atom.type
        dielectric_atoms.append(dielectric_atom)
    
    dielectric_state.atoms = dielectric_atoms
    dielectric_state.residues = state.residues
    dielectric_state.activeAtomCount = state.activeAtomCount
    dielectric_state.activeResidueCount = state.activeResidueCount
    
    # Calculate with "dielectric" system
    pygcmc.setPMEParameters(alpha_high, mesh_size_high, spline_order_high)  # Use high precision settings
    pygcmc.initializePMEParameters(cutoff, box, alpha_high, mesh_size_high, spline_order_high)
    
    dielectric_elec, dielectric_vdw, dielectric_dict = pygcmc.computeSystemEnergyPME(dielectric_state)
    energy_dielectric = dielectric_dict["total"]
    
    ratio_with_dielectric = energy_high / energy_dielectric
    
    print(f"Energy with dielectric 2.0: {energy_dielectric:.6f}")
    print(f"Energy ratio with dielectric 2.0 (should be ~2.0): {ratio_with_dielectric:.6f}")
    
    # Verify results meet reasonable accuracy requirements
    print("\nVerifying accuracy thresholds...")
    
    # Alpha test - should be reasonably accurate
    assert rel_error_low_alpha < 0.05, f"Low alpha relative error too high: {rel_error_low_alpha:.6f}"
    
    # Grid test - should be reasonably accurate
    assert rel_error_coarse < 0.05, f"Coarse grid relative error too high: {rel_error_coarse:.6f}"
    
    # Spline order test - should be reasonably accurate
    assert rel_error_low_order < 0.10, f"Low spline order relative error too high: {rel_error_low_order:.6f}"
    
    # Dielectric test - energy should scale approximately with dielectric
    assert abs(ratio_with_dielectric - 2.0) < 0.2, f"Dielectric scaling off, ratio: {ratio_with_dielectric:.6f}"
    
    print("test_pme_parameters completed successfully")
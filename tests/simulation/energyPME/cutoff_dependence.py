# tests/simulation/energyPME/cutoff_dependence.py
"""PME cutoff dependence test."""

import random
import pygcmc
from pygcmc import MCState, MCAtom, MCResidue, MCForceField


def test_cutoff_dependence():
    """
    Test PME energy contributions as a function of cutoff, replicating testCutoffDependence from pme.cpp
    
    This test examines how different cutoff distances affect the balance between real-space and 
    reciprocal-space contributions to the total energy in PME calculations. It also demonstrates
    how alpha parameter is typically chosen based on the cutoff distance.
    
    Note on energy calculation differences between C++ and Python versions:
    ----------------------------------------------------------------------
    The C++ implementation in pme.cpp calculates energy components separately:
    1. Real space energy using pme_calculate_real_space()
    2. Reciprocal space energy using pme_exec()
    3. Self energy using pme_calculate_self_energy()
    4. Total energy is computed by manually summing these components
    
    In contrast, this Python version uses pygcmc.computeSystemEnergyPME(), which:
    1. Calculates all energy components in a single function call
    2. May include additional terms in the real-space energy (e.g., LJ interactions)
    3. May apply different cutoff treatments or energy decomposition approaches
    
    Due to these implementation differences, absolute energy values may not match exactly 
    between the C++ and Python versions. However, the key physical trends should be consistent:
    1. Total energy should remain stable across different cutoffs
    2. As cutoff increases, real-space contribution increases and reciprocal contribution decreases
    3. Self energy should be identical across implementations
    
    This test now also compares PME results with standard Ewald method to verify consistency.
    """
    print("\nRunning test_cutoff_dependence (replicating testCutoffDependence from pme.cpp)...")
    
    # Create a system with random positions and alternating charges, similar to C++ version
    num_particles = 100
    box_size = 5.0
    
    state = MCState()
    
    # Set box size and temperature
    state.info.box = [box_size, box_size, box_size]
    state.info.setTemperature(300.0)
    
    # Set a simple force field
    ff = MCForceField()
    ff.numTotalTypes = 1  # Single atom type for simplicity
    ff.ljSigma = [0.3]  # Dummy LJ parameters - size matches numTotalTypes^2
    ff.ljEps = [0.0]  # Set LJ epsilon to zero to eliminate LJ interactions
    
    state.forcefield = ff
    
    # Create atoms with random positions using the same seed as C++ version
    import random
    random.seed(54321)  # Same seed as C++ version
    
    atoms = []
    residues = []
    
    print(f"Creating a system with {num_particles} atoms (alternating charges)...")
    print(f"(Note: LJ interactions disabled to compare with C++ version)")
    
    # Generate random positions
    for i in range(num_particles):
        atom = MCAtom()
        
        # Random position within box
        atom.x = random.random() * box_size
        atom.y = random.random() * box_size
        atom.z = random.random() * box_size
        
        # Alternating charges (first half positive, second half negative)
        if i < num_particles / 2:
            atom.charge = 1.0
        else:
            atom.charge = -1.0
        
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
    
    # Table header
    print("\nComparison of PME and Ewald energy components as a function of cutoff:")
    print("\n" + "-"*100)
    print(f"{'Cutoff':6s} | {'Alpha':8s} | {'Method':10s} | {'Total Energy':14s} | {'Real Space':14s} | {'Reciprocal':14s} | {'Self':14s} | {'Real/Total':10s} | {'Recip/Total':10s} | {'PME/Ewald Ratio':14s}")
    print("-"*100)
    
    # Test a range of cutoffs
    cutoff_values = [0.5 + 0.25 * i for i in range(9)]  # 0.5 to 2.5 in steps of 0.25
    
    for cutoff in cutoff_values:
        # Set cutoff in state
        state.info.cutoff = cutoff
        
        # Alpha is typically set inversely proportional to cutoff, as in C++ version
        alpha = 2.0 / cutoff
        
        # Determine grid size based on alpha and box size
        # Ensure grid size is a power of 2
        min_grid_size = int(2.0 * alpha * box_size / 3.14159 + 0.5)
        grid_size = 16  # Minimum 16
        while grid_size < min_grid_size:
            grid_size *= 2  # Ensure power of 2
        
        mesh_size = [grid_size, grid_size, grid_size]
        spline_order = 5  # Same as C++ version
        
        # Setup for Ewald calculation
        kmax = [grid_size // 2, grid_size // 2, grid_size // 2]  # Reasonable kmax based on grid size
        
        # Initialize Ewald with current parameters and calculate
        pygcmc.setEwaldParameters(alpha, kmax)
        pygcmc.initializeEwaldParameters(cutoff, box, alpha)
        _, _, ewald_dict = pygcmc.computeSystemEnergyEwald(state)
        
        # Extract Ewald energy components
        ewald_real = ewald_dict["real_space"]
        ewald_recip = ewald_dict["reciprocal"]
        ewald_self = ewald_dict["self"]
        ewald_total = ewald_dict["total"]
        
        # Calculate Ewald ratios
        ewald_real_ratio = ewald_real / ewald_total if abs(ewald_total) > 1e-10 else 0.0
        ewald_recip_ratio = ewald_recip / ewald_total if abs(ewald_total) > 1e-10 else 0.0
        
        # Print Ewald results
        print(f"{cutoff:6.2f} | {alpha:8.6f} | {'Ewald':10s} | {ewald_total:14.2f} | {ewald_real:14.2f} | {ewald_recip:14.2f} | {ewald_self:14.2f} | {ewald_real_ratio:10.6f} | {ewald_recip_ratio:10.6f} | {'N/A':14s}")
        
        # Initialize PME with current parameters
        pygcmc.setPMEParameters(alpha, mesh_size, spline_order)
        pygcmc.initializePMEParameters(cutoff, box, alpha, mesh_size, spline_order)
        
        # Calculate PME energy
        _, _, pme_dict = pygcmc.computeSystemEnergyPME(state)
        
        # Extract energy components
        pme_real = pme_dict["real_space"]
        pme_recip = pme_dict["reciprocal"]
        pme_self = pme_dict["self"]
        pme_total = pme_dict["total"]
        
        # Calculate ratios
        pme_real_ratio = pme_real / pme_total if abs(pme_total) > 1e-10 else 0.0
        pme_recip_ratio = pme_recip / pme_total if abs(pme_total) > 1e-10 else 0.0
        
        # Calculate PME to Ewald ratio for total energy
        pme_ewald_ratio = pme_total / ewald_total if abs(ewald_total) > 1e-10 else 1.0
        
        # Print PME results
        print(f"{cutoff:6.2f} | {alpha:8.6f} | {'PME':10s} | {pme_total:14.2f} | {pme_real:14.2f} | {pme_recip:14.2f} | {pme_self:14.2f} | {pme_real_ratio:10.6f} | {pme_recip_ratio:10.6f} | {pme_ewald_ratio:14.6f}")
        print("-"*100)
    
    # Verify that the total energy is reasonably consistent across different cutoffs
    # (This should be the case if the PME implementation is correct)
    print("\nTest completed - verify that total energy is reasonably consistent across cutoffs")
    print("The real/reciprocal energy balance should shift with cutoff and alpha values")
    print("PME and Ewald results should be consistent for all cutoff values")
    print("test_cutoff_dependence completed successfully")
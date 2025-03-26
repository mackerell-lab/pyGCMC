# tests/simulation/test_energy_PGP.py

import pytest
import math
import random
import pygcmc
from pygcmc import MCState, MCInfo, MCAtom, MCResidue, MCForceField, MCMovementResidueInfo
import sys

# Set log level to INFO or lower to ensure detailed log output
# System log settings
pygcmc.System.set_log_level(pygcmc.LogLevel.INFO)
pygcmc.System.set_verbose(True)

# Platform log settings (for log output in energyPGP.cpp)
pygcmc.set_platform_verbose(True)  # Enable platform log output
pygcmc.set_platform_log_level(pygcmc.PlatformLogLevel.INFO)
pygcmc.set_platform_debug_mode(True)  # Enable debug mode for testing

# For more detailed logs, set to DEBUG
# pygcmc.System.set_log_level(pygcmc.LogLevel.DEBUG)
# pygcmc.set_platform_log_level(pygcmc.PlatformLogLevel.DEBUG)

# Ensure output buffer is flushed immediately
sys.stdout.flush()
print("Log level settings completed")
sys.stdout.flush()

# Direct copy of create_nacl_crystal function from test_energy_PME.py
def create_nacl_crystal(box_size, n_cells):
    """
    Create a NaCl crystal model
    
    Args:
        box_size: box size (nm)
        n_cells: number of unit cells in each dimension
    """
    state = MCState()
    
    # Set box size and temperature
    state.info.box = [box_size, box_size, box_size]
    state.info.setTemperature(300.0)  # 300K
    state.info.cutoff = 1.2  # 1.2 nm cutoff
    
    # Set force field parameters
    ff = MCForceField()
    ff.numTotalTypes = 2  # Na+ and Cl-
    
    # LJ parameters (from OPLS-AA force field)
    sigma_na = 0.333  # nm
    sigma_cl = 0.442  # nm
    eps_na = 0.0115  # kJ/mol
    eps_cl = 0.4184  # kJ/mol
    
    # Set LJ parameter matrix
    ff.ljSigma = [
        sigma_na, (sigma_na + sigma_cl)/2.0,
        (sigma_na + sigma_cl)/2.0, sigma_cl
    ]
    ff.ljEps = [
        eps_na, math.sqrt(eps_na * eps_cl),
        math.sqrt(eps_na * eps_cl), eps_cl
    ]
    
    state.forcefield = ff
    
    # NaCl lattice constant (0.564 nm)
    a = 0.564  
    atoms = []
    residues = []
    
    # Create NaCl lattice
    print(f"\nCreating {n_cells}x{n_cells}x{n_cells} NaCl crystal...")
    for i in range(n_cells):
        for j in range(n_cells):
            for k in range(n_cells):
                # Na+ ion
                na = MCAtom()
                na.x = i * a
                na.y = j * a
                na.z = k * a
                na.charge = 1.0
                na.type = 0
                atoms.append(na)
                
                # Cl- ion
                cl = MCAtom()
                cl.x = i * a + a/2
                cl.y = j * a + a/2
                cl.z = k * a + a/2
                cl.charge = -1.0
                cl.type = 1
                atoms.append(cl)
                
                # Create a residue for each ion pair
                res = MCResidue()
                res.atomStart = len(atoms) - 2
                res.atomCount = 2
                res.active = True
                res.fixed = False
                residues.append(res)
                
    print(f"Creation complete, added a total of {len(atoms)} atoms and {len(residues)} residues.")
    state.atoms = atoms
    state.residues = residues
    state.activeAtomCount = len(atoms)
    state.activeResidueCount = len(residues)
    
    return state

def create_long_distance_system(box_size):
    """
    Create a specialized system for testing reciprocal space calculations
    
    Atoms are placed far apart, beyond the real space cutoff distance, 
    so that reciprocal space calculations will dominate
    
    Args:
        box_size: Box size (nm)
    """
    state = MCState()
    
    # Set box size and temperature
    state.info.box = [box_size, box_size, box_size]
    state.info.setTemperature(300.0)  # 300K
    state.info.cutoff = 1.2  # 1.2 nm cutoff
    
    # Set force field parameters
    ff = MCForceField()
    ff.numTotalTypes = 2  # Two ion types
    
    # LJ parameters
    sigma_na = 0.333  # nm
    sigma_cl = 0.442  # nm
    eps_na = 0.0115  # kJ/mol
    eps_cl = 0.4184  # kJ/mol
    
    # Set LJ parameter matrix
    ff.ljSigma = [
        sigma_na, (sigma_na + sigma_cl)/2.0,
        (sigma_na + sigma_cl)/2.0, sigma_cl
    ]
    ff.ljEps = [
        eps_na, math.sqrt(eps_na * eps_cl),
        math.sqrt(eps_na * eps_cl), eps_cl
    ]
    
    state.forcefield = ff
    
    atoms = []
    residues = []
    
    # Create fixed part: two charged ions placed at opposite corners of the box
    print(f"\nCreating system with long-distance interactions...")
    
    # Ion 1: placed at one corner of the box
    ion1 = MCAtom()
    ion1.x = 0.1
    ion1.y = 0.1
    ion1.z = 0.1
    ion1.charge = 1.0
    ion1.type = 0
    atoms.append(ion1)
    
    # Ion 2: placed at the opposite corner
    ion2 = MCAtom()
    ion2.x = box_size - 0.1
    ion2.y = box_size - 0.1
    ion2.z = box_size - 0.1
    ion2.charge = -1.0
    ion2.type = 1
    atoms.append(ion2)
    
    # Create fixed residue
    fixed_res = MCResidue()
    fixed_res.atomStart = 0
    fixed_res.atomCount = 2
    fixed_res.active = True
    fixed_res.fixed = True
    residues.append(fixed_res)
    
    # Create an ion for movement, placed in the center of the box
    ion3 = MCAtom()
    ion3.x = box_size / 2.0
    ion3.y = box_size / 2.0
    ion3.z = box_size / 2.0
    ion3.charge = 1.0
    ion3.type = 0
    atoms.append(ion3)
    
    # Create moving residue
    move_res = MCResidue()
    move_res.atomStart = 2
    move_res.atomCount = 1
    move_res.active = True
    move_res.fixed = False
    residues.append(move_res)
    
    print(f"System creation complete, total of {len(atoms)} atoms and {len(residues)} residues.")
    state.atoms = atoms
    state.residues = residues
    state.activeAtomCount = len(atoms)
    state.activeResidueCount = len(residues)
    
    return state

def test_pgp_parameter_setting():
    """
    Test that PGP parameters can be set properly
    """
    # Just test that the parameter setting doesn't throw an exception
    alpha = 0.29  # 1/nm
    mesh_size = [32, 32, 32]
    potential_grid_size = [16, 16, 16]
    spline_order = 4
    tolerance = 1e-5
    potential_cutoff = 0.5  # nm
    
    # Set the parameters
    pygcmc.setPGPParameters(
        alpha=alpha,
        meshSize=mesh_size,
        potential_cutoff=potential_cutoff,
        potentialGridSize=potential_grid_size,
        splineOrder=spline_order,
        tolerance=tolerance
    )
    
    # This test passes if setPGPParameters doesn't throw an exception
    assert True, "Parameters set successfully"
        
def test_precompute_grid_potential():
    """
    Test precomputing the grid potential
    """
    # Set basic parameters
    box_size = 2.82  # nm, approximately 28.2 Å
    n_cells = 2      # 2x2x2 supercell
    cutoff = 1.0   # nm
    potential_cutoff = 0.5  # nm
    
    # Set PGP parameters
    alpha = 0.29  # 1/nm
    mesh_size = [32, 32, 32]
    potential_grid_size = [16, 16, 16]
    spline_order = 4
    tolerance = 1e-5
    box = [box_size, box_size, box_size]
    
    # Initialize parameters
    pygcmc.setPMEParameters(
        alpha=alpha,
        meshSize=mesh_size,
        splineOrder=spline_order,
        tolerance=tolerance
    )
    
    # Initialize PME parameters - this is a critical step
    pygcmc.initializePMEParameters(cutoff, box, alpha)
    
    pygcmc.setPGPParameters(
        alpha=alpha,
        meshSize=mesh_size,
        potential_cutoff=potential_cutoff,
        potentialGridSize=potential_grid_size,
        splineOrder=spline_order,
        tolerance=tolerance
    )
    
    # Create a model containing fixed and moving parts
    system = create_nacl_crystal(box_size, n_cells)
    
    # Mark half of the residues as fixed
    n_residues = len(system.residues)
    for i in range(0, n_residues, 2):
        system.residues[i].fixed = True
    
    # Precompute grid potential for the fixed part
    pygcmc.precomputeGridPotential(system, fixed_only=True)
    
    # Test successful execution without crashing
    assert True, "Grid potential precomputation succeeded"
    print("Grid potential precomputation succeeded")
        
def test_interpolate_molecule_energy():
    """
    Test interpolating molecule energy from the precomputed grid
    """
    # Set basic parameters
    box_size = 2.82  # nm, approximately 28.2 Å
    n_cells = 2      # 2x2x2 supercell
    cutoff = 1.0   # nm
    potential_cutoff = 0.5  # nm
    
    # Set PGP parameters
    alpha = 0.29  # 1/nm
    mesh_size = [32, 32, 32]
    potential_grid_size = [16, 16, 16]
    spline_order = 4
    tolerance = 1e-5
    box = [box_size, box_size, box_size]
    
    # Initialize parameters
    pygcmc.setPMEParameters(
        alpha=alpha,
        meshSize=mesh_size,
        splineOrder=spline_order,
        tolerance=tolerance
    )
    
    # Initialize PME parameters - this is a critical step
    pygcmc.initializePMEParameters(cutoff, box, alpha)
    
    pygcmc.setPGPParameters(
        alpha=alpha,
        meshSize=mesh_size,
        potential_cutoff=potential_cutoff,
        potentialGridSize=potential_grid_size,
        splineOrder=spline_order,
        tolerance=tolerance
    )
    
    # Create a model containing fixed and moving parts
    system = create_nacl_crystal(box_size, n_cells)
    
    # Mark half of the residues as fixed, half as moving
    n_residues = len(system.residues)
    fixed_residues = []
    moving_residues = []
    
    for i in range(n_residues):
        if i % 2 == 0:
            system.residues[i].fixed = True
            fixed_residues.append(i)
        else:
            system.residues[i].fixed = False
            moving_residues.append(i)
    
    # Set moving residues - using MCMovementResidueInfo correctly
    movement_info = MCMovementResidueInfo()
    movement_info.startIndex = moving_residues[0]  # Index of moving residue
    movement_info.activeCount = len(moving_residues)  # Number of moving residues
    system.movementResidues.append(movement_info)
    
    # Precompute grid potential for the fixed part
    pygcmc.precomputeGridPotential(system, fixed_only=True)
    
    # Calculate interpolated energy - using new function name
    energy = pygcmc.calculateMoleculeEnergy(system)
    
    # Check if energy value is reasonable
    # Note: We don't check the specific energy value here, as the calculation result depends on many factors
    # Just check if the energy value is finite and not exactly zero
    assert math.isfinite(energy), "Energy value should be finite"
    assert energy != 0.0, "Energy value should not be exactly zero"
    
    print(f"Interpolated energy: {energy} kJ/mol")

# Move test functions to module level
def test_compare_pme_pgp_energy():
    """
    Compare whether the energy values calculated by PME and PGP are consistent before and after movement
    
    This test verifies:
    1. In the initial state, energy values calculated by PME and PGP should be the same
    2. After moving the molecule, the energy changes calculated by PME and PGP should be the same
    """
    # Set parameters - ensure PME and PGP use the same parameters
    box_size = 5.0  # nm - use a larger box
    cutoff = 1.0   # nm
    potential_cutoff = 1.0  # nm - same as cutoff
    box = [box_size, box_size, box_size]
    
    alpha = 0.29  # 1/nm
    mesh_size = [32, 32, 32]
    potential_grid_size = [32, 32, 32]  # Use the same grid size as PME for accurate comparison
    spline_order = 4
    tolerance = 1e-5
    
    # Set test timeout to avoid long running
    timeout = 10  # seconds
    
    print("Creating long-distance test system...")
    sys.stdout.flush()
    
    # Create test system - using specially designed long-distance system
    system = create_long_distance_system(box_size)
    
    # Ensure box size is set correctly
    system.info.box = box
    system.info.cutoff = cutoff
    
    print("Setting PME parameters...")
    sys.stdout.flush()
    
    # Set PME and PGP parameters
    pygcmc.setPMEParameters(
        alpha=alpha,
        meshSize=mesh_size,
        splineOrder=spline_order,
        tolerance=tolerance
    )
    
    # Initialize PME parameters - critical step
    pygcmc.initializePMEParameters(cutoff, box, alpha)
    
    print("Setting PGP parameters...")
    sys.stdout.flush()
    
    pygcmc.setPGPParameters(
        alpha=alpha,
        meshSize=mesh_size,
        potential_cutoff=potential_cutoff,
        potentialGridSize=potential_grid_size,
        splineOrder=spline_order,
        tolerance=tolerance
    )
    
    print("Setting up moving residues...")
    sys.stdout.flush()
    
    # Moving residues already set up in create_long_distance_system
    # Just need to prepare movementResidues list
    moving_residues = [1]  # Second residue is the moving residue
    
    # Verify fixed residue information
    fixed_count = sum(1 for res in system.residues if res.fixed)
    print(f"Number of fixed residues: {fixed_count}")
    print(f"Number of moving residues: {len(moving_residues)}")
    sys.stdout.flush()
    
    # Set up moving residue info
    system.movementResidues.clear()
    
    # Create moving residue info
    movement_info = pygcmc.MCMovementResidueInfo()
    movement_info.startIndex = moving_residues[0]  # Index of moving residue
    movement_info.activeCount = 1  # Only one moving residue
    system.movementResidues.append(movement_info)
    
    print(f"Added movement info: startIndex={movement_info.startIndex}, activeCount={movement_info.activeCount}")
    print(f"System has {system.activeResidueCount} active residues and {len(system.movementResidues)} movement residue groups")
    sys.stdout.flush()
    
    # Step 1: Calculate initial system energy using PME
    print("Calculating initial PME energy...")
    sys.stdout.flush()
    initial_pme_result = pygcmc.computeMovementEnergyPME(system)
    initial_pme_energy = initial_pme_result[0]  # PME electrostatic energy
    initial_pme_dict = initial_pme_result[2]  # PME energy details dictionary
    initial_pme_reciprocal = initial_pme_dict['reciprocal']  # Only take reciprocal space part
    print(f"Initial PME energy result: {initial_pme_result}")
    print(f"Initial PME reciprocal energy: {initial_pme_reciprocal}")
    sys.stdout.flush()
    
    # Step 2: Precompute grid potential using PGP and calculate moving residue energy
    print("Precomputing PGP grid potential...")
    sys.stdout.flush()
    pygcmc.precomputeGridPotential(system, fixed_only=True)
    
    print("Calculating initial PGP energy...")
    sys.stdout.flush()
    initial_pgp_energy = pygcmc.calculateMoleculeEnergy(system)
    
    # Print initial energy
    print(f"Initial PME reciprocal energy: {initial_pme_reciprocal}")
    print(f"Initial PGP energy: {initial_pgp_energy}")
    sys.stdout.flush()
    
    # Step 3: Move moving residue (e.g., translate 0.1 nm)
    translation = [0.1, 0.1, 0.1]  # nm
    
    print("Moving residue...")
    sys.stdout.flush()
    for res_idx in moving_residues:
        residue = system.residues[res_idx]
        for atom_idx in range(residue.atomCount):
            atom_index = residue.atomStart + atom_idx
            atom = system.atoms[atom_index]
            atom.x += translation[0]
            atom.y += translation[1]
            atom.z += translation[2]
    
    # Step 4: Calculate energy using PME after movement
    print("Calculating PME after movement energy...")
    sys.stdout.flush()
    moved_pme_result = pygcmc.computeMovementEnergyPME(system)
    moved_pme_energy = moved_pme_result[0]  # PME electrostatic energy
    moved_pme_dict = moved_pme_result[2]  # PME energy details dictionary
    moved_pme_reciprocal = moved_pme_dict['reciprocal']  # Only take reciprocal space part
    print(f"Moved PME energy result: {moved_pme_result}")
    print(f"Moved PME reciprocal energy: {moved_pme_reciprocal}")
    sys.stdout.flush()
    
    # Step 5: Calculate energy using PGP after movement
    print("Calculating PGP after movement energy...")
    sys.stdout.flush()
    moved_pgp_energy = pygcmc.calculateMoleculeEnergy(system)
    
    # Print moved energy
    print(f"Moved PME reciprocal energy: {moved_pme_reciprocal}")
    print(f"Moved PGP energy: {moved_pgp_energy}")
    sys.stdout.flush()
    
    # Calculate energy change - only use reciprocal space part
    pme_energy_change = moved_pme_reciprocal - initial_pme_reciprocal
    pgp_energy_change = moved_pgp_energy - initial_pgp_energy
    
    print(f"PME reciprocal energy change: {pme_energy_change}")
    print(f"PGP energy change: {pgp_energy_change}")
    sys.stdout.flush()
    
    if abs(pme_energy_change) < 1e-10:
        print("PME energy change is too small, cannot compute relative error")
        assert abs(pgp_energy_change) < 1e-10, f"PGP energy should also be close to zero"
    else:
        # Calculate relative error, allowing some error range (e.g., 10%)
        relative_error = abs((pgp_energy_change - pme_energy_change) / pme_energy_change)
        print(f"Relative error: {relative_error * 100:.4f}%")
        sys.stdout.flush()
        
        # Verify PGP and PME calculated energy changes are consistent within error range
        assert relative_error < 0.1, f"Relative error too large: {relative_error*100:.2f}%"  # Allow 10% error

def test_compare_ewald_pme_pgp_complex():
    """
    Use a more complex system to compare Ewald, PME, and PGP algorithms
    
    This test:
    1. Creates a complex system with fixed part containing multiple charged particles
    2. Moving part contains 2-3 atoms, distance exceeds cutoff
    3. Compare three methods calculated reciprocal space energy change
    """
    # Set system parameters
    box_size = 8.0  # nm - use a larger box
    cutoff = 1.0    # nm
    potential_cutoff = 1.0  # nm
    box = [box_size, box_size, box_size]
    
    # Ewald and PME parameters
    alpha = 0.29    # 1/nm
    mesh_size = [32, 32, 32]
    potential_grid_size = [32, 32, 32]
    spline_order = 4
    tolerance = 1e-5
    
    print("Creating complex test system...")
    sys.stdout.flush()
    
    # Create test system
    system = MCState()
    system.info.box = box
    system.info.setTemperature(300.0)
    system.info.cutoff = cutoff
    
    # Set force field
    ff = MCForceField()
    ff.numTotalTypes = 2 
    ff.ljSigma = [0.333, 0.3875, 0.3875, 0.442]
    ff.ljEps = [0.0115, 0.0693, 0.0693, 0.4184]
    system.forcefield = ff
    
    atoms = []
    residues = []
    
    # Create fixed part - 8 ions forming a cube
    fixed_positions = [
        (1.0, 1.0, 1.0),
        (1.0, 1.0, box_size-1.0),
        (1.0, box_size-1.0, 1.0),
        (1.0, box_size-1.0, box_size-1.0),
        (box_size-1.0, 1.0, 1.0),
        (box_size-1.0, 1.0, box_size-1.0),
        (box_size-1.0, box_size-1.0, 1.0),
        (box_size-1.0, box_size-1.0, box_size-1.0)
    ]
    
    # Add fixed ions
    for i, pos in enumerate(fixed_positions):
        ion = MCAtom()
        ion.x, ion.y, ion.z = pos
        ion.charge = 1.0 if i % 2 == 0 else -1.0  # Alternating positive/negative charges
        ion.type = 0 if i % 2 == 0 else 1
        atoms.append(ion)
    
    # Create fixed residue
    fixed_res = MCResidue()
    fixed_res.atomStart = 0
    fixed_res.atomCount = len(fixed_positions)
    fixed_res.active = True
    fixed_res.fixed = True
    residues.append(fixed_res)
    
    # Create moving part - 3 atoms forming a water molecule
    # Position in the center of the box, ensuring distance greater than cutoff from fixed parts
    mobile_atoms = [
        (box_size/2.0, box_size/2.0, box_size/2.0),       # Central oxygen atom
        (box_size/2.0 + 0.1, box_size/2.0, box_size/2.0), # Hydrogen atom 1
        (box_size/2.0, box_size/2.0 + 0.1, box_size/2.0)  # Hydrogen atom 2
    ]
    
    mobile_charges = [-0.8, 0.4, 0.4]  # Water molecule charges
    
    # Add moving atoms
    for pos, q in zip(mobile_atoms, mobile_charges):
        atom = MCAtom()
        atom.x, atom.y, atom.z = pos
        atom.charge = q
        atom.type = 0
        atoms.append(atom)
    
    # Create moving residue
    mobile_res = MCResidue()
    mobile_res.atomStart = len(fixed_positions)
    mobile_res.atomCount = len(mobile_atoms)
    mobile_res.active = True
    mobile_res.fixed = False
    residues.append(mobile_res)
    
    # Set up system
    system.atoms = atoms
    system.residues = residues
    system.activeAtomCount = len(atoms)
    system.activeResidueCount = len(residues)
    
    # Set moving residue info
    system.movementResidues.clear()
    movement_info = pygcmc.MCMovementResidueInfo()
    movement_info.startIndex = 1  # Index of second residue (moving residue)
    movement_info.activeCount = 1  # Only one moving residue
    system.movementResidues.append(movement_info)
    
    print(f"System created: {system.activeAtomCount} atoms, {system.activeResidueCount} residues")
    print(f"Fixed atoms: {len(fixed_positions)}, Moving atoms: {len(mobile_atoms)}")
    sys.stdout.flush()
    
    # Initialize different charging methods
    print("Setting calculation parameters...")
    
    # PME parameters
    pygcmc.setPMEParameters(alpha, mesh_size, spline_order, tolerance)
    pygcmc.initializePMEParameters(cutoff, box, alpha)
    
    # PGP parameters
    pygcmc.setPGPParameters(alpha, mesh_size, potential_cutoff, 
                           potential_grid_size, spline_order, tolerance)
    
    # Ewald parameters (using computeEwaldEnergy function)
    # Note: This assumes your library has a function to calculate standard Ewald energy
    
    # Step 1: Calculate initial energy
    print("Calculating initial energy...")
    
    # PME energy
    initial_pme_result = pygcmc.computeMovementEnergyPME(system)
    initial_pme_energy = initial_pme_result[0]
    initial_pme_dict = initial_pme_result[2]
    initial_pme_reciprocal = initial_pme_dict['reciprocal']
    
    # Ewald energy (if available)
    # initial_ewald_energy = pygcmc.computeEwaldEnergy(system)
    
    # PGP energy
    pygcmc.precomputeGridPotential(system, fixed_only=True)
    initial_pgp_energy = pygcmc.calculateMoleculeEnergy(system)
    
    print(f"Initial PME reciprocal energy: {initial_pme_reciprocal}")
    # print(f"Initial Ewald energy: {initial_ewald_energy}")
    print(f"Initial PGP energy: {initial_pgp_energy}")
    
    # Step 2: Move the residue (translate by 0.2 nm)
    translation = [0.2, 0.2, 0.2]
    print(f"Moving residue: translating {translation}...")
    
    # Move water molecule
    for i in range(mobile_res.atomCount):
        atom_index = mobile_res.atomStart + i
        system.atoms[atom_index].x += translation[0]
        system.atoms[atom_index].y += translation[1]
        system.atoms[atom_index].z += translation[2]
    
    # Step 3: Calculate energy after movement
    print("Calculating after movement energy...")
    
    # PME energy
    moved_pme_result = pygcmc.computeMovementEnergyPME(system)
    moved_pme_energy = moved_pme_result[0]
    moved_pme_dict = moved_pme_result[2]
    moved_pme_reciprocal = moved_pme_dict['reciprocal']
    
    # Ewald energy (if available)
    # moved_ewald_energy = pygcmc.computeEwaldEnergy(system)
    
    # PGP energy
    moved_pgp_energy = pygcmc.calculateMoleculeEnergy(system)
    
    print(f"Moved PME reciprocal energy: {moved_pme_reciprocal}")
    # print(f"Moved Ewald energy: {moved_ewald_energy}")
    print(f"Moved PGP energy: {moved_pgp_energy}")
    
    # Calculate energy change
    pme_energy_change = moved_pme_reciprocal - initial_pme_reciprocal
    # ewald_energy_change = moved_ewald_energy - initial_ewald_energy
    pgp_energy_change = moved_pgp_energy - initial_pgp_energy
    
    print(f"PME reciprocal energy change: {pme_energy_change}")
    # print(f"Ewald energy change: {ewald_energy_change}")
    print(f"PGP energy change: {pgp_energy_change}")
    
    # Calculate relative error
    if abs(pme_energy_change) > 1e-10:
        pgp_relative_error = abs((pgp_energy_change - pme_energy_change) / pme_energy_change)
        print(f"PGP relative error: {pgp_relative_error*100:.4f}%")
        assert pgp_relative_error < 0.1, f"PGP relative error too large: {pgp_relative_error*100:.2f}%"
    
    # If Ewald calculation is available, also compare with PME
    # if abs(pme_energy_change) > 1e-10:
    #     ewald_relative_error = abs((ewald_energy_change - pme_energy_change) / pme_energy_change)
    #     print(f"Ewald relative error: {ewald_relative_error*100:.4f}%")
    #     assert ewald_relative_error < 0.1, f"Ewald relative error too large: {ewald_relative_error*100:.2f}%"

def test_compare_ewald_pme_pgp_planar():
    """
    Use planar system to compare Ewald, PME and PGP algorithms
    
    Features:
    1. All atoms are on a plane at z=4.0 nm
    2. Moving molecule is at a distance greater than cutoff from fixed parts
    3. The plane is parallel to grid planes
    4. Atom positions deliberately deviate from grid points
    5. Using a coarse 8x8x8 grid for better observation
    """
    # Set system parameters
    box_size = 8.0  # nm
    cutoff = 1.0    # nm
    potential_cutoff = 1.0  # nm
    box = [box_size, box_size, box_size]
    
    # Calculation parameters
    alpha = 0.29    # 1/nm
    mesh_size = [8, 8, 8]  # Changed to 8x8x8 coarse grid
    potential_grid_size = [8, 8, 8]  # Also changed to 8x8x8
    spline_order = 4
    tolerance = 1e-5
    
    # Calculate grid spacing
    grid_spacing = box_size / mesh_size[0]
    print(f"Grid spacing: {grid_spacing:.3f} nm")
    
    # Choose plane z-coordinate (ensure parallel to grid planes)
    z_plane = 4.0  # nm
    
    print("Creating planar test system...")
    sys.stdout.flush()
    
    # Create test system
    system = MCState()
    system.info.box = box
    system.info.setTemperature(300.0)
    system.info.cutoff = cutoff
    
    # Set force field
    ff = MCForceField()
    ff.numTotalTypes = 2 
    ff.ljSigma = [0.333, 0.3875, 0.3875, 0.442]
    ff.ljEps = [0.0115, 0.0693, 0.0693, 0.4184]
    system.forcefield = ff
    
    atoms = []
    residues = []
    
    # Create fixed part - 4 ions forming a square
    # Positions deliberately deviate from grid points
    fixed_positions = [
        (1.0 + grid_spacing/3, 1.0 + grid_spacing/3, z_plane),  # Bottom left
        (1.0 + grid_spacing/3, box_size-1.0 + grid_spacing/3, z_plane),  # Top left
        (box_size-1.0 + grid_spacing/3, 1.0 + grid_spacing/3, z_plane),  # Bottom right
        (box_size-1.0 + grid_spacing/3, box_size-1.0 + grid_spacing/3, z_plane)  # Top right
    ]
    
    # Add fixed ions
    for i, pos in enumerate(fixed_positions):
        ion = MCAtom()
        ion.x, ion.y, ion.z = pos
        ion.charge = 1.0 if i % 2 == 0 else -1.0  # Alternating positive/negative charges
        ion.type = 0 if i % 2 == 0 else 1
        atoms.append(ion)
        
        # Print atom position and nearest grid point
        grid_x = round(pos[0] / grid_spacing)
        grid_y = round(pos[1] / grid_spacing)
        grid_z = round(pos[2] / grid_spacing)
        print(f"Fixed ion {i}: position=({pos[0]:.3f}, {pos[1]:.3f}, {pos[2]:.3f})")
        print(f"   Nearest grid point: ({grid_x}, {grid_y}, {grid_z})")
        print(f"   Deviation from grid point: ({pos[0]-grid_x*grid_spacing:.3f}, {pos[1]-grid_y*grid_spacing:.3f}, {pos[2]-grid_z*grid_spacing:.3f})")
    
    # Create fixed residue
    fixed_res = MCResidue()
    fixed_res.atomStart = 0
    fixed_res.atomCount = len(fixed_positions)
    fixed_res.active = True
    fixed_res.fixed = True
    residues.append(fixed_res)
    
    # Create moving part - 3 atoms forming a water molecule
    # Position in the center of the box, ensuring distance greater than cutoff from fixed part
    mobile_atoms = [
        (box_size/2.0 + grid_spacing/3, box_size/2.0 + grid_spacing/3, z_plane),  # Oxygen atom
        (box_size/2.0 + grid_spacing/3 + 0.1, box_size/2.0 + grid_spacing/3, z_plane),  # Hydrogen atom 1
        (box_size/2.0 + grid_spacing/3, box_size/2.0 + grid_spacing/3 + 0.1, z_plane)  # Hydrogen atom 2
    ]
    
    mobile_charges = [-0.8, 0.4, 0.4]  # Water molecule charges
    
    # Add moving atoms
    for pos, q in zip(mobile_atoms, mobile_charges):
        atom = MCAtom()
        atom.x, atom.y, atom.z = pos
        atom.charge = q
        atom.type = 0
        atoms.append(atom)
    
    # Create moving residue
    mobile_res = MCResidue()
    mobile_res.atomStart = len(fixed_positions)
    mobile_res.atomCount = len(mobile_atoms)
    mobile_res.active = True
    mobile_res.fixed = False
    residues.append(mobile_res)
    
    # Set up system
    system.atoms = atoms
    system.residues = residues
    system.activeAtomCount = len(atoms)
    system.activeResidueCount = len(residues)
    
    # Set moving residue info
    system.movementResidues.clear()
    movement_info = pygcmc.MCMovementResidueInfo()
    movement_info.startIndex = 1  # Index of second residue (moving residue)
    movement_info.activeCount = 1  # Only one moving residue
    system.movementResidues.append(movement_info)
    
    print(f"System created: {system.activeAtomCount} atoms, {system.activeResidueCount} residues")
    print(f"Fixed atoms: {len(fixed_positions)}, Moving atoms: {len(mobile_atoms)}")
    sys.stdout.flush()
    
    # Initialize different charging methods
    print("Setting calculation parameters...")
    
    # PME parameters
    pygcmc.setPMEParameters(alpha, mesh_size, spline_order, tolerance)
    pygcmc.initializePMEParameters(cutoff, box, alpha)
    
    # PGP parameters
    pygcmc.setPGPParameters(alpha, mesh_size, potential_cutoff, 
                           potential_grid_size, spline_order, tolerance)
    
    # Ewald parameters (using computeEwaldEnergy function)
    # Note: This assumes your library has a function to calculate standard Ewald energy
    
    # Step 1: Calculate initial energy
    print("Calculating initial energy...")
    
    # PME energy
    initial_pme_result = pygcmc.computeMovementEnergyPME(system)
    initial_pme_energy = initial_pme_result[0]
    initial_pme_dict = initial_pme_result[2]
    initial_pme_reciprocal = initial_pme_dict['reciprocal']
    
    # Ewald energy (if available)
    # initial_ewald_energy = pygcmc.computeEwaldEnergy(system)
    
    # PGP energy
    pygcmc.precomputeGridPotential(system, fixed_only=True)
    initial_pgp_energy = pygcmc.calculateMoleculeEnergy(system)
    
    print(f"Initial PME reciprocal energy: {initial_pme_reciprocal}")
    # print(f"Initial Ewald energy: {initial_ewald_energy}")
    print(f"Initial PGP energy: {initial_pgp_energy}")
    
    # Step 2: Move the residue (translate by 0.2 nm)
    translation = [0.2, 0.2, 0.2]
    print(f"Moving residue: translating {translation}...")
    
    # Move water molecule
    for i in range(mobile_res.atomCount):
        atom_index = mobile_res.atomStart + i
        system.atoms[atom_index].x += translation[0]
        system.atoms[atom_index].y += translation[1]
        system.atoms[atom_index].z += translation[2]
    
    # Step 3: Calculate energy after movement
    print("Calculating after movement energy...")
    
    # PME energy
    moved_pme_result = pygcmc.computeMovementEnergyPME(system)
    moved_pme_energy = moved_pme_result[0]
    moved_pme_dict = moved_pme_result[2]
    moved_pme_reciprocal = moved_pme_dict['reciprocal']
    
    # Ewald energy (if available)
    # moved_ewald_energy = pygcmc.computeEwaldEnergy(system)
    
    # PGP energy
    moved_pgp_energy = pygcmc.calculateMoleculeEnergy(system)
    
    print(f"Moved PME reciprocal energy: {moved_pme_reciprocal}")
    # print(f"Moved Ewald energy: {moved_ewald_energy}")
    print(f"Moved PGP energy: {moved_pgp_energy}")
    
    # Calculate energy change
    pme_energy_change = moved_pme_reciprocal - initial_pme_reciprocal
    # ewald_energy_change = moved_ewald_energy - initial_ewald_energy
    pgp_energy_change = moved_pgp_energy - initial_pgp_energy
    
    print(f"PME reciprocal energy change: {pme_energy_change}")
    # print(f"Ewald energy change: {ewald_energy_change}")
    print(f"PGP energy change: {pgp_energy_change}")
    
    # Calculate relative error
    if abs(pme_energy_change) > 1e-10:
        pgp_relative_error = abs((pgp_energy_change - pme_energy_change) / pme_energy_change)
        print(f"PGP relative error: {pgp_relative_error*100:.4f}%")
        assert pgp_relative_error < 0.1, f"PGP relative error too large: {pgp_relative_error*100:.2f}%"
    
    # If Ewald calculation is available, also compare error with PME
    # if abs(pme_energy_change) > 1e-10:
    #     ewald_relative_error = abs((ewald_energy_change - pme_energy_change) / pme_energy_change)
    #     print(f"Ewald relative error: {ewald_relative_error*100:.4f}%")
    #     assert ewald_relative_error < 0.1, f"Ewald relative error too large: {ewald_relative_error*100:.2f}%"

def test_compare_ewald_pme_pgp_asymmetric():
    """
    Test whether PGP calculation is accurate in systems with asymmetric charge distribution
    
    Features:
    1. Fixed part contains multiple charged particles with asymmetric distribution
    2. Moving residue is a water molecule, far from fixed part
    3. Execute multiple random movements to ensure distance from fixed part is always > cutoff
    4. Verify accuracy by comparing energy calculated by Ewald, PME, and PGP
    """
    # Set system parameters
    box_size = 8.0  # nm - use a larger box
    cutoff = 1.0    # nm
    potential_cutoff = 1.0  # nm
    box = [box_size, box_size, box_size]
    
    # Ewald and PME parameters
    alpha = 0.29    # 1/nm
    kmax = [8, 8, 8]  # Number of k-space vectors for Ewald calculation
    mesh_size = [32, 32, 32]
    potential_grid_size = [32, 32, 32]
    spline_order = 4
    tolerance = 1e-5
    
    # Define function to calculate distance in periodic boundary conditions
    def calculate_pbc_distance(pos1, pos2, box_size):
        """Calculate distance between two points in periodic boundary conditions"""
        dx = abs(pos1[0] - pos2[0])
        dy = abs(pos1[1] - pos2[1])
        dz = abs(pos1[2] - pos2[2])

        # Apply periodic boundary conditions
        if dx > box_size/2:
            dx = box_size - dx
        if dy > box_size/2:
            dy = box_size - dy
        if dz > box_size/2:
            dz = box_size - dz

        return math.sqrt(dx*dx + dy*dy + dz*dz)
    
    # Check if movement is safe (all distances are > cutoff)
    def is_safe_position(mobile_positions, fixed_positions, cutoff, box_size):
        """Check if all distances between mobile particles and fixed particles are > cutoff"""
        for mobile_pos in mobile_positions:
            for fixed_pos, _ in fixed_positions:
                distance = calculate_pbc_distance(mobile_pos, fixed_pos, box_size)
                if distance <= cutoff:
                    return False, distance
        return True, None
    
    # Generate safe random movement vector
    def generate_safe_move(current_positions, fixed_positions, cutoff, box_size, max_step=0.3):
        """Generate a safe random movement vector, ensuring all particles are at distance > cutoff from fixed particles after movement"""
        for attempt in range(100):  # Try up to 100 times
            # Generate random displacement
            dx = (random.random() - 0.5) * 2 * max_step
            dy = (random.random() - 0.5) * 2 * max_step
            dz = (random.random() - 0.5) * 2 * max_step
            
            # Calculate new positions
            new_positions = []
            for pos in current_positions:
                new_pos = (
                    (pos[0] + dx) % box_size,
                    (pos[1] + dy) % box_size,
                    (pos[2] + dz) % box_size
                )
                new_positions.append(new_pos)
            
            # Check if new positions are safe
            is_safe, min_distance = is_safe_position(new_positions, fixed_positions, cutoff, box_size)
            if is_safe:
                return (dx, dy, dz), new_positions
        
        # If 100 attempts fail, return a smaller movement
        print("Warning: 100 attempts failed to find a safe position, using smaller movement")
        dx = 0.05
        dy = 0.05
        dz = 0.05
        new_positions = []
        for pos in current_positions:
            new_pos = (
                (pos[0] + dx) % box_size,
                (pos[1] + dy) % box_size,
                (pos[2] + dz) % box_size
            )
            new_positions.append(new_pos)
        return (dx, dy, dz), new_positions
    
    print("Creating asymmetric charge distribution test system...")
    sys.stdout.flush()
    
    # Create test system
    system = MCState()
    system.info.box = box
    system.info.setTemperature(300.0)
    system.info.cutoff = cutoff
    
    # Set up force field
    ff = MCForceField()
    ff.numTotalTypes = 2 
    ff.ljSigma = [0.333, 0.3875, 0.3875, 0.442]
    ff.ljEps = [0.0115, 0.0693, 0.0693, 0.4184]
    system.forcefield = ff
    
    atoms = []
    residues = []
    
    # Create fixed part - asymmetric charged particle distribution
    fixed_particles = [
        # Position (x, y, z)                 Charge
        ((1.0, 1.0, 1.0),               1.0),  # Positive charge
        ((1.5, 1.0, 1.2),              -0.8),  # Negative charge
        ((1.3, 1.7, 1.5),               0.6),  # Positive charge
        ((0.8, 1.6, 0.9),              -0.7),  # Negative charge
        ((box_size-1.0, 1.0, 1.0),      1.0),  # Positive charge
        ((box_size-1.5, 1.2, 1.3),     -0.5),  # Negative charge
        ((1.0, box_size-1.0, 1.0),      0.8),  # Positive charge
        ((1.2, box_size-1.5, 0.9),     -0.6),  # Negative charge
        ((1.0, 1.0, box_size-1.0),      0.9),  # Positive charge
        ((1.4, 1.3, box_size-1.4),     -0.7),  # Negative charge
        # Add more particles to increase system complexity
        ((box_size-2.0, box_size-2.0, 2.0),  0.7),  # Positive charge
        ((box_size-2.5, box_size-2.3, 2.2), -0.4),  # Negative charge
        ((2.0, box_size-2.0, box_size-2.0),  0.5),  # Positive charge
        ((2.2, box_size-2.2, box_size-2.4), -0.3),  # Negative charge
        ((box_size-2.0, 2.0, box_size-2.0), -1.5),  # Negative charge, roughly balancing system
    ]
    
    # Confirm particle count
    print(f"Number of fixed particles: {len(fixed_particles)}")
    
    # Calculate total charge of fixed part
    total_fixed_charge = sum(charge for _, charge in fixed_particles)
    print(f"Total fixed charge: {total_fixed_charge}")
    
    # Add fixed particles
    for i, (pos, charge) in enumerate(fixed_particles):
        atom = MCAtom()
        atom.x, atom.y, atom.z = pos
        atom.charge = charge
        atom.type = 0 if charge > 0 else 1  # Positive charge use type 0, negative charge use type 1
        atoms.append(atom)
    
    # Create fixed residue
    fixed_res = MCResidue()
    fixed_res.atomStart = 0
    fixed_res.atomCount = len(fixed_particles)
    fixed_res.active = True
    fixed_res.fixed = True
    residues.append(fixed_res)
    
    # Calculate box center position to ensure moving residue is far from fixed particles
    center_x = box_size / 2.0
    center_y = box_size / 2.0
    center_z = box_size / 2.0
    
    # Create moving part - a water molecule placed in the center of the box
    # Water molecule charge parameters: oxygen atom -0.8, two hydrogen atoms each +0.4
    # Water molecule bond length: O-H about 0.1 nm
    
    # Set water molecule position
    oxygen_pos = (center_x, center_y, center_z)
    hydrogen1_pos = (center_x + 0.1, center_y, center_z)  # First H atom
    hydrogen2_pos = (center_x, center_y + 0.1, center_z)  # Second H atom
    
    # Water molecule charge
    oxygen_charge = -0.8
    hydrogen_charge = 0.4  # Each hydrogen atom
    
    # Record starting index of moving part
    mobile_start_idx = len(atoms)
    
    # Create oxygen atom
    o_atom = MCAtom()
    o_atom.x, o_atom.y, o_atom.z = oxygen_pos
    o_atom.charge = oxygen_charge
    o_atom.type = 1  # Oxygen atom type
    atoms.append(o_atom)
    
    # Create first hydrogen atom
    h1_atom = MCAtom()
    h1_atom.x, h1_atom.y, h1_atom.z = hydrogen1_pos
    h1_atom.charge = hydrogen_charge
    h1_atom.type = 0  # Hydrogen atom type
    atoms.append(h1_atom)
    
    # Create second hydrogen atom
    h2_atom = MCAtom()
    h2_atom.x, h2_atom.y, h2_atom.z = hydrogen2_pos
    h2_atom.charge = hydrogen_charge
    h2_atom.type = 0  # Hydrogen atom type
    atoms.append(h2_atom)
    
    print(f"Moving water molecule: O position=({oxygen_pos[0]}, {oxygen_pos[1]}, {oxygen_pos[2]}), charge={oxygen_charge}")
    print(f"  H1 position=({hydrogen1_pos[0]}, {hydrogen1_pos[1]}, {hydrogen1_pos[2]}), charge={hydrogen_charge}")
    print(f"  H2 position=({hydrogen2_pos[0]}, {hydrogen2_pos[1]}, {hydrogen2_pos[2]}), charge={hydrogen_charge}")
    print(f"   Water molecule total charge: {oxygen_charge + 2*hydrogen_charge}")
    
    # Create moving residue (water molecule)
    mobile_res = MCResidue()
    mobile_res.atomStart = mobile_start_idx
    mobile_res.atomCount = 3  # Water molecule has 3 atoms
    mobile_res.active = True
    mobile_res.fixed = False
    residues.append(mobile_res)
    
    # Set system
    system.atoms = atoms
    system.residues = residues
    system.activeAtomCount = len(atoms)
    system.activeResidueCount = len(residues)
    
    # Calculate system total charge
    total_system_charge = sum(atom.charge for atom in atoms)
    print(f"System total charge: {total_system_charge}")
    # Allow system to have small charge, no strict assertion needed
    
    # Set moving residue info
    system.movementResidues.clear()
    movement_info = pygcmc.MCMovementResidueInfo()
    movement_info.startIndex = 1  # Second residue (moving residue) index
    movement_info.activeCount = 1  # Only one moving residue
    system.movementResidues.append(movement_info)
    
    print(f"System created: {system.activeAtomCount} atoms, {system.activeResidueCount} residues")
    print(f"Fixed atoms: {len(fixed_particles)}, Moving atoms: 3 (Water molecule)")
    sys.stdout.flush()
    
    # Initialize Ewald, PME, and PGP
    print("Setting calculation parameters...")
    
    # Ewald parameters initialization
    print("Initializing Ewald parameters...")
    pygcmc.setEwaldParameters(alpha, kmax)
    pygcmc.initializeEwaldParameters(cutoff, box, alpha)
    
    # PME parameters initialization
    print("Initializing PME parameters...")
    pygcmc.setPMEParameters(alpha, mesh_size, spline_order, tolerance)
    pygcmc.initializePMEParameters(cutoff, box, alpha)
    
    # PGP parameters initialization
    print("Initializing PGP parameters...")
    pygcmc.setPGPParameters(alpha, mesh_size, potential_cutoff, 
                         potential_grid_size, spline_order, tolerance)
    
    # Precompute grid potential for fixed parts
    print("Precomputing fixed parts grid potential...")
    pygcmc.precomputeGridPotential(system, fixed_only=True)
    
    # Set number of random movements to execute
    num_moves = 5  # Reduce test times to speed up test
    print(f"Executing {num_moves} random movement tests")
    
    # Store relative errors between methods
    pgp_pme_errors = []
    ewald_pme_errors = []
    pgp_ewald_errors = []
    
    # Verify initial position is safe
    current_positions = []
    for i in range(mobile_res.atomCount):
        atom_idx = mobile_res.atomStart + i
        current_positions.append((
            system.atoms[atom_idx].x,
            system.atoms[atom_idx].y,
            system.atoms[atom_idx].z
        ))
    
    is_safe, min_dist = is_safe_position(current_positions, fixed_particles, cutoff, box_size)
    if not is_safe:
        print(f"Warning: Initial position not safe, minimum distance is {min_dist} nm")
    else:
        print(f"Initial position safe, minimum distance from fixed particles > {cutoff} nm")
    
    # Execute multiple random movements
    for move_idx in range(num_moves):
        print(f"\nExecuting {move_idx+1}/{num_moves} random movement test")
        
        # Step 1: Calculate initial energy
        print("Calculating initial energy...")
        
        # Ewald energy calculation
        initial_ewald_result = pygcmc.computeSystemEnergyEwald(system)
        initial_ewald_elec = initial_ewald_result[0]
        initial_ewald_dict = initial_ewald_result[2]
        initial_ewald_reciprocal = initial_ewald_dict['reciprocal']
        
        # PME energy calculation
        initial_pme_result = pygcmc.computeMovementEnergyPME(system)
        initial_pme_energy = initial_pme_result[0]
        initial_pme_dict = initial_pme_result[2]
        initial_pme_reciprocal = initial_pme_dict['reciprocal']
        
        # PGP energy calculation
        initial_pgp_energy = pygcmc.calculateMoleculeEnergy(system)
        
        print(f"Initial Ewald reciprocal energy: {initial_ewald_reciprocal}")
        print(f"Initial PME reciprocal energy: {initial_pme_reciprocal}")
        print(f"Initial PGP energy: {initial_pgp_energy}")
        
        # Step 2: Generate safe random movement and apply
        current_positions = []
        for i in range(mobile_res.atomCount):
            atom_idx = mobile_res.atomStart + i
            current_positions.append((
                system.atoms[atom_idx].x,
                system.atoms[atom_idx].y,
                system.atoms[atom_idx].z
            ))
        
        delta, new_positions = generate_safe_move(current_positions, fixed_particles, cutoff, box_size)
        print(f"Random movement vector: {delta}")
        
        # Move all particles in moving residue
        for i in range(mobile_res.atomCount):
            atom_idx = mobile_res.atomStart + i
            system.atoms[atom_idx].x = new_positions[i][0]
            system.atoms[atom_idx].y = new_positions[i][1]
            system.atoms[atom_idx].z = new_positions[i][2]
            print(f"Moved atom {i} to: ({new_positions[i][0]:.4f}, {new_positions[i][1]:.4f}, {new_positions[i][2]:.4f})")
        
        # Verify moved position is safe
        is_safe, min_dist = is_safe_position(new_positions, fixed_particles, cutoff, box_size)
        if not is_safe:
            print(f"Warning: Moved position not safe, minimum distance is {min_dist} nm")
            assert min_dist > cutoff, f"Moved distance ({min_dist} nm) less than cutoff ({cutoff} nm), will introduce real space energy"
        else:
            min_dist = float('inf')
            for mobile_pos in new_positions:
                for fixed_pos, _ in fixed_particles:
                    dist = calculate_pbc_distance(mobile_pos, fixed_pos, box_size)
                    min_dist = min(min_dist, dist)
            print(f"Moved position safe, minimum distance from fixed particles: {min_dist:.4f} nm (cutoff={cutoff} nm)")
        
        # Step 3: Calculate moved energy
        print("Calculating moved energy...")
        
        # Ewald energy calculation
        moved_ewald_result = pygcmc.computeSystemEnergyEwald(system)
        moved_ewald_elec = moved_ewald_result[0]
        moved_ewald_dict = moved_ewald_result[2]
        moved_ewald_reciprocal = moved_ewald_dict['reciprocal']
        
        # PME energy calculation
        moved_pme_result = pygcmc.computeMovementEnergyPME(system)
        moved_pme_energy = moved_pme_result[0]
        moved_pme_dict = moved_pme_result[2]
        moved_pme_reciprocal = moved_pme_dict['reciprocal']
        
        # Check direct space energy is 0
        moved_pme_direct = moved_pme_dict.get('direct', 0.0)
        if abs(moved_pme_direct) > 1e-10:
            print(f"Warning: PME direct space energy not zero: {moved_pme_direct}")
            print("This means there are particles < cutoff!")
        
        # PGP energy calculation
        moved_pgp_energy = pygcmc.calculateMoleculeEnergy(system)
        
        print(f"Moved Ewald reciprocal energy: {moved_ewald_reciprocal}")
        print(f"Moved PME reciprocal energy: {moved_pme_reciprocal}")
        print(f"Moved PGP energy: {moved_pgp_energy}")
        
        # Calculate energy change
        ewald_energy_change = moved_ewald_reciprocal - initial_ewald_reciprocal
        pme_energy_change = moved_pme_reciprocal - initial_pme_reciprocal
        pgp_energy_change = moved_pgp_energy - initial_pgp_energy
        
        print(f"Ewald reciprocal energy change: {ewald_energy_change}")
        print(f"PME reciprocal energy change: {pme_energy_change}")
        print(f"PGP energy change: {pgp_energy_change}")
        
        # Calculate relative error
        if abs(pme_energy_change) > 1e-6:
            # PGP relative error compared to PME
            pgp_pme_error = abs((pgp_energy_change - pme_energy_change) / pme_energy_change)
            print(f"PGP relative error: {pgp_pme_error*100:.4f}%")
            pgp_pme_errors.append(pgp_pme_error)
            
            # Ewald relative error compared to PME
            ewald_pme_error = abs((ewald_energy_change - pme_energy_change) / pme_energy_change)
            print(f"Ewald relative error: {ewald_pme_error*100:.4f}%")
            ewald_pme_errors.append(ewald_pme_error)
            
            # PGP relative error compared to Ewald
            if abs(ewald_energy_change) > 1e-6:
                pgp_ewald_error = abs((pgp_energy_change - ewald_energy_change) / ewald_energy_change)
                print(f"PGP relative error: {pgp_ewald_error*100:.4f}%")
                pgp_ewald_errors.append(pgp_ewald_error)
        else:
            print("PME energy change near zero, skipping relative error calculation")
    
    # Calculate average error
    print("\nError analysis statistics:")
    
    if pgp_pme_errors:
        avg_pgp_pme_error = sum(pgp_pme_errors) / len(pgp_pme_errors)
        print(f"Average PGP relative error: {avg_pgp_pme_error*100:.4f}%")
        
        # Use more lenient error tolerance because PGP is an approximation method
        acceptable_error = 0.5  # Allow 50% error
        assert avg_pgp_pme_error < acceptable_error, f"Average PGP relative error too large: {avg_pgp_pme_error*100:.2f}%"
    else:
        print("PGP relative error: No valid error data for statistics")
    
    if ewald_pme_errors:
        avg_ewald_pme_error = sum(ewald_pme_errors) / len(ewald_pme_errors)
        print(f"Average Ewald relative error: {avg_ewald_pme_error*100:.4f}%")
        
        # Ewald and PME should theoretically have very high consistency
        acceptable_error = 0.2  # Allow 20% error
        assert avg_ewald_pme_error < acceptable_error, f"Average Ewald relative error too large: {avg_ewald_pme_error*100:.2f}%"
    else:
        print("Ewald relative error: No valid error data for statistics")
    
    if pgp_ewald_errors:
        avg_pgp_ewald_error = sum(pgp_ewald_errors) / len(pgp_ewald_errors)
        print(f"Average PGP relative error: {avg_pgp_ewald_error*100:.4f}%")
        
        # Use more lenient error tolerance because PGP is an approximation method
        acceptable_error = 0.5  # Allow 50% error
        assert avg_pgp_ewald_error < acceptable_error, f"Average PGP relative error too large: {avg_pgp_ewald_error*100:.2f}%"
    else:
        print("PGP relative error: No valid error data for statistics")

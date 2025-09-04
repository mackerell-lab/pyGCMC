# tests/simulation/pgp/planar_systems.py
"""PGP planar systems tests."""

import pytest
import pygcmc
from pygcmc import MCState, MCAtom, MCResidue, MCForceField, MCMovementResidueInfo
import sys


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
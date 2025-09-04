# tests/simulation/energyPME/crystal_file_helpers.py
"""Helper functions for reading crystal data files."""

import math
import os
import re
from pygcmc import MCState, MCAtom, MCResidue, MCForceField


def create_nacl_crystal_from_file(data_file_paths=None, box_size=2.82, cutoff=1.0, num_particles=1000):
    """
    Create a NaCl crystal system from nacl_crystal.dat file.
    
    This function reads atomic positions from the nacl_crystal.dat file
    (same as used in pme.cpp testEwaldExact function) and creates a MCState
    with the same configuration.
    
    Args:
        data_file_paths: List of possible paths to find nacl_crystal.dat, default None will use standard paths
        box_size: Box size in nm (default 2.82, same as in pme.cpp)
        cutoff: Cutoff distance in nm (default 1.0, same as in pme.cpp)
        num_particles: Expected number of particles (default 1000, same as in pme.cpp)
    
    Returns:
        MCState object with the NaCl crystal configuration
    """
    # Set up the system
    state = MCState()
    state.info.box = [box_size, box_size, box_size]
    state.info.setTemperature(300.0)
    state.info.cutoff = cutoff
    
    # Set force field parameters - same as in pme.cpp
    ff = MCForceField()
    ff.numTotalTypes = 2  # Na+ and Cl-
    
    # LJ parameters
    sigma_na = 0.333  # nm
    sigma_cl = 0.442  # nm
    eps_na = 0.0115   # kJ/mol
    eps_cl = 0.4184   # kJ/mol
    
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
    
    # If no data file paths provided, use default paths
    if data_file_paths is None:
        data_file_paths = [
            '../pygcmc_dev/tests/data/nacl_crystal.dat',  # Relative to build directory
            '../tests/data/nacl_crystal.dat',             # Relative to current directory
            'tests/data/nacl_crystal.dat',                # From project root
            '/home/zhaomt/gcmc/test100/pygcmc_dev/tests/data/nacl_crystal.dat'  # Absolute path
        ]
    
    # Find the data file
    data_file_path = None
    for path in data_file_paths:
        if os.path.exists(path):
            data_file_path = path
            break
    
    if data_file_path is None:
        raise FileNotFoundError("Could not find nacl_crystal.dat file in any of the expected locations")
    
    print(f"Found data file at: {data_file_path}")
    
    # Parse positions from nacl_crystal.dat
    positions = []
    with open(data_file_path, 'r') as f:
        for line in f:
            # Look for lines like: positions[0] = Vec3(0.141000,0.141000,0.141000);
            match = re.search(r'Vec3\(([^)]+)\)', line)
            if match:
                coords_str = match.group(1)
                x, y, z = map(float, coords_str.split(','))
                positions.append((x, y, z))
    
    print(f"Read {len(positions)} positions from file")
    
    if len(positions) != num_particles:
        print(f"Warning: Expected {num_particles} particles but found {len(positions)} in the file")
        # We'll still proceed with what we have
    
    # Create atoms with the exact positions from the file
    atoms = []
    residues = []
    
    # Assign the first half as Na+ and second half as Cl-
    half_count = len(positions) // 2
    
    for i, (x, y, z) in enumerate(positions):
        atom = MCAtom()
        atom.x = x
        atom.y = y
        atom.z = z
        
        if i < half_count:
            atom.charge = 1.0  # Na+
            atom.type = 0
        else:
            atom.charge = -1.0  # Cl-
            atom.type = 1
        
        atoms.append(atom)
        
        # Create residues (one per atom or one per pair)
        if i % 2 == 0:
            res = MCResidue()
            res.atomStart = i
            res.atomCount = 2 if i < len(positions) - 1 else 1  # Last atom might be alone
            res.active = True
            res.fixed = False
            residues.append(res)
    
    # Print charges summary
    total_charge = sum(1.0 if i < half_count else -1.0 for i in range(len(positions)))
    print(f"Total system charge: {total_charge}")
    print(f"Created {len(atoms)} atoms and {len(residues)} residues")
    
    state.atoms = atoms
    state.residues = residues
    state.activeAtomCount = len(atoms)
    state.activeResidueCount = len(residues)
    
    return state
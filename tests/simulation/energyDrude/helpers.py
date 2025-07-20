# tests/simulation/energyDrude/helpers.py
"""Helper functions for Drude energy tests."""

import math
import pygcmc
from pygcmc import MCState, MCInfo, MCAtom, MCResidue, MCForceField
from pygcmc import (
    addDrudeParticle,
    addDrudeScreenedPair,
    setDrudeSCFParameters,
    setDrudeSCFTolerance,
    computeSystemEnergyDrude,
    getNumDrudeParticles,
    initializeDrudeForce,
    clearDrudeForce,
    createDrudeSWM4Water,
    DrudeSCFParams,
    PSFParser
)

# Physical constants
ONE_4PI_EPS0 = 138.935456  # kJ/mol·nm·e^-2


def create_simple_drude_system():
    """
    Create a simple system with one atom and its Drude particle
    """
    state = MCState()
    
    # Set box size
    state.info.box = [10.0, 10.0, 10.0]
    state.info.setTemperature(300.0)
    state.info.cutoff = 5.0
    
    # Set force field parameters
    ff = MCForceField()
    ff.numTotalTypes = 2  # Parent atom and Drude
    ff.numMovementTypes = 2
    
    # LJ parameters (minimal, just for completeness)
    ff.ljSigma = [0.3, 0.1, 0.1, 0.1]
    ff.ljEps = [0.1, 0.0, 0.0, 0.0]
    
    state.forcefield = ff
    
    # Create parent atom
    parent = MCAtom()
    parent.x = 5.0
    parent.y = 5.0
    parent.z = 5.0
    parent.charge = 0.0  # Neutral parent (charge is on Drude)
    parent.type = 0
    # parent.mass = 12.0  # Mass not needed for SCF
    
    # Create Drude particle
    drude = MCAtom()
    drude.x = 5.1  # Slightly displaced
    drude.y = 5.0
    drude.z = 5.0
    drude.charge = -1.0  # Negative charge
    drude.type = 1
    # drude.mass = 0.4  # Mass not needed for SCF
    
    state.atoms = [parent, drude]
    state.activeAtomCount = 2
    
    # Create residues
    res = MCResidue()
    res.atomStart = 0
    res.atomCount = 2
    res.active = True
    res.fixed = False
    res.type = 0
    
    state.residues = [res]
    state.activeResidueCount = 1
    
    return state


def create_swm4_water(x=0.0, y=0.0, z=0.0):
    """
    Create a single SWM4-NDP water molecule
    
    Args:
        x, y, z: Position of oxygen atom
        
    Returns:
        atoms: List of 5 MCAtom objects (O, D, H1, H2, M)
        residue: MCResidue for the water molecule
    """
    atoms = []
    
    # Oxygen atom
    o = MCAtom()
    o.x = x
    o.y = y
    o.z = z
    o.charge = 1.71636
    o.type = 0
    # o.mass = 15.6  # Mass not needed for SCF
    atoms.append(o)
    
    # Drude particle on oxygen
    d = MCAtom()
    d.x = x
    d.y = y
    d.z = z + 0.01  # Small displacement
    d.charge = -1.71636
    d.type = 1
    # d.mass = 0.4  # Mass not needed for SCF
    atoms.append(d)
    
    # Hydrogen 1
    h1 = MCAtom()
    h1.x = x + 0.09572
    h1.y = y
    h1.z = z
    h1.charge = 0.55733
    h1.type = 2
    # h1.mass = 1.0  # Mass not needed for SCF
    atoms.append(h1)
    
    # Hydrogen 2
    h2 = MCAtom()
    h2.x = x - 0.023999
    h2.y = y + 0.092663
    h2.z = z
    h2.charge = 0.55733
    h2.type = 2
    # h2.mass = 1.0  # Mass not needed for SCF
    atoms.append(h2)
    
    # Virtual site M
    m = MCAtom()
    # Position is average of O, H1, H2 with weights
    m.x = x + 0.786646558 * 0 + 0.106676721 * 0.09572 + 0.106676721 * (-0.023999)
    m.y = y + 0.786646558 * 0 + 0.106676721 * 0 + 0.106676721 * 0.092663
    m.z = z
    m.charge = -1.11466
    m.type = 3
    # m.mass = 0.0  # Mass not needed for SCF
    atoms.append(m)
    
    # Create residue
    res = MCResidue()
    res.atomStart = 0  # Will be adjusted when adding to state
    res.atomCount = 5
    res.active = True
    res.fixed = False
    res.type = 0  # Water residue type
    
    return atoms, res


def create_water_box(n_waters, box_size=3.0):
    """
    Create a box of SWM4-NDP water molecules
    
    Args:
        n_waters: Number of water molecules
        box_size: Box size in nm
        
    Returns:
        state: MCState with water molecules
    """
    state = MCState()
    
    # Set box parameters
    state.info.box = [box_size, box_size, box_size]
    state.info.setTemperature(300.0)
    state.info.cutoff = 1.0
    
    # Set force field parameters
    ff = MCForceField()
    ff.numTotalTypes = 4  # O, D, H, M
    ff.numMovementTypes = 4
    
    # LJ parameters (simplified)
    # Only oxygen has LJ interactions in SWM4-NDP
    sigma_o = 0.318395
    eps_o = 0.21094 * 4.184  # Convert kcal/mol to kJ/mol
    
    ff.ljSigma = [sigma_o] * 16  # 4x4 matrix
    ff.ljEps = [0.0] * 16
    ff.ljEps[0] = eps_o  # O-O interaction
    
    state.forcefield = ff
    
    # Create water molecules in a grid
    atoms = []
    residues = []
    
    # Simple cubic arrangement
    n_per_dim = int(math.ceil(n_waters ** (1.0/3.0)))
    spacing = box_size / n_per_dim
    
    water_count = 0
    for i in range(n_per_dim):
        for j in range(n_per_dim):
            for k in range(n_per_dim):
                if water_count >= n_waters:
                    break
                    
                x = i * spacing + spacing/2
                y = j * spacing + spacing/2
                z = k * spacing + spacing/2
                
                water_atoms, water_res = create_swm4_water(x, y, z)
                
                # Update residue atom start
                water_res.atomStart = len(atoms)
                
                atoms.extend(water_atoms)
                residues.append(water_res)
                water_count += 1
    
    state.atoms = atoms
    state.residues = residues
    state.activeAtomCount = len(atoms)
    state.activeResidueCount = len(residues)
    
    return state


def calculate_drude_polarizability(charge, force_constant):
    """
    Calculate polarizability from charge and force constant
    
    α = q² / (4πε₀ * k)
    
    Args:
        charge: Drude particle charge
        force_constant: Harmonic force constant
        
    Returns:
        polarizability in atomic units
    """
    return charge * charge / (ONE_4PI_EPS0 * force_constant)


def setup_drude_from_psf(psf_file, drude_force=None):
    """
    Automatically set up DrudeForce from PSF topology information.
    
    This function:
    1. Identifies Drude particles by mass (~0.4 amu) and type containing 'D'
    2. Finds corresponding parent atoms (usually preceding atom)
    3. Uses alpha parameter from PSF to calculate polarizability
    4. Sets up Thole screening pairs based on PSF thole parameters
    
    Parameters
    ----------
    psf_file : str
        Path to PSF file containing Drude parameters
    drude_force : DrudeForce, optional
        Existing DrudeForce object. If None, uses global force.
    
    Returns
    -------
    dict
        Dictionary containing:
        - 'drude_particles': List of (drude_idx, parent_idx, polarizability)
        - 'screened_pairs': List of (idx1, idx2, thole)
        - 'topology': The parsed topology
    """
    # Parse PSF file
    topology = PSFParser.parse_file(psf_file)
    
    # Initialize DrudeForce if needed
    if drude_force is None:
        initializeDrudeForce()
    
    drude_particles = []
    parent_atoms = {}
    
    # 1. Identify Drude particles and their parents
    for i in range(topology.atom_count):
        atom = topology.get_atom(i)
        
        # Check if it's a Drude particle (mass ~0.4 and type contains 'D')
        if 0.3 < atom.get_mass() < 0.5 and 'D' in atom.get_type():
            # Parent is usually the previous atom
            if i > 0:
                parent_idx = i - 1
                parent = topology.get_atom(parent_idx)
                
                # Use parent's alpha parameter if available
                if abs(parent.get_alpha()) > 0.0:
                    # Convert alpha (Å³) to polarizability for force constant calculation
                    # In PSF, alpha is in Å³, we need to convert to appropriate units
                    # The force constant k = q²/(4πε₀α) where α is polarizability
                    
                    # Get Drude charge
                    drude_charge = atom.get_charge()
                    
                    # Convert alpha from Å³ to nm³ (1 Å = 0.1 nm)
                    alpha_nm3 = abs(parent.get_alpha()) * 0.001
                    
                    # Calculate polarizability from force constant
                    # k = q²/(4πε₀α) => α = q²/(4πε₀k)
                    # For typical Drude oscillators, k ~ 1000 kJ/mol/nm²
                    # We use alpha directly as polarizability parameter
                    polarizability = alpha_nm3
                    
                    # Add Drude particle
                    idx = addDrudeParticle(
                        drudeIndex=i,
                        parentIndex=parent_idx,
                        charge=drude_charge,
                        polarizability=polarizability
                    )
                    
                    drude_particles.append((i, parent_idx, polarizability))
                    parent_atoms[parent_idx] = parent.get_thole()
    
    # 2. Set up Thole screening pairs
    screened_pairs = []
    
    # For each pair of parent atoms with Drude particles
    parent_indices = list(parent_atoms.keys())
    for i in range(len(parent_indices)):
        for j in range(i + 1, len(parent_indices)):
            idx1 = parent_indices[i]
            idx2 = parent_indices[j]
            
            # Average thole parameters
            thole1 = parent_atoms[idx1]
            thole2 = parent_atoms[idx2]
            
            if thole1 > 0 or thole2 > 0:
                # Use average of non-zero thole parameters
                if thole1 > 0 and thole2 > 0:
                    thole_avg = (thole1 + thole2) / 2.0
                else:
                    thole_avg = max(thole1, thole2)
                
                # Find which Drude particles correspond to these parents
                drude_idx1 = None
                drude_idx2 = None
                
                for drude_idx, parent_idx, _ in drude_particles:
                    if parent_idx == idx1:
                        drude_idx1 = drude_idx
                    elif parent_idx == idx2:
                        drude_idx2 = drude_idx
                
                if drude_idx1 is not None and drude_idx2 is not None:
                    # Add screened pair between the two dipoles
                    # The dipole index in DrudeForce is the order they were added
                    dipole1 = next(i for i, (d, _, _) in enumerate(drude_particles) if d == drude_idx1)
                    dipole2 = next(i for i, (d, _, _) in enumerate(drude_particles) if d == drude_idx2)
                    
                    addDrudeScreenedPair(dipole1, dipole2, thole_avg)
                    screened_pairs.append((dipole1, dipole2, thole_avg))
    
    return {
        'drude_particles': drude_particles,
        'screened_pairs': screened_pairs,
        'topology': topology
    }


def create_system_from_psf_pdb(psf_file, pdb_file=None, box_dims=None):
    """
    Create a complete system from PSF and PDB files with automatic Drude setup.
    
    Parameters
    ----------
    psf_file : str
        Path to PSF file
    pdb_file : str, optional
        Path to PDB file (not used yet, requires PDB parser)
    box_dims : tuple of float, optional
        Box dimensions (x, y, z) in nm. Default is (10, 10, 10)
    
    Returns
    -------
    dict
        Dictionary with Drude setup details from PSF
    """
    # Set up Drude automatically from PSF
    drude_info = setup_drude_from_psf(psf_file)
    
    # Add box dimensions to the info
    if box_dims is None:
        box_dims = (10.0, 10.0, 10.0)
    
    drude_info['box_dims'] = box_dims
    
    # TODO: When PDB parser is available, read coordinates from PDB
    # and create MCState with proper atom positions
    
    return drude_info
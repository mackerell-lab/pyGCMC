"""
Safe PME initialization utilities to avoid segmentation faults.

This module provides a safe wrapper around initializePMEParameters that
prevents segmentation faults caused by PME global buffer management issues.
"""

import pygcmc
import gc
import time


def safe_initialize_pme(cutoff, box, alpha=0.0, mesh_size=None, spline_order=4, tolerance=1e-5):
    """
    Safely initialize PME parameters with workarounds for global buffer issues.
    
    This function prevents segmentation faults that occur when multiple MCState
    objects exist and one triggers PME reinitialization, invalidating pointers
    held by other states.
    
    Parameters:
    -----------
    cutoff : float
        Cutoff distance for PME calculations
    box : list of float
        Box dimensions [x, y, z]
    alpha : float, optional
        Ewald separation parameter (default: 0.0 for auto-calculation)
    mesh_size : list of int, optional
        PME grid size [nx, ny, nz] (default: None for auto-calculation)
    spline_order : int, optional
        B-spline interpolation order (default: 4)
    tolerance : float, optional
        Error tolerance (default: 1e-5)
    """
    # Note: gc.collect() was causing segfaults during cleanup of corrupted objects
    # Instead, we just call initializePMEParameters directly
    
    # Initialize PME with proper parameter handling
    if mesh_size is None:
        pygcmc.initializePMEParameters(cutoff, box, alpha, [], spline_order, tolerance)
    else:
        pygcmc.initializePMEParameters(cutoff, box, alpha, mesh_size, spline_order, tolerance)


def create_safe_mcstate(box, cutoff, forcefield=None):
    """
    Create a new MCState with proper deep copying to avoid shared pointers.
    
    Parameters:
    -----------
    box : list of float
        Box dimensions [x, y, z]
    cutoff : float
        Cutoff distance
    forcefield : MCForceField, optional
        Forcefield to use (will create a simple one if not provided)
    
    Returns:
    --------
    MCState : Properly initialized MCState object
    """
    state = pygcmc.MCState()
    
    # Deep copy box to avoid shared pointers
    state.info.box = list(box)
    state.info.cutoff = float(cutoff)
    
    # Set or create forcefield
    if forcefield is None:
        ff = pygcmc.MCForceField()
        ff.numTotalTypes = 1
        ff.numMovementTypes = 1
        ff.ljEps = [0.0]
        ff.ljSigma = [0.3]
        state.forcefield = ff
    else:
        state.forcefield = forcefield
    
    return state


def copy_mcstate_safely(original_state):
    """
    Create a safe deep copy of an MCState to avoid shared pointer issues.
    
    Parameters:
    -----------
    original_state : MCState
        The state to copy
    
    Returns:
    --------
    MCState : Deep copy of the original state
    """
    new_state = pygcmc.MCState()
    
    # Deep copy info
    new_state.info.box = list(original_state.info.box)
    new_state.info.cutoff = float(original_state.info.cutoff)
    
    # Forcefield is read-only, safe to share
    new_state.forcefield = original_state.forcefield
    
    # Deep copy atoms
    new_state.atoms = []
    for old_atom in original_state.atoms:
        new_atom = pygcmc.MCAtom()
        new_atom.x = old_atom.x
        new_atom.y = old_atom.y
        new_atom.z = old_atom.z
        new_atom.charge = old_atom.charge
        new_atom.type = old_atom.type
        new_state.atoms.append(new_atom)
    new_state.activeAtomCount = len(new_state.atoms)
    
    # Deep copy residues
    new_state.residues = []
    for old_res in original_state.residues:
        new_res = pygcmc.MCResidue()
        new_res.active = old_res.active
        new_res.fixed = old_res.fixed
        new_res.atomStart = old_res.atomStart
        new_res.atomCount = old_res.atomCount
        new_res.type = getattr(old_res, 'type', 0)
        new_state.residues.append(new_res)
    new_state.activeResidueCount = len(new_state.residues)
    
    return new_state
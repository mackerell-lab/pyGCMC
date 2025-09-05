"""Helper functions for movement tests to ensure isolation and reproducibility."""

import numpy as np
import pygcmc


def create_ideal_gas_state(base_state=None, box_size=(5.0, 5.0, 5.0), seed=42):
    """Create an ideal gas MCState with zero interactions.
    
    Args:
        base_state: Optional base state to use (will modify forcefield in place)
        box_size: Tuple of (Lx, Ly, Lz) for box dimensions
        seed: Random seed for reproducibility
        
    Returns:
        MCState configured as ideal gas (same object if base_state provided)
    """
    if base_state is not None:
        # Modify the existing state's forcefield for ideal gas
        state = base_state
        ff = state.forcefield
        
        # Zero out interactions if they exist
        if hasattr(ff, 'ljEps') and ff.ljEps:
            # Direct modification of existing forcefield
            for i in range(len(ff.ljEps)):
                ff.ljEps[i] = 0.0
        
        # Zero charges if they exist
        if hasattr(ff, 'charges') and ff.charges:
            for i in range(len(ff.charges)):
                ff.charges[i] = 0.0
    else:
        # Create fresh state - but this path may not work with current PyGCMC
        # Better to always pass a base_state
        state = base_state if base_state else pygcmc.MCState()
        if not hasattr(state, 'info'):
            state.info = pygcmc.SimInfo()
            state.info.box = box_size
            state.info.temperature = 298.15
    
    # Set random seed if needed
    if seed is not None:
        np.random.seed(seed)
    
    return state


def create_weak_interaction_state(base_state=None, epsilon_scale=0.01, box_size=(5.0, 5.0, 5.0), seed=42):
    """Create an MCState with weak interactions for testing.
    
    Args:
        base_state: Optional base state to use (will modify forcefield in place)
        epsilon_scale: Scaling factor for LJ epsilon (0.01 = 1% of original)
        box_size: Tuple of (Lx, Ly, Lz) for box dimensions
        seed: Random seed for reproducibility
        
    Returns:
        MCState with weak interactions (same object if base_state provided)
    """
    if base_state is not None:
        # Modify the existing state's forcefield for weak interactions
        state = base_state
        ff = state.forcefield
        
        # Scale down interactions if they exist
        if hasattr(ff, 'ljEps') and ff.ljEps:
            # Store original values for potential restoration
            original_eps = list(ff.ljEps)
            # Scale down epsilon
            for i in range(len(ff.ljEps)):
                ff.ljEps[i] = original_eps[i] * epsilon_scale
        
        # Zero charges if they exist
        if hasattr(ff, 'charges') and ff.charges:
            for i in range(len(ff.charges)):
                ff.charges[i] = 0.0
    else:
        # Use base_state if provided
        state = base_state if base_state else pygcmc.MCState()
        if not hasattr(state, 'info'):
            state.info = pygcmc.SimInfo()
            state.info.box = box_size
            state.info.temperature = 298.15
    
    # Set random seed if needed
    if seed is not None:
        np.random.seed(seed)
    
    return state


def save_forcefield(state):
    """Save forcefield parameters for later restoration.
    
    Args:
        state: MCState whose forcefield to save
        
    Returns:
        Dictionary of saved forcefield parameters
    """
    ff = state.forcefield
    saved = {}
    
    if hasattr(ff, 'ljEps'):
        saved['ljEps'] = list(ff.ljEps)
    if hasattr(ff, 'ljSigma'):
        saved['ljSigma'] = list(ff.ljSigma)
    if hasattr(ff, 'charges'):
        saved['charges'] = list(ff.charges)
    if hasattr(ff, 'masses'):
        saved['masses'] = list(ff.masses)
        
    return saved


def restore_forcefield(state, saved):
    """Restore forcefield parameters from saved values.
    
    Args:
        state: MCState whose forcefield to restore
        saved: Dictionary of saved forcefield parameters
    """
    ff = state.forcefield
    
    if 'ljEps' in saved and hasattr(ff, 'ljEps'):
        for i, val in enumerate(saved['ljEps']):
            if i < len(ff.ljEps):
                ff.ljEps[i] = val
    
    if 'ljSigma' in saved and hasattr(ff, 'ljSigma'):
        for i, val in enumerate(saved['ljSigma']):
            if i < len(ff.ljSigma):
                ff.ljSigma[i] = val
    
    if 'charges' in saved and hasattr(ff, 'charges'):
        for i, val in enumerate(saved['charges']):
            if i < len(ff.charges):
                ff.charges[i] = val
    
    if 'masses' in saved and hasattr(ff, 'masses'):
        for i, val in enumerate(saved['masses']):
            if i < len(ff.masses):
                ff.masses[i] = val


def residue_centroid_pbc(state, ridx):
    """Calculate centroid of residue with PBC wrapping.
    
    Args:
        state: MCState containing the residue
        ridx: Residue index
        
    Returns:
        numpy array [cx, cy, cz] wrapped to primary box, or None if invalid
    """
    if ridx < 0 or ridx >= len(state.residues):
        return None
        
    r = state.residues[ridx]
    if not hasattr(r, 'atomStart') or not hasattr(r, 'atomCount') or r.atomCount == 0:
        return None
        
    try:
        Lx, Ly, Lz = state.info.box
        xs = [state.atoms[i].x for i in range(r.atomStart, r.atomStart + r.atomCount)]
        ys = [state.atoms[i].y for i in range(r.atomStart, r.atomStart + r.atomCount)]
        zs = [state.atoms[i].z for i in range(r.atomStart, r.atomStart + r.atomCount)]
        
        # Calculate centroid with proper PBC wrapping
        cx = float(np.mean(xs)) % Lx
        cy = float(np.mean(ys)) % Ly
        cz = float(np.mean(zs)) % Lz
        
        return np.array([cx, cy, cz])
    except Exception:
        return None
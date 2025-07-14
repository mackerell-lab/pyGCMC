"""
PGP wrapper module for backward compatibility.
This module provides wrapper functions that use PGPContext internally.
"""

import pygcmc
from pygcmc import PGPContext

# Global PGPContext instance
_pgp_context = None
_pgp_initialized = False
_stored_pgp_params = {}
_stored_pme_params = {}

def reset_pgp_state():
    """Reset PGP state."""
    global _pgp_context, _pgp_initialized, _stored_pgp_params, _stored_pme_params
    _pgp_context = None
    _pgp_initialized = False
    _stored_pgp_params = {}
    _stored_pme_params = {}

def setPGPParameters(alpha, meshSize, potential_cutoff, potentialGridSize, splineOrder=4, tolerance=1e-5):
    """Set PGP parameters - pass through to original function."""
    # Just pass through to the original function
    pygcmc.setPGPParameters(alpha, meshSize, potential_cutoff, potentialGridSize, splineOrder, tolerance)

def setPMEParameters(alpha, meshSize, splineOrder=4, tolerance=1e-5):
    """Store PME parameters (for compatibility)."""
    # Call the original PME function for tests that use PME directly
    pygcmc.setPMEParameters(alpha, meshSize, splineOrder, tolerance)

def initializePMEParameters(cutoff, box, alpha):
    """Initialize PME parameters - pass through to original function."""
    # Just pass through to the original function
    # The C++ binding will handle all the initialization
    # IMPORTANT: Call with all 6 parameters to ensure the PME version is called
    # which properly sets global PME parameters
    pygcmc.initializePMEParameters(cutoff, box, alpha, [], 4, 1e-5)

def _ensure_pgp_initialized():
    """Ensure PGPContext is initialized with stored parameters."""
    global _pgp_context, _pgp_initialized
    
    # Debug print
    print(f"Debug: pgp_initialized = {_pgp_initialized}")
    print(f"Debug: global_pgp_context = {'null' if _pgp_context is None else 'set'}")
    print(f"Debug: stored_pme_params.set = {bool(_stored_pme_params)}")
    print(f"Debug: stored_pgp_params.set = {bool(_stored_pgp_params)}")
    
    if not _pgp_initialized or _pgp_context is None:
        if not _stored_pme_params or not _stored_pgp_params:
            raise RuntimeError("PGP not properly initialized. Call initializePMEParameters and setPGPParameters first.")
        
        _pgp_context = PGPContext()
        _pgp_context.initialize(
            cutoff=_stored_pme_params['cutoff'],
            box=_stored_pme_params['box'],
            alpha=_stored_pgp_params['alpha'],
            meshSize=_stored_pgp_params['meshSize'],
            potential_cutoff=_stored_pgp_params['potential_cutoff'],
            potentialGridSize=_stored_pgp_params['potentialGridSize'],
            splineOrder=_stored_pgp_params['splineOrder'],
            tolerance=_stored_pgp_params['tolerance']
        )
        _pgp_initialized = True

def precomputeGridPotential(state, fixed_only=True):
    """Precompute grid potential using standard PGP implementation."""
    # The C++ binding will handle initialization automatically
    # when it sees that both PME and PGP parameters have been set
    pygcmc.precomputeGridPotential(state, fixed_only)

def calculateMoleculeEnergy(state):
    """Calculate molecule energy using PGP interpolation."""
    # Call the interpolateMoleculeEnergy function directly
    return pygcmc.interpolateMoleculeEnergy(state)

def interpolateMoleculeEnergy(state):
    """Interpolate molecule energy - delegate to calculateMoleculeEnergy."""
    return calculateMoleculeEnergy(state)

def computeSystemEnergyPGP(state):
    """Compute system energy using standard PGP implementation."""
    # Use the real PGP implementation
    return pygcmc.computeSystemEnergyPGP(state)

def computeMovementEnergyPGP(state):
    """Compute movement energy using standard PGP implementation."""
    # Use the real PGP implementation
    return pygcmc.computeMovementEnergyPGP(state)


# Re-export original functions that don't need wrapping
from pygcmc import (
    resetPGPState,
    # PGP Complete functions
    computeSystemEnergyPGPComplete,
    computeMovementEnergyPGPComplete,
    # PME functions that tests might use for comparison
    computeSystemEnergyPME,
    computeMovementEnergyPME,
    computeSystemEnergyPMEComplete,
    computeMovementEnergyPMEFixed,
    computeSystemEnergyPMEFixed,
    # PGP Independent functions
    setPGPParametersIndependent,
    initializePGPParametersIndependent,
    computeSystemEnergyPGPIndependent,
    computeMovementEnergyPGPIndependent,
    # Ewald functions
    computeSystemEnergyEwald,
    # VDW functions
    computeSystemVdwEnergyCutoff
)

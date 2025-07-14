#!/usr/bin/env python3
"""
Script to convert PGP tests from global functions to PGPContext usage.
"""

import re
import sys
from pathlib import Path

def create_pgp_context_wrapper():
    """Create a wrapper module for backward compatibility."""
    wrapper_content = '''"""
PGP wrapper module for backward compatibility.
This module provides wrapper functions that use PGPContext internally.
"""

import pygcmc
from . import pgp_wrapper
from pygcmc import PGPContext
from .pgp_wrapper import (

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
    """Store PGP parameters for later initialization."""
    global _stored_pgp_params
    _stored_pgp_params = {
        'alpha': alpha,
        'meshSize': meshSize,
        'potential_cutoff': potential_cutoff,
        'potentialGridSize': potentialGridSize,
        'splineOrder': splineOrder,
        'tolerance': tolerance
    }

def setPMEParameters(alpha, meshSize, splineOrder=4, tolerance=1e-5):
    """Store PME parameters (for compatibility)."""
    # This is often called but not needed for PGPContext
    pass

def initializePMEParameters(cutoff, box, alpha):
    """Store PME parameters for PGPContext initialization."""
    global _stored_pme_params
    _stored_pme_params = {
        'cutoff': cutoff):
    """Ensure PGPContext is initialized with stored parameters."""
    global _pgp_context)
        
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
    """Precompute grid potential using PGPContext."""
    _ensure_pgp_initialized()
    
    # Precompute for all atom types
    for atom_type in range(state.forcefield.numTotalTypes):
        _pgp_context.precompute_grid_potential(state, atom_type=atom_type)

def calculateMoleculeEnergy(state):
    """Calculate molecule energy using PGPContext."""
    _ensure_pgp_initialized()
    return _pgp_context.compute_system_energy(state).total

def interpolateMoleculeEnergy(state):
    """Interpolate molecule energy using PGPContext."""
    return calculateMoleculeEnergy(state)

def computeSystemEnergyPGP(state):
    """Compute system energy using PGPContext."""
    _ensure_pgp_initialized()
    
    energy = _pgp_context.compute_system_energy(state)
    
    # Calculate VDW energy
    vdw = sum(res.energy_vdw for res in state.residues if res.active)
    
    # Calculate electrostatic energy
    electrostatic = energy.total - vdw
    
    # Create PGP dict for compatibility
    pgp_dict = {
        "real_space": energy.real_space,
        "reciprocal": energy.reciprocal,
        "self": energy.self,
        "total": energy.total
    }
    
    return electrostatic, vdw, pgp_dict

def computeMovementEnergyPGP(state):
    """Compute movement energy using PGPContext."""
    _ensure_pgp_initialized()
    
    # Extract movement indices
    movement_indices = []
    for info in state.movementResidues:
        for i in range(info.startIndex, info.startIndex + info.activeCount):
            movement_indices.append(i)
    
    energy = _pgp_context.compute_movement_energy(state, movement_indices)
    
    # Calculate VDW energy for movement residues only
    vdw = sum(state.residues[i].energy_vdw for i in movement_indices if state.residues[i].active)
    
    # Calculate electrostatic energy
    electrostatic = energy.total - vdw
    
    # Create PGP dict for compatibility
    pgp_dict = {
        "real_space": energy.real_space,
        "reciprocal": energy.reciprocal,
        "self": energy.self,
        "total": energy.total
    }
    
    return electrostatic, vdw, pgp_dict

# Re-export original functions that don't need wrapping
    computeSystemEnergyPGPComplete,
    computeMovementEnergyPGPComplete,
    resetPGPState
)
'''
    
    return wrapper_content

def fix_imports_in_file(file_path):
    """Fix imports in a test file to use the wrapper module."""
    
    with open(file_path, 'r') as f:
        content = f.read()
    
    original_content = content
    
    # Replace direct pygcmc PGP function imports
    patterns = [
        # Replace specific PGP function imports
        (r'from\s+pygcmc\s+import\s+\([^)]*?(setPGPParameters|precomputeGridPotential|computeSystemEnergyPGP|computeMovementEnergyPGP|calculateMoleculeEnergy|interpolateMoleculeEnergy)[^)]*\)',
         lambda m: m.group(0).replace('from pygcmc import', 'from .pgp_wrapper import')),
        
        # Replace pygcmc.function calls with wrapper calls
        (r'pygcmc\.(setPGPParameters|precomputeGridPotential|computeSystemEnergyPGP|computeMovementEnergyPGP|calculateMoleculeEnergy|interpolateMoleculeEnergy|setPMEParameters|initializePMEParameters)',
         r'pgp_wrapper.\1'),
        
        # Add wrapper import if using pygcmc.function syntax
        (r'(import pygcmc\n)(?!.*pgp_wrapper)',
         r'\1from . import pgp_wrapper
from .pgp_wrapper import (\n')
    ]
    
    # Apply replacements
    for pattern, replacement in patterns:
        if callable(replacement):
            content = re.sub(pattern, replacement, content)
        else:
            content = re.sub(pattern, replacement, content)
    
    # Write back if changed
    if content != original_content:
        with open(file_path, 'w') as f:
            f.write(content)
        return True
    
    return False

def main():
    """Main function to fix all PGP tests."""
    
    test_dir = Path("/home/zhaomt/gcmc/test107/pygcmc_dev/tests/simulation/energyPGP")
    
    # First, create the wrapper module
    wrapper_path = test_dir / "pgp_wrapper.py"
    with open(wrapper_path, 'w') as f:
        f.write(create_pgp_context_wrapper())
    print(f"Created wrapper module: {wrapper_path}")
    
    # Fix all test files
    skip_files = {'__init__.py', 'helpers.py', 'pgp_wrapper.py', 'README_PGP_COMPLETE_TESTS.md'}
    
    fixed_count = 0
    for py_file in test_dir.glob("*.py"):
        if py_file.name in skip_files:
            continue
        
        if fix_imports_in_file(py_file):
            print(f"Fixed imports in: {py_file.name}")
            fixed_count += 1
    
    print(f"\nFixed {fixed_count} files")
    print("\nNext steps:")
    print("1. Rebuild the project: cd build && make -j8")
    print("2. Run tests to verify: PYTHONPATH=$PYTHONPATH:./modules/bindings pytest ../tests/simulation/energyPGP/basic_operations.py -v")

if __name__ == "__main__":
    main()

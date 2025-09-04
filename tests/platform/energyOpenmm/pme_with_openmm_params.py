"""
Test PME with OpenMM-style parameter calculation

This test uses the same parameter calculation as OpenMM to verify
if the energy differences are reduced.
"""

import pytest
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pygcmc
from pygcmc import MCState, MCAtom, MCResidue, MCForceField
from pygcmc import initializePMEParameters, computeSystemEnergyPMEComplete

from energyOpenmm.pme_medium_complexity_helpers import create_medium_complexity_system, calculate_openmm_energy_components

try:
    from openmm import *
    from openmm.app import *
    from openmm.unit import *
    OPENMM_AVAILABLE = True
except ImportError:
    OPENMM_AVAILABLE = False


def calculate_openmm_style_params(cutoff, box_size, error_tol=5e-4):
    """Calculate PME parameters using OpenMM's formulas"""
    
    # Alpha calculation (from NonbondedForceImpl.cpp line 205)
    alpha = (1.0 / cutoff) * math.sqrt(-math.log(2.0 * error_tol))
    
    # Grid size calculation for electrostatics (from lines 212-214)
    # xsize = (int) ceil(2*alpha*boxVectors[0][0]/(3*pow(tol, 0.2)));
    grid_size = int(math.ceil(2 * alpha * box_size / (3 * math.pow(error_tol, 0.2))))
    grid_size = max(grid_size, 6)  # Minimum grid size is 6
    
    return alpha, grid_size


@pytest.mark.skipif(not OPENMM_AVAILABLE, reason="OpenMM not available")
def test_pme_with_openmm_params():
    """Test PME using OpenMM's parameter calculation"""
    
    # Create test system
    state = create_medium_complexity_system()
    
    print("\n" + "="*70)
    print("PME Test with OpenMM Parameters")
    print("="*70)
    
    # System info
    cutoff = state.info.cutoff
    box_size = state.info.box[0]
    
    print(f"\nSystem parameters:")
    print(f"  Atoms: {state.activeAtomCount}")
    print(f"  Box size: {box_size} nm")
    print(f"  Cutoff: {cutoff} nm")
    
    # Test with different error tolerances
    tolerances = [5e-4, 1e-4, 5e-5, 1e-5, 1e-6]
    
    print(f"\n{'Error Tol':>10} | {'Alpha':>8} | {'Grid':>6} | {'PyGCMC':>12} | {'OpenMM':>12} | {'Diff':>10} | {'Rel Error':>10}")
    print("-"*90)
    
    best_params = None
    best_error = float('inf')
    
    for tol in tolerances:
        # Calculate parameters using OpenMM's formula
        alpha, grid_size = calculate_openmm_style_params(cutoff, box_size, tol)
        
        # Use these parameters in PyGCMC
        mesh_size = [grid_size, grid_size, grid_size]
        spline_order = 4
        
        pygcmc.setPMEParameters(alpha, mesh_size, spline_order)
        initializePMEParameters(cutoff, state.info.box, alpha, mesh_size, spline_order)
        
        # Calculate with PyGCMC
        elec_pygcmc, vdw_pygcmc, total_pygcmc = computeSystemEnergyPMEComplete(state)
        
        # Calculate with OpenMM
        elec_openmm, vdw_openmm, total_openmm = calculate_openmm_energy_components(state, alpha)
        
        # Compare electrostatic energies (VdW should be 0 in this test)
        diff = elec_pygcmc - elec_openmm
        rel_error = abs(diff) / abs(elec_openmm) * 100 if elec_openmm != 0 else 0
        
        print(f"{tol:10.1e} | {alpha:8.3f} | {grid_size:6d} | {elec_pygcmc:12.4f} | {elec_openmm:12.4f} | {diff:10.4f} | {rel_error:9.3f}%")
        
        if rel_error < best_error:
            best_error = rel_error
            best_params = (tol, alpha, grid_size)
    
    print(f"\nBest parameters: tol={best_params[0]:.1e}, alpha={best_params[1]:.3f}, grid={best_params[2]}")
    print(f"Best relative error: {best_error:.3f}%")
    
    # Now test with exact OpenMM default
    print("\n" + "="*70)
    print("Using exact OpenMM defaults")
    print("="*70)
    
    # OpenMM's actual default from constructor
    default_tol = 5e-4
    alpha_default, grid_default = calculate_openmm_style_params(cutoff, box_size, default_tol)
    
    print(f"\nOpenMM default parameters:")
    print(f"  Error tolerance: {default_tol}")
    print(f"  Calculated alpha: {alpha_default:.3f}")
    print(f"  Calculated grid: {grid_default}")
    
    # Also try with the exact alpha that minimizes error from our previous test
    test_alphas = [2.0, 2.2, 2.4, 2.6, 2.8, 3.0, alpha_default]
    
    print(f"\n{'Alpha':>8} | {'Grid':>6} | {'PyGCMC':>12} | {'OpenMM':>12} | {'Diff':>10} | {'Rel Error':>10}")
    print("-"*80)
    
    for alpha in test_alphas:
        # Use a reasonable grid size
        grid_size = 32
        mesh_size = [grid_size, grid_size, grid_size]
        
        pygcmc.setPMEParameters(alpha, mesh_size, spline_order)
        initializePMEParameters(cutoff, state.info.box, alpha, mesh_size, spline_order)
        
        elec_pygcmc, vdw_pygcmc, total_pygcmc = computeSystemEnergyPMEComplete(state)
        elec_openmm, vdw_openmm, total_openmm = calculate_openmm_energy_components(state, alpha)
        
        diff = elec_pygcmc - elec_openmm
        rel_error = abs(diff) / abs(elec_openmm) * 100 if elec_openmm != 0 else 0
        
        print(f"{alpha:8.3f} | {grid_size:6d} | {elec_pygcmc:12.4f} | {elec_openmm:12.4f} | {diff:10.4f} | {rel_error:9.3f}%")


@pytest.mark.skipif(not OPENMM_AVAILABLE, reason="OpenMM not available")
def test_simple_two_charge_system():
    """Test with a simple two-charge system for clearer analysis"""
    
    print("\n" + "="*70)
    print("Simple Two-Charge System Test")
    print("="*70)
    
    # Create simple system
    state = MCState()
    state.info.box = [3.0, 3.0, 3.0]
    state.info.cutoff = 1.2
    
    ff = MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.ljEps = [0.0]
    ff.ljSigma = [0.1]
    state.forcefield = ff
    
    # Two opposite charges
    atoms = []
    
    atom1 = MCAtom()
    atom1.x, atom1.y, atom1.z = 1.5, 1.5, 1.2
    atom1.charge = 1.0
    atom1.type = 0
    atoms.append(atom1)
    
    atom2 = MCAtom()
    atom2.x, atom2.y, atom2.z = 1.5, 1.5, 1.8
    atom2.charge = -1.0
    atom2.type = 0
    atoms.append(atom2)
    
    residues = []
    for i in range(2):
        res = MCResidue()
        res.active = True
        res.fixed = True
        res.atomStart = i
        res.atomCount = 1
        res.type = 0
        residues.append(res)
    
    state.atoms = atoms
    state.activeAtomCount = 2
    state.residues = residues
    state.activeResidueCount = 2
    
    # Calculate OpenMM-style parameters
    alpha_openmm, grid_openmm = calculate_openmm_style_params(state.info.cutoff, state.info.box[0])
    
    print(f"\nOpenMM-style parameters:")
    print(f"  Alpha: {alpha_openmm:.3f}")
    print(f"  Grid: {grid_openmm}")
    
    # Test with these parameters
    mesh_size = [grid_openmm, grid_openmm, grid_openmm]
    spline_order = 4
    
    pygcmc.setPMEParameters(alpha_openmm, mesh_size, spline_order)
    initializePMEParameters(state.info.cutoff, state.info.box, alpha_openmm, mesh_size, spline_order)
    
    elec_pygcmc, _, _ = computeSystemEnergyPMEComplete(state)
    
    # Calculate with OpenMM
    if OPENMM_AVAILABLE:
        system = System()
        for atom in state.atoms:
            system.addParticle(1.0 * dalton)
        
        box = state.info.box
        system.setDefaultPeriodicBoxVectors(
            Vec3(box[0], 0, 0) * nanometer,
            Vec3(0, box[1], 0) * nanometer,
            Vec3(0, 0, box[2]) * nanometer
        )
        
        nonbonded = NonbondedForce()
        nonbonded.setNonbondedMethod(NonbondedForce.PME)
        nonbonded.setCutoffDistance(state.info.cutoff * nanometer)
        
        for atom in state.atoms:
            nonbonded.addParticle(
                atom.charge * elementary_charge,
                0.1 * nanometer,
                0.0 * kilojoule_per_mole
            )
        
        system.addForce(nonbonded)
        
        integrator = VerletIntegrator(1.0 * femtosecond)
        platform = Platform.getPlatformByName('Reference')
        context = Context(system, integrator, platform)
        
        positions = []
        for atom in state.atoms:
            positions.append(Vec3(atom.x, atom.y, atom.z) * nanometer)
        context.setPositions(positions)
        
        # Get actual PME parameters used by OpenMM
        actual_alpha, nx, ny, nz = nonbonded.getPMEParametersInContext(context)
        print(f"\nActual OpenMM PME parameters in context:")
        print(f"  Alpha: {actual_alpha}")
        print(f"  Grid: {nx}x{ny}x{nz}")
        
        energy_state = context.getState(getEnergy=True)
        elec_openmm = energy_state.getPotentialEnergy().value_in_unit(kilojoule_per_mole)
        
        print(f"\nEnergy comparison:")
        print(f"  PyGCMC: {elec_pygcmc:.6f} kJ/mol")
        print(f"  OpenMM: {elec_openmm:.6f} kJ/mol")
        print(f"  Difference: {elec_pygcmc - elec_openmm:.6f} kJ/mol")
        print(f"  Relative error: {abs(elec_pygcmc - elec_openmm) / abs(elec_openmm) * 100:.3f}%")


if __name__ == "__main__":
    test_pme_with_openmm_params()
    test_simple_two_charge_system()
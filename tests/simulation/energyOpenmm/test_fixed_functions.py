"""
Test the fixed energy calculation functions that correct LJ double-counting
"""

import pytest
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pygcmc
from pygcmc import MCState, MCAtom, MCResidue, MCForceField
from pygcmc import initializePMEParameters, computeSystemEnergyPME, computeSystemEnergyCutoff
from pygcmc import computeSystemEnergyPMEFixed, computeSystemEnergyCutoffFixed

# Only import OpenMM if available
try:
    from openmm import *
    from openmm.app import *
    from openmm.unit import *
    OPENMM_AVAILABLE = True
except ImportError:
    OPENMM_AVAILABLE = False


def create_lj_test_system():
    """Create a simple system with LJ interactions"""
    
    state = MCState()
    state.info.box = [4.0, 4.0, 4.0]
    state.info.cutoff = 1.5
    
    # Force field with LJ parameters
    ff = MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.ljEps = [1.0]     # kJ/mol
    ff.ljSigma = [0.35]  # nm
    state.forcefield = ff
    
    # Create two atoms at 0.5 nm distance
    atoms = []
    residues = []
    
    positions = [
        [2.0, 2.0, 2.0],
        [2.5, 2.0, 2.0]
    ]
    
    for i, pos in enumerate(positions):
        atom = MCAtom()
        atom.x, atom.y, atom.z = pos
        atom.charge = 0.0  # No charge for LJ-only test
        atom.type = 0
        atoms.append(atom)
        
        res = MCResidue()
        res.active = True
        res.fixed = True
        res.atomStart = i
        res.atomCount = 1
        res.type = 0
        residues.append(res)
    
    state.atoms = atoms
    state.activeAtomCount = len(atoms)
    state.residues = residues
    state.activeResidueCount = len(residues)
    
    return state


def calculate_expected_lj(r, sigma, epsilon):
    """Calculate expected LJ energy for a pair"""
    sr6 = (sigma/r)**6
    return 4.0 * epsilon * (sr6*sr6 - sr6)


def test_cutoff_fixed():
    """Test that the fixed cutoff function corrects double-counting"""
    
    state = create_lj_test_system()
    
    print("\nTesting cutoff energy calculation (fixed version):")
    
    # Calculate with original function
    computeSystemEnergyCutoff(state)
    
    # Get original energies
    original_vdw = 0.0
    for res in state.residues:
        if res.active:
            original_vdw += res.energy_vdw
    print(f"  Original VDW energy: {original_vdw:.6f} kJ/mol")
    
    # Reset and calculate with fixed function
    state = create_lj_test_system()
    elec, vdw, total = computeSystemEnergyCutoffFixed(state)
    
    print(f"  Fixed VDW energy: {vdw:.6f} kJ/mol")
    print(f"  Fixed electrostatic: {elec:.6f} kJ/mol")
    
    # Calculate expected
    expected = calculate_expected_lj(0.5, 0.35, 1.0)
    print(f"  Expected LJ energy: {expected:.6f} kJ/mol")
    
    # Check that fixed version gives correct result
    rel_diff = abs(vdw - expected) / abs(expected) if expected != 0 else 0
    print(f"  Relative difference: {rel_diff*100:.3f}%")
    
    assert rel_diff < 0.01, f"Fixed VDW energy differs from expected by {rel_diff*100:.3f}%"
    
    # Verify that original was indeed double-counted
    assert abs(original_vdw - 2*expected) < 0.001, "Original doesn't show double-counting"
    
    print("  ✓ Fixed cutoff function corrects double-counting!")


def test_pme_fixed():
    """Test that the fixed PME function corrects LJ double-counting"""
    
    # Create system with charges and LJ
    state = MCState()
    state.info.box = [4.0, 4.0, 4.0]
    state.info.cutoff = 1.5
    
    ff = MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.ljEps = [1.0]
    ff.ljSigma = [0.35]
    state.forcefield = ff
    
    atoms = []
    residues = []
    
    # Two atoms with charges
    positions = [
        ([2.0, 2.0, 2.0], 0.5),
        ([2.5, 2.0, 2.0], -0.5)
    ]
    
    for i, (pos, charge) in enumerate(positions):
        atom = MCAtom()
        atom.x, atom.y, atom.z = pos
        atom.charge = charge
        atom.type = 0
        atoms.append(atom)
        
        res = MCResidue()
        res.active = True
        res.fixed = True
        res.atomStart = i
        res.atomCount = 1
        res.type = 0
        residues.append(res)
    
    state.atoms = atoms
    state.activeAtomCount = len(atoms)
    state.residues = residues
    state.activeResidueCount = len(residues)
    
    print("\nTesting PME energy calculation (fixed version):")
    
    # Initialize PME
    alpha = 3.2
    mesh_size = [32, 32, 32]
    spline_order = 4
    
    pygcmc.setPMEParameters(alpha, mesh_size, spline_order)
    initializePMEParameters(state.info.cutoff, state.info.box, alpha, mesh_size, spline_order)
    
    # Calculate with original function
    state_copy = state  # Save for comparison
    computeSystemEnergyPME(state)
    original_vdw = 0.0
    for res in state.residues:
        if res.active:
            original_vdw += res.energy_vdw
    
    print(f"  Original VDW energy: {original_vdw:.6f} kJ/mol")
    
    # Calculate with fixed function
    elec, vdw, pme_dict = computeSystemEnergyPMEFixed(state_copy)
    
    print(f"  Fixed VDW energy: {vdw:.6f} kJ/mol")
    print(f"  Fixed electrostatic: {elec:.6f} kJ/mol")
    
    # Expected LJ
    expected_lj = calculate_expected_lj(0.5, 0.35, 1.0)
    print(f"  Expected LJ energy: {expected_lj:.6f} kJ/mol")
    
    # Check LJ correction
    rel_diff = abs(vdw - expected_lj) / abs(expected_lj) if expected_lj != 0 else 0
    print(f"  Relative difference: {rel_diff*100:.3f}%")
    
    assert rel_diff < 0.01, f"Fixed VDW energy differs from expected by {rel_diff*100:.3f}%"
    
    print("  ✓ Fixed PME function corrects LJ double-counting!")


@pytest.mark.skipif(not OPENMM_AVAILABLE, reason="OpenMM not available")
def test_fixed_vs_openmm():
    """Test that fixed functions match OpenMM"""
    
    print("\nTesting fixed functions vs OpenMM:")
    
    # Create test system
    state = create_lj_test_system()
    
    # Calculate with fixed function
    elec, vdw, total = computeSystemEnergyCutoffFixed(state)
    
    # Calculate with OpenMM
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
    nonbonded.setNonbondedMethod(NonbondedForce.CutoffPeriodic)
    nonbonded.setCutoffDistance(state.info.cutoff * nanometer)
    
    ff = state.forcefield
    for atom in state.atoms:
        nonbonded.addParticle(
            0.0 * elementary_charge,
            ff.ljSigma[0] * nanometer,
            ff.ljEps[0] * kilojoule_per_mole
        )
    
    system.addForce(nonbonded)
    
    integrator = VerletIntegrator(1.0 * femtosecond)
    platform = Platform.getPlatformByName('Reference')
    context = Context(system, integrator, platform)
    
    positions = []
    for atom in state.atoms:
        positions.append(Vec3(atom.x, atom.y, atom.z) * nanometer)
    context.setPositions(positions)
    
    energy_state = context.getState(getEnergy=True)
    openmm_energy = energy_state.getPotentialEnergy().value_in_unit(kilojoule_per_mole)
    
    print(f"  PyGCMC Fixed: {vdw:.6f} kJ/mol")
    print(f"  OpenMM: {openmm_energy:.6f} kJ/mol")
    
    diff = abs(vdw - openmm_energy)
    rel_diff = diff / abs(openmm_energy) if openmm_energy != 0 else 0
    
    print(f"  Difference: {rel_diff*100:.3f}%")
    
    assert rel_diff < 0.01, f"Fixed function differs from OpenMM by {rel_diff*100:.3f}%"
    
    print("  ✓ Fixed function matches OpenMM!")


if __name__ == "__main__":
    test_cutoff_fixed()
    test_pme_fixed()
    if OPENMM_AVAILABLE:
        test_fixed_vs_openmm()
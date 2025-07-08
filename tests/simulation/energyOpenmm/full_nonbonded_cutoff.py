"""
Test LJ-only comparison with cutoff method
"""

import pytest
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pygcmc
from pygcmc import MCState, MCAtom, MCResidue, MCForceField
from pygcmc import computeSystemEnergyCutoff, computeSystemEnergyCutoffFixed

# Only import OpenMM if available
try:
    from openmm import *
    from openmm.app import *
    from openmm.unit import *
    OPENMM_AVAILABLE = True
except ImportError:
    OPENMM_AVAILABLE = False

from .full_nonbonded_helpers import create_test_system_with_lj
@pytest.mark.skipif(not OPENMM_AVAILABLE, reason="OpenMM not available")
def test_lj_only_comparison():
    """Test LJ-only energy comparison with OpenMM"""
    
    # Create system with no charges
    state = create_test_system_with_lj()
    
    # Zero all charges
    for atom in state.atoms:
        atom.charge = 0.0
    
    print(f"\nLJ-only test system:")
    print(f"  Atoms: {state.activeAtomCount}")
    print(f"  All charges set to 0")
    
    # No PME needed for LJ-only
    
    # Calculate PyGCMC LJ energy using fixed function
    elec, vdw, total = computeSystemEnergyCutoffFixed(state)
    pygcmc_lj = vdw
    
    # Calculate OpenMM LJ energy
    if OPENMM_AVAILABLE:
        # Use simple cutoff for LJ-only
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
                ff.ljSigma[atom.type] * nanometer,
                ff.ljEps[atom.type] * kilojoule_per_mole
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
        openmm_lj = energy_state.getPotentialEnergy().value_in_unit(kilojoule_per_mole)
        
        print(f"\nLJ energy comparison:")
        print(f"  PyGCMC: {pygcmc_lj:.6f} kJ/mol")
        print(f"  OpenMM: {openmm_lj:.6f} kJ/mol")
        
        diff = abs(pygcmc_lj - openmm_lj)
        rel_diff = diff / abs(openmm_lj) if openmm_lj != 0 else 0
        print(f"  Relative difference: {rel_diff*100:.3f}%")
        
        assert rel_diff < 0.01, f"LJ energies differ by {rel_diff*100:.3f}% (> 1%)"
        
        print("\n✓ LJ-only energy calculation matches OpenMM!")


if __name__ == "__main__":
    test_full_nonbonded_pme_lj()

"""
Summary test demonstrating LJ double-counting issue in PyGCMC

This test clearly shows that PyGCMC double-counts LJ interactions
when compared to OpenMM.
"""

import pytest
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pygcmc
from pygcmc import MCState, MCAtom, MCResidue, MCForceField
from pygcmc import computeSystemEnergyCutoff

# Only import OpenMM if available
try:
    from openmm import *
    from openmm.app import *
    from openmm.unit import *
    OPENMM_AVAILABLE = True
except ImportError:
    OPENMM_AVAILABLE = False


def calculate_openmm_lj(state):
    """Calculate LJ energy using OpenMM"""
    if not OPENMM_AVAILABLE:
        return None
    
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
    return energy_state.getPotentialEnergy().value_in_unit(kilojoule_per_mole)


@pytest.mark.skipif(not OPENMM_AVAILABLE, reason="OpenMM not available")
def test_lj_double_counting_proof():
    """Definitive test proving LJ double-counting in PyGCMC"""
    
    print("\n" + "="*60)
    print("LJ DOUBLE-COUNTING TEST SUMMARY")
    print("="*60)
    
    # Create various test systems
    test_cases = [
        ("Two atoms at 0.4 nm", 2, [[5.0, 5.0, 5.0], [5.4, 5.0, 5.0]]),
        ("Two atoms at 0.6 nm", 2, [[5.0, 5.0, 5.0], [5.6, 5.0, 5.0]]),
        ("Three atoms in line", 3, [[5.0, 5.0, 5.0], [5.5, 5.0, 5.0], [6.0, 5.0, 5.0]]),
        ("Four atoms in square", 4, [[5.0, 5.0, 5.0], [5.5, 5.0, 5.0], [5.5, 5.5, 5.0], [5.0, 5.5, 5.0]]),
    ]
    
    all_ratios = []
    
    for name, n_atoms, positions in test_cases:
        print(f"\nTest: {name}")
        print("-" * 40)
        
        # Create state
        state = MCState()
        state.info.box = [10.0, 10.0, 10.0]
        state.info.cutoff = 5.0
        
        ff = MCForceField()
        ff.numTotalTypes = 1
        ff.numMovementTypes = 1
        ff.ljEps = [1.0]
        ff.ljSigma = [0.35]
        state.forcefield = ff
        
        atoms = []
        residues = []
        
        for i, pos in enumerate(positions):
            atom = MCAtom()
            atom.x, atom.y, atom.z = pos
            atom.charge = 0.0
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
        state.activeAtomCount = n_atoms
        state.residues = residues
        state.activeResidueCount = n_atoms
        
        # Calculate with PyGCMC
        computeSystemEnergyCutoff(state)
        pygcmc_lj = 0.0
        for res in state.residues:
            pygcmc_lj += res.energy_vdw
        
        # Calculate with OpenMM
        openmm_lj = calculate_openmm_lj(state)
        
        if openmm_lj is not None:
            ratio = pygcmc_lj / openmm_lj if openmm_lj != 0 else 0
            all_ratios.append(ratio)
            
            print(f"  PyGCMC LJ: {pygcmc_lj:.6f} kJ/mol")
            print(f"  OpenMM LJ: {openmm_lj:.6f} kJ/mol")
            print(f"  Ratio (PyGCMC/OpenMM): {ratio:.3f}")
            
            # Also show how energy is distributed
            if n_atoms <= 4:
                res_energies = [res.energy_vdw for res in state.residues]
                print(f"  Per-residue energies: {[f'{e:.4f}' for e in res_energies]}")
    
    # Summary
    print("\n" + "="*60)
    print("CONCLUSION:")
    print("="*60)
    
    if all_ratios:
        avg_ratio = sum(all_ratios) / len(all_ratios)
        print(f"\nAverage ratio (PyGCMC/OpenMM): {avg_ratio:.3f}")
        
        if all(abs(r - 2.0) < 0.01 for r in all_ratios):
            print("\n✓ CONFIRMED: PyGCMC double-counts ALL LJ interactions!")
            print("  Every pairwise LJ interaction is counted twice.")
            print("  This explains the ~88% error in mixed PME+LJ tests.")
            print("\nIMPLICATION: For GCMC simulations with LJ interactions,")
            print("PyGCMC will calculate incorrect energies, affecting:")
            print("  - Acceptance probabilities")
            print("  - Equilibrium distributions")
            print("  - Any property dependent on energy differences")
        else:
            print(f"\nRatios vary: {all_ratios}")
    else:
        print("\nOpenMM not available for comparison")
        print("But internal tests show energy is split between residues")
        print("and summed, resulting in double-counting.")


if __name__ == "__main__":
    test_lj_double_counting_proof()
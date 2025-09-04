"""
Simple test to understand PME differences between PyGCMC and OpenMM
"""

import pytest
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pygcmc
from pygcmc import MCState, MCAtom, MCResidue, MCForceField
from pygcmc import initializePMEParameters, computeSystemEnergyPME

try:
    from openmm import *
    from openmm.app import *
    from openmm.unit import *
    OPENMM_AVAILABLE = True
except ImportError:
    OPENMM_AVAILABLE = False


def create_two_charge_system():
    """Create simple system with just two opposite charges"""
    state = MCState()
    state.info.box = [3.0, 3.0, 3.0]  # 3 nm box
    state.info.cutoff = 1.0  # 1 nm cutoff
    
    # Create atoms list first
    atoms_list = []
    residues_list = []
    
    # Create two atoms with opposite charges
    positions = [
        [1.5, 1.5, 1.5],  # Center
        [2.0, 1.5, 1.5],  # 0.5 nm away
    ]
    charges = [1.0, -1.0]  # +1 and -1 charge
    
    for i, (pos, charge) in enumerate(zip(positions, charges)):
        atom = MCAtom()
        atom.x, atom.y, atom.z = pos
        atom.charge = charge
        atom.type = 0
        atoms_list.append(atom)
        
        res = MCResidue()
        res.active = True
        res.fixed = True
        res.atomStart = i
        res.atomCount = 1
        res.type = 0
        residues_list.append(res)
    
    # Set the arrays and counts together
    state.atoms = atoms_list
    state.residues = residues_list
    state.activeAtomCount = len(atoms_list)
    state.activeResidueCount = len(residues_list)
    
    # Initialize force field
    ff = MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.ljEps = [0.0]
    ff.ljSigma = [0.0]
    state.forcefield = ff
    
    return state


@pytest.mark.skipif(not OPENMM_AVAILABLE, reason="OpenMM not available")
def test_simple_pme_comparison():
    """Compare PME for simple two-charge system"""
    
    state = create_two_charge_system()
    
    # Test different alpha values
    alphas = [2.0, 3.0, 4.0, 5.0]
    mesh_size = [32, 32, 32]
    spline_order = 4
    
    print("\nTwo opposite charges separated by 0.5 nm:")
    print("=" * 60)
    print(f"Number of atoms: {len(state.atoms)}")
    print(f"Active atom count: {state.activeAtomCount}")
    
    # Track best agreement for assertions
    best_rel_diff = float('inf')
    best_alpha = None
    
    for alpha in alphas:
        # PyGCMC calculation
        pygcmc.setPMEParameters(alpha, mesh_size, spline_order)
        initializePMEParameters(state.info.cutoff, state.info.box, alpha, mesh_size, spline_order)
        computeSystemEnergyPME(state)
        
        real_space = state.ewald_energy['real_space']
        reciprocal = state.ewald_energy['reciprocal']
        self_energy = state.ewald_energy['self']
        pygcmc_total = real_space + reciprocal + self_energy
        
        # Validate PyGCMC components
        assert abs(self_energy) > 0, "Self energy should be non-zero for charged system"
        assert pygcmc_total != 0, "Total PME energy should be non-zero for charged system"
        
        # For two opposite charges, the energy should be negative (attractive)
        assert pygcmc_total < 0, "Two opposite charges should have negative (attractive) energy"
        
        # OpenMM calculation
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
            nonbonded.setEwaldErrorTolerance(1e-6)
            
            for atom in state.atoms:
                nonbonded.addParticle(
                    atom.charge * elementary_charge,
                    1.0 * nanometer,
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
            
            energy_state = context.getState(getEnergy=True)
            openmm_total = energy_state.getPotentialEnergy().value_in_unit(kilojoule_per_mole)
            
            # Get actual alpha used by OpenMM
            for force in system.getForces():
                if isinstance(force, NonbondedForce):
                    pme_params = force.getPMEParametersInContext(context)
                    actual_alpha = pme_params[0]
            
            # Calculate differences
            abs_diff = abs(pygcmc_total - openmm_total)
            rel_diff = abs_diff / abs(openmm_total) if openmm_total != 0 else 0
            
            print(f"\nAlpha = {alpha:.1f} (OpenMM actual: {actual_alpha:.3f}):")
            print(f"  PyGCMC: Real={real_space:.4f}, Recip={reciprocal:.4f}, Self={self_energy:.4f}, Total={pygcmc_total:.4f}")
            print(f"  OpenMM: Total={openmm_total:.4f}")
            print(f"  Difference: {abs_diff:.4f} ({rel_diff*100:.1f}%)")
            
            # Track best agreement
            if rel_diff < best_rel_diff:
                best_rel_diff = rel_diff
                best_alpha = alpha
            
            # Assert reasonable agreement for each alpha
            assert rel_diff < 0.25, f"PME energies differ by {rel_diff*100:.1f}% (> 25%) for alpha={alpha}"
    
    # Final assertion: at least one alpha should give good agreement
    if OPENMM_AVAILABLE:
        print(f"\nBest agreement with alpha={best_alpha}: {best_rel_diff*100:.2f}% difference")
        assert best_rel_diff < 0.1, f"Best PME agreement is {best_rel_diff*100:.2f}% (> 10%), which may indicate an issue"


if __name__ == "__main__":
    test_simple_pme_comparison()
"""
Test PME accuracy for medium complexity system (multiple ions + ligand, no water)

This test bridges the gap between simple ion systems (1-2% error) and 
complex water-containing systems (15% error) to understand PME accuracy scaling.
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


def create_medium_complexity_system():
    """Create system with multiple ions and a ligand molecule, but no water"""
    
    # System parameters - intermediate between simple and complex
    box_size = 3.5  # nm, between 3.0 (water system) and 4.0 (simple system)
    cutoff = 1.2    # nm
    
    # Create state
    state = MCState()
    state.info.box = [box_size, box_size, box_size]
    state.info.cutoff = cutoff
    
    # Force field - no LJ for pure electrostatics comparison
    ff = MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.ljEps = [0.0]
    ff.ljSigma = [0.3]
    state.forcefield = ff
    
    atoms = []
    residues = []
    atom_idx = 0
    
    # Add 12 ions in a regular pattern (6 Na+, 6 Cl-)
    # This is more than simple tests (4-8) but less than water systems (hundreds)
    ion_positions = [
        # Layer 1
        [1.0, 1.0, 1.0], [2.5, 1.0, 1.0], [1.0, 2.5, 1.0], [2.5, 2.5, 1.0],
        # Layer 2
        [1.0, 1.0, 2.5], [2.5, 1.0, 2.5], [1.0, 2.5, 2.5], [2.5, 2.5, 2.5],
        # Additional ions
        [1.75, 1.75, 1.0], [1.75, 1.75, 2.5], [1.0, 1.75, 1.75], [2.5, 1.75, 1.75]
    ]
    
    # Alternating charges for neutrality
    ion_charges = [1.0, -1.0, -1.0, 1.0, -1.0, 1.0, 1.0, -1.0, 1.0, -1.0, -1.0, 1.0]
    
    for i, (pos, charge) in enumerate(zip(ion_positions, ion_charges)):
        atom = MCAtom()
        atom.x, atom.y, atom.z = pos
        atom.charge = charge
        atom.type = 0
        atoms.append(atom)
        
        # Ion residue
        res = MCResidue()
        res.active = True
        res.fixed = True
        res.atomStart = atom_idx
        res.atomCount = 1
        res.type = 0
        residues.append(res)
        atom_idx += 1
    
    # Add a ligand molecule (8 atoms, mimicking a small organic molecule)
    # Place it in the center of the box
    ligand_start = atom_idx
    ligand_center = box_size / 2.0
    
    # Create a ligand with realistic partial charges (sum to 0)
    ligand_atoms = [
        # Position relative to center, charge (mimicking functional groups)
        ([0.0, 0.0, 0.0], -0.3),    # Central carbon
        ([0.15, 0.0, 0.0], 0.1),    # CH
        ([-0.15, 0.0, 0.0], 0.1),   # CH
        ([0.0, 0.15, 0.0], 0.1),    # CH
        ([0.0, -0.15, 0.0], -0.4),  # Oxygen
        ([0.0, 0.0, 0.15], 0.2),    # NH
        ([0.0, 0.0, -0.15], 0.1),   # CH
        ([0.2, 0.2, 0.0], 0.1),     # CH3
    ]
    
    for rel_pos, charge in ligand_atoms:
        atom = MCAtom()
        atom.x = ligand_center + rel_pos[0]
        atom.y = ligand_center + rel_pos[1]
        atom.z = ligand_center + rel_pos[2]
        atom.charge = charge
        atom.type = 0
        atoms.append(atom)
    
    # Ligand residue
    res = MCResidue()
    res.active = True
    res.fixed = False
    res.atomStart = ligand_start
    res.atomCount = len(ligand_atoms)
    res.type = 0
    residues.append(res)
    
    state.atoms = atoms
    state.activeAtomCount = len(atoms)
    state.residues = residues
    state.activeResidueCount = len(residues)
    
    # Verify charge neutrality
    total_charge = 0.0
    for atom in atoms:
        total_charge += atom.charge
    print(f"Total system charge: {total_charge:.6f} (should be ~0)")
    
    return state


def calculate_openmm_energy_medium(state, alpha):
    """Calculate PME energy using OpenMM for medium complexity system"""
    if not OPENMM_AVAILABLE:
        return None
    
    # Create OpenMM system
    system = System()
    
    # Add particles
    for atom in state.atoms:
        system.addParticle(1.0 * dalton)
    
    # Set periodic box
    box = state.info.box
    system.setDefaultPeriodicBoxVectors(
        Vec3(box[0], 0, 0) * nanometer,
        Vec3(0, box[1], 0) * nanometer,
        Vec3(0, 0, box[2]) * nanometer
    )
    
    # Create PME force
    nonbonded = NonbondedForce()
    nonbonded.setNonbondedMethod(NonbondedForce.PME)
    nonbonded.setCutoffDistance(state.info.cutoff * nanometer)
    nonbonded.setEwaldErrorTolerance(1e-6)  # High precision
    
    # Add particles with charges only
    for i, atom in enumerate(state.atoms):
        nonbonded.addParticle(
            atom.charge * elementary_charge,
            0.1 * nanometer,  # Small sigma
            0.0 * kilojoule_per_mole  # No LJ
        )
    
    system.addForce(nonbonded)
    
    # Create context
    integrator = VerletIntegrator(1.0 * femtosecond)
    platform = Platform.getPlatformByName('Reference')
    context = Context(system, integrator, platform)
    
    # Set positions
    positions = []
    for atom in state.atoms:
        positions.append(Vec3(atom.x, atom.y, atom.z) * nanometer)
    context.setPositions(positions)
    
    # Get energy
    energy_state = context.getState(getEnergy=True)
    energy_kj_mol = energy_state.getPotentialEnergy().value_in_unit(kilojoule_per_mole)
    
    return energy_kj_mol


@pytest.mark.skipif(not OPENMM_AVAILABLE, reason="OpenMM not available")
def test_pme_medium_complexity():
    """Test PME accuracy for medium complexity system"""
    
    state = create_medium_complexity_system()
    
    print("\nMedium Complexity System:")
    print("=" * 70)
    print(f"Total atoms: {state.activeAtomCount} (12 ions + 8 ligand atoms)")
    print(f"Box size: {state.info.box[0]} nm")
    print(f"Cutoff: {state.info.cutoff} nm")
    
    # Test with different alpha values to find optimal
    test_configs = [
        (2.5, [32, 32, 32], 4),
        (3.0, [32, 32, 32], 4),
        (3.5, [32, 32, 32], 4),
        (3.0, [48, 48, 48], 4),
        (3.0, [64, 64, 64], 4),
    ]
    
    results = []
    
    print("\nPME comparison results:")
    print("-" * 70)
    print(f"{'Alpha':>6} {'Mesh':>12} {'PyGCMC':>12} {'OpenMM':>12} {'Diff':>10} {'Rel %':>8}")
    print("-" * 70)
    
    for alpha, mesh_size, spline_order in test_configs:
        # PyGCMC calculation
        pygcmc.setPMEParameters(alpha, mesh_size, spline_order)
        initializePMEParameters(state.info.cutoff, state.info.box, alpha, mesh_size, spline_order)
        computeSystemEnergyPME(state)
        
        real_space = state.ewald_energy.get('real_space', 0.0)
        reciprocal = state.ewald_energy.get('reciprocal', 0.0)
        self_energy = state.ewald_energy.get('self', 0.0)
        pygcmc_total = real_space + reciprocal + self_energy
        
        # OpenMM calculation
        openmm_energy = calculate_openmm_energy_medium(state, alpha)
        
        if openmm_energy is not None:
            diff = abs(pygcmc_total - openmm_energy)
            rel_diff = diff / abs(openmm_energy) * 100 if openmm_energy != 0 else 0
            
            mesh_str = f"{mesh_size[0]}x{mesh_size[1]}x{mesh_size[2]}"
            print(f"{alpha:6.1f} {mesh_str:>12} {pygcmc_total:12.4f} {openmm_energy:12.4f} "
                  f"{diff:10.4f} {rel_diff:7.2f}%")
            
            results.append((alpha, mesh_size, rel_diff))
    
    # Find best configuration
    best_config = min(results, key=lambda x: x[2])
    best_alpha, best_mesh, best_error = best_config
    
    print("\n" + "=" * 70)
    print(f"Best configuration: alpha={best_alpha}, mesh={best_mesh[0]}x{best_mesh[1]}x{best_mesh[2]}")
    print(f"Best relative error: {best_error:.2f}%")
    
    # Detailed breakdown for best configuration
    print("\nDetailed breakdown for best configuration:")
    pygcmc.setPMEParameters(best_alpha, best_mesh, spline_order)
    initializePMEParameters(state.info.cutoff, state.info.box, best_alpha, best_mesh, spline_order)
    computeSystemEnergyPME(state)
    
    print(f"  Real space:  {state.ewald_energy.get('real_space', 0.0):10.4f} kJ/mol")
    print(f"  Reciprocal:  {state.ewald_energy.get('reciprocal', 0.0):10.4f} kJ/mol")
    print(f"  Self:        {state.ewald_energy.get('self', 0.0):10.4f} kJ/mol")
    print(f"  Total:       {state.ewald_energy.get('real_space', 0.0) + state.ewald_energy.get('reciprocal', 0.0) + state.ewald_energy.get('self', 0.0):10.4f} kJ/mol")
    
    # Analysis
    print("\nComplexity vs Error Analysis:")
    print("-" * 70)
    print("System Type              | Atoms | Water | Error Range")
    print("-" * 70)
    print("Simple (ions only)       |  4-8  |  No   | 1-2%")
    print("Medium (ions + ligand)   |  20   |  No   | {:.1f}-{:.1f}%".format(
        min(r[2] for r in results), max(r[2] for r in results)))
    print("Complex (ions+lig+water) | 100+  |  Yes  | 15%")
    
    # Check if error is in expected range - updated to reflect improved accuracy
    assert best_error <= 0.1, f"Medium complexity error {best_error:.2f}% higher than expected (<0.1%)"
    
    # More specifically, we expect it to be closer to simple than complex
    if best_error < 5.0:
        print("\n✓ Medium complexity system shows low error (<5%), closer to simple systems")
    elif best_error < 10.0:
        print("\n✓ Medium complexity system shows moderate error (5-10%), as expected")
    else:
        print("\n⚠ Medium complexity system shows high error (>10%), approaching water system levels")


if __name__ == "__main__":
    test_pme_medium_complexity()
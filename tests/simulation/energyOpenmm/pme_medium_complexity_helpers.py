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


def calculate_openmm_energy_components(state, alpha):
    """Calculate PME energy components separately using OpenMM"""
    if not OPENMM_AVAILABLE:
        return None, None, None
    
    # Create two systems - one for electrostatics, one for LJ
    system_elec = System()
    system_lj = System()
    
    # Add particles to both
    for atom in state.atoms:
        system_elec.addParticle(1.0 * dalton)
        system_lj.addParticle(1.0 * dalton)
    
    # Set periodic box for both
    box = state.info.box
    box_vectors = [
        Vec3(box[0], 0, 0) * nanometer,
        Vec3(0, box[1], 0) * nanometer,
        Vec3(0, 0, box[2]) * nanometer
    ]
    system_elec.setDefaultPeriodicBoxVectors(*box_vectors)
    system_lj.setDefaultPeriodicBoxVectors(*box_vectors)
    
    # Create forces
    nonbonded_elec = NonbondedForce()
    nonbonded_lj = NonbondedForce()
    
    # Configure PME for both
    for nb in [nonbonded_elec, nonbonded_lj]:
        nb.setNonbondedMethod(NonbondedForce.PME)
        nb.setCutoffDistance(state.info.cutoff * nanometer)
        nb.setEwaldErrorTolerance(1e-6)
    
    # Add particles with appropriate parameters
    for i, atom in enumerate(state.atoms):
        # Electrostatics only
        nonbonded_elec.addParticle(
            atom.charge * elementary_charge,
            0.1 * nanometer,  # Small sigma
            0.0 * kilojoule_per_mole  # No LJ
        )
        # LJ only
        nonbonded_lj.addParticle(
            0.0 * elementary_charge,  # No charge
            0.1 * nanometer,  # Small sigma (since ljEps is 0 anyway)
            0.0 * kilojoule_per_mole  # No LJ in this test system
        )
    
    system_elec.addForce(nonbonded_elec)
    system_lj.addForce(nonbonded_lj)
    
    # Create contexts
    integrator_elec = VerletIntegrator(1.0 * femtosecond)
    integrator_lj = VerletIntegrator(1.0 * femtosecond)
    platform = Platform.getPlatformByName('Reference')
    context_elec = Context(system_elec, integrator_elec, platform)
    context_lj = Context(system_lj, integrator_lj, platform)
    
    # Set positions
    positions = []
    for atom in state.atoms:
        positions.append(Vec3(atom.x, atom.y, atom.z) * nanometer)
    context_elec.setPositions(positions)
    context_lj.setPositions(positions)
    
    # Get energies
    energy_state_elec = context_elec.getState(getEnergy=True)
    energy_state_lj = context_lj.getState(getEnergy=True)
    
    elec_energy = energy_state_elec.getPotentialEnergy().value_in_unit(kilojoule_per_mole)
    lj_energy = energy_state_lj.getPotentialEnergy().value_in_unit(kilojoule_per_mole)
    total_energy = elec_energy + lj_energy
    
    return elec_energy, lj_energy, total_energy



"""
Test PME energy calculations for ligand insertion in charged systems

This test validates that PyGCMC correctly calculates PME electrostatic energies
for ligand insertion scenarios, which is critical for GCMC simulations.
"""

import pytest
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pygcmc
from pygcmc import MCState, MCAtom, MCResidue, MCForceField, MCMovementResidueInfo
from pygcmc import initializePMEParameters, computeSystemEnergyPME, computeMovementEnergyPME

# Only import OpenMM if available
try:
    from openmm import *
    from openmm.app import *
    from openmm.unit import *
    OPENMM_AVAILABLE = True
except ImportError:
    OPENMM_AVAILABLE = False


def create_protein_like_system():
    """Create a system mimicking a protein binding site with water and ions"""
    
    # Larger box to simulate protein environment
    box_size = 6.0  # nm
    cutoff = 1.8    # nm, typical for protein simulations
    
    state = MCState()
    state.info.box = [box_size, box_size, box_size]
    state.info.cutoff = cutoff
    
    # Force field
    ff = MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.ljEps = [0.0]  # Pure electrostatics
    ff.ljSigma = [0.3]
    state.forcefield = ff
    
    atoms = []
    residues = []
    atom_idx = 0
    
    # Create a binding site with charged residues
    # Mimic lysine (positive) and glutamate (negative) side chains
    binding_site_atoms = [
        # Lysine-like (positive)
        ([2.5, 3.0, 3.0], 1.0),   # NZ
        ([3.5, 3.0, 3.0], 1.0),   # NZ
        # Glutamate-like (negative)
        ([3.0, 2.5, 3.0], -0.5),  # OE1
        ([3.0, 3.5, 3.0], -0.5),  # OE2
        ([3.0, 3.0, 2.5], -0.5),  # OE1
        ([3.0, 3.0, 3.5], -0.5),  # OE2
    ]
    
    for pos, charge in binding_site_atoms:
        atom = MCAtom()
        atom.x, atom.y, atom.z = pos
        atom.charge = charge
        atom.type = 0
        atoms.append(atom)
        
        res = MCResidue()
        res.active = True
        res.fixed = True
        res.atomStart = atom_idx
        res.atomCount = 1
        res.type = 0
        residues.append(res)
        atom_idx += 1
    
    # Add water molecules around the binding site
    water_positions = [
        [2.0, 2.0, 3.0], [4.0, 4.0, 3.0],
        [2.0, 4.0, 3.0], [4.0, 2.0, 3.0],
        [3.0, 3.0, 4.0], [3.0, 3.0, 2.0],
    ]
    
    for pos in water_positions:
        # Simplified water as single point with TIP3P oxygen charge
        atom = MCAtom()
        atom.x, atom.y, atom.z = pos
        atom.charge = -0.834
        atom.type = 0
        atoms.append(atom)
        
        res = MCResidue()
        res.active = True
        res.fixed = True
        res.atomStart = atom_idx
        res.atomCount = 1
        res.type = 0
        residues.append(res)
        atom_idx += 1
    
    # Add ions in the bulk
    ion_positions = [
        ([1.0, 1.0, 1.0], 1.0),   # Na+
        ([5.0, 5.0, 5.0], -1.0),  # Cl-
        ([1.0, 5.0, 1.0], 1.0),   # Na+
        ([5.0, 1.0, 5.0], -1.0),  # Cl-
    ]
    
    for pos, charge in ion_positions:
        atom = MCAtom()
        atom.x, atom.y, atom.z = pos
        atom.charge = charge
        atom.type = 0
        atoms.append(atom)
        
        res = MCResidue()
        res.active = True
        res.fixed = True
        res.atomStart = atom_idx
        res.atomCount = 1
        res.type = 0
        residues.append(res)
        atom_idx += 1
    
    state.atoms = atoms
    state.activeAtomCount = len(atoms)
    state.residues = residues
    state.activeResidueCount = len(residues)
    
    return state, atom_idx


def add_ligand_to_system(state, ligand_position, ligand_charges):
    """Add a ligand to the system at specified position"""
    
    start_idx = state.activeAtomCount
    
    # Add ligand atoms
    for i, charge in enumerate(ligand_charges):
        atom = MCAtom()
        # Place atoms in a line with 0.15 nm spacing
        atom.x = ligand_position[0] + i * 0.15
        atom.y = ligand_position[1]
        atom.z = ligand_position[2]
        atom.charge = charge
        atom.type = 0
        state.atoms.append(atom)
    
    # Add ligand residue
    res = MCResidue()
    res.active = True
    res.fixed = False
    res.atomStart = start_idx
    res.atomCount = len(ligand_charges)
    res.type = 0
    state.residues.append(res)
    
    state.activeAtomCount = len(state.atoms)
    state.activeResidueCount = len(state.residues)
    
    # Set as movement residue
    movement_info = MCMovementResidueInfo()
    movement_info.startIndex = res.atomStart
    movement_info.activeCount = res.atomCount
    state.movementResidues = [movement_info]
    
    return res


def calculate_openmm_pme_for_system(state, alpha):
    """Calculate PME energy using OpenMM with specific alpha"""
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
    nonbonded.setNonbondedMethod(NonbondedForce.PME)
    nonbonded.setCutoffDistance(state.info.cutoff * nanometer)
    nonbonded.setEwaldErrorTolerance(1e-6)
    
    for i, atom in enumerate(state.atoms):
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
    
    energy_state = context.getState(getEnergy=True)
    return energy_state.getPotentialEnergy().value_in_unit(kilojoule_per_mole)

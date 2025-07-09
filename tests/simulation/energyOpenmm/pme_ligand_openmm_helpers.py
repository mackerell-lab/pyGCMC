"""
PME energy comparison with OpenMM for ligands in charged system

This test places multiple ligands in a system with charged molecules (like ions)
and compares the PME electrostatic energy calculation between pygcmc and OpenMM.
This is a critical test for GCMC simulations where ligands are inserted into
protein-water-ion systems.
"""

# Helper functions for PME ligand energy comparison tests

import pytest
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pygcmc
from pygcmc import MCState, MCAtom, MCResidue, MCForceField
from pygcmc import initializePMEParameters, computeSystemEnergyPME

# Only import OpenMM if available
try:
    from openmm import *
    from openmm.app import *
    from openmm.unit import *
    OPENMM_AVAILABLE = True
except ImportError:
    OPENMM_AVAILABLE = False


def create_charged_system_with_ligands():
    """Create a system with ions and water-like molecules, then add ligands"""
    
    # System parameters
    box_size = 3.0  # nm, small box for testing
    cutoff = 1.2    # nm
    
    # Create state
    state = MCState()
    state.info.box = [box_size, box_size, box_size]
    state.info.cutoff = cutoff
    
    # Force field parameters - simplified to single type for now
    ff = MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.ljEps = [0.0]  # No LJ, only electrostatics
    ff.ljSigma = [0.3]
    state.forcefield = ff
    
    atoms = []
    residues = []
    
    # Add ions (Na+ and Cl-) - all same type now
    ion_positions = [
        ([1.0, 1.0, 1.0], 1.0, 0),   # Na+ at position, charge, type
        ([2.0, 2.0, 2.0], -1.0, 0),  # Cl-
        ([1.0, 2.0, 1.5], 1.0, 0),   # Na+
        ([2.0, 1.0, 1.5], -1.0, 0),  # Cl-
    ]
    
    atom_idx = 0
    for pos, charge, atom_type in ion_positions:
        atom = MCAtom()
        atom.x, atom.y, atom.z = pos
        atom.charge = charge
        atom.type = atom_type
        atoms.append(atom)
        
        # Each ion is its own residue
        res = MCResidue()
        res.active = True
        res.fixed = True  # Ions are fixed
        res.atomStart = atom_idx
        res.atomCount = 1
        res.type = atom_type
        residues.append(res)
        atom_idx += 1
    
    # Add water molecules (simplified as single point)
    water_positions = [
        [1.5, 1.5, 1.5],
        [2.5, 2.5, 0.5],
        [0.5, 2.5, 1.5],
        [2.5, 0.5, 1.5],
    ]
    
    for pos in water_positions:
        # Oxygen
        atom = MCAtom()
        atom.x, atom.y, atom.z = pos
        atom.charge = -0.834  # TIP3P oxygen charge
        atom.type = 0  # Same type as others
        atoms.append(atom)
        
        # Water residue
        res = MCResidue()
        res.active = True
        res.fixed = True  # Water is fixed
        res.atomStart = atom_idx
        res.atomCount = 1
        res.type = 0
        residues.append(res)
        atom_idx += 1
    
    # Add ligands (simplified molecules with partial charges)
    ligand_configs = [
        # Each ligand: list of (position, charge, type)
        [  # Ligand 1 at center
            ([1.5, 1.5, 2.0], 0.2, 0),
            ([1.6, 1.5, 2.0], -0.2, 0),
        ],
        [  # Ligand 2
            ([0.8, 0.8, 2.2], 0.15, 0),
            ([0.9, 0.8, 2.2], -0.15, 0),
        ],
        [  # Ligand 3
            ([2.2, 2.2, 0.8], 0.1, 0),
            ([2.3, 2.2, 0.8], -0.1, 0),
        ],
    ]
    
    ligand_start_indices = []
    for ligand_atoms in ligand_configs:
        ligand_start_indices.append(atom_idx)
        
        # Add ligand atoms
        for pos, charge, atom_type in ligand_atoms:
            atom = MCAtom()
            atom.x, atom.y, atom.z = pos
            atom.charge = charge
            atom.type = atom_type
            atoms.append(atom)
        
        # Ligand residue
        res = MCResidue()
        res.active = True
        res.fixed = False  # Ligands are movable
        res.atomStart = atom_idx
        res.atomCount = len(ligand_atoms)
        res.type = 0  # Same type
        residues.append(res)
        atom_idx += len(ligand_atoms)
    
    state.atoms = atoms
    state.activeAtomCount = len(atoms)
    state.residues = residues
    state.activeResidueCount = len(residues)
    
    return state, ligand_start_indices


def calculate_openmm_pme_energy(state):
    """Calculate PME energy using OpenMM for comparison"""
    if not OPENMM_AVAILABLE:
        return None
    
    # Create OpenMM system
    system = System()
    
    # Add particles
    for atom in state.atoms:
        # Use hydrogen mass for all atoms (doesn't affect energy)
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
    nonbonded.setEwaldErrorTolerance(1e-5)
    
    # Add particles with charges (set epsilon=0 to get only electrostatic)
    for i, atom in enumerate(state.atoms):
        nonbonded.addParticle(
            atom.charge * elementary_charge,
            1.0 * nanometer,  # dummy sigma
            0.0 * kilojoule_per_mole  # epsilon = 0 for electrostatic only
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



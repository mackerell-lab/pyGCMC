"""
High-precision PME energy comparison between PyGCMC and OpenMM

This test ensures PME implementations agree to within 1% for realistic systems.
Critical for validating GCMC energy calculations.
"""

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


def create_test_system(n_ions=4, n_ligand_atoms=4):
    """Create a simple test system with ions and ligands"""
    
    # System parameters - use values that work well for PME
    box_size = 4.0  # nm, larger box for better PME convergence
    cutoff = 1.4    # nm, standard cutoff
    
    # Create state
    state = MCState()
    state.info.box = [box_size, box_size, box_size]
    state.info.cutoff = cutoff
    
    # Simple force field
    ff = MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.ljEps = [0.0]  # No LJ, pure electrostatics
    ff.ljSigma = [0.3]
    state.forcefield = ff
    
    atoms = []
    residues = []
    
    # Add ions in a regular pattern with better spacing
    # For 8 ions, use 2x2x2 grid
    grid_size = int(math.ceil(n_ions ** (1/3)))
    ion_spacing = box_size / (grid_size + 1)
    charges = [1.0, -1.0] * (n_ions // 2)  # Alternating charges
    
    atom_idx = 0
    for i in range(n_ions):
        atom = MCAtom()
        # Place ions in a 3D grid
        ix = i % grid_size
        iy = (i // grid_size) % grid_size
        iz = i // (grid_size * grid_size)
        atom.x = ion_spacing * (1 + ix)
        atom.y = ion_spacing * (1 + iy)
        atom.z = ion_spacing * (1 + iz)
        atom.charge = charges[i]
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
    
    # Add ligand atoms in the center
    ligand_start = atom_idx
    ligand_center = box_size / 2.0
    ligand_radius = 0.3  # nm, larger radius to avoid overlaps
    
    for i in range(n_ligand_atoms):
        atom = MCAtom()
        angle = 2.0 * math.pi * i / n_ligand_atoms
        atom.x = ligand_center + ligand_radius * math.cos(angle)
        atom.y = ligand_center + ligand_radius * math.sin(angle)
        atom.z = ligand_center
        # Small charges that sum to zero
        atom.charge = 0.1 if i % 2 == 0 else -0.1
        atom.type = 0
        atoms.append(atom)
    
    # Ligand residue
    res = MCResidue()
    res.active = True
    res.fixed = False
    res.atomStart = ligand_start
    res.atomCount = n_ligand_atoms
    res.type = 0
    residues.append(res)
    
    state.atoms = atoms
    state.activeAtomCount = len(atoms)
    state.residues = residues
    state.activeResidueCount = len(residues)
    
    return state


def calculate_openmm_pme_precise(state, alpha, mesh_size):
    """Calculate PME energy using OpenMM with precise parameters"""
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
    
    # Create PME force with specific parameters
    nonbonded = NonbondedForce()
    nonbonded.setNonbondedMethod(NonbondedForce.PME)
    nonbonded.setCutoffDistance(state.info.cutoff * nanometer)
    
    # Set PME parameters to match PyGCMC
    # OpenMM uses error tolerance, we need to convert from alpha
    # Smaller tolerance = more accurate
    nonbonded.setEwaldErrorTolerance(1e-6)
    
    # Try to set PME parameters directly if available
    try:
        # Set alpha (Ewald parameter) 
        nonbonded.setPMEParameters(alpha, mesh_size[0], mesh_size[1], mesh_size[2])
    except:
        # If direct setting not available, rely on error tolerance
        pass
    
    # Add particles with charges only
    for i, atom in enumerate(state.atoms):
        nonbonded.addParticle(
            atom.charge * elementary_charge,
            0.1 * nanometer,  # Small sigma to avoid numerical issues
            0.0 * kilojoule_per_mole  # epsilon = 0
        )
    
    system.addForce(nonbonded)
    
    # Create context with Reference platform for consistency
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
    
    # Get PME parameters that were actually used
    forces = system.getForces()
    for force in forces:
        if isinstance(force, NonbondedForce):
            print(f"  OpenMM cutoff: {force.getCutoffDistance()}")
            print(f"  OpenMM error tolerance: {force.getEwaldErrorTolerance()}")
            break
    
    return energy_kj_mol



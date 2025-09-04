"""
Helper functions for PME LJ-only tests

Contains shared setup functions for creating test systems with only LJ interactions
(no charges) to isolate LJ calculation testing.
"""

import math
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pygcmc
from pygcmc import MCState, MCAtom, MCResidue, MCForceField

# Only import OpenMM if available
try:
    from openmm import *
    from openmm.app import *
    from openmm.unit import *
    OPENMM_AVAILABLE = True
except ImportError:
    OPENMM_AVAILABLE = False


def create_lj_only_system(n_atoms=6):
    """Create a system with only LJ interactions (no charges)"""
    
    box_size = 3.0  # nm
    cutoff = 1.2    # nm
    
    state = MCState()
    state.info.box = [box_size, box_size, box_size]
    state.info.cutoff = cutoff
    
    # Force field with real LJ parameters
    ff = MCForceField()
    ff.numTotalTypes = 2  # Two types for testing mixing rules
    ff.numMovementTypes = 2
    # Need to provide full interaction matrix (2x2 = 4 values)
    # Using Lorentz-Berthelot mixing rules:
    # eps_ij = sqrt(eps_i * eps_j)
    # sigma_ij = (sigma_i + sigma_j) / 2
    eps0 = 0.996   # Argon-like
    eps1 = 1.230   # Methane-like
    sigma0 = 0.340
    sigma1 = 0.373
    
    # Full matrix: [00, 01, 10, 11]
    ff.ljEps = [
        eps0,                          # 0-0
        math.sqrt(eps0 * eps1),        # 0-1
        math.sqrt(eps0 * eps1),        # 1-0
        eps1                           # 1-1
    ]
    ff.ljSigma = [
        sigma0,                        # 0-0
        (sigma0 + sigma1) / 2,         # 0-1
        (sigma0 + sigma1) / 2,         # 1-0
        sigma1                         # 1-1
    ]
    state.forcefield = ff
    
    atoms = []
    residues = []
    
    # Create atoms in a regular pattern with safe distances (> 0.4 nm apart)
    positions = [
        ([1.0, 1.5, 1.5], 0),  # Type 0
        ([1.6, 1.5, 1.5], 0),  # Type 0 - increased spacing
        ([2.2, 1.5, 1.5], 0),  # Type 0 - increased spacing
        ([1.5, 0.9, 1.5], 1),  # Type 1 - increased spacing
        ([1.5, 2.1, 1.5], 1),  # Type 1 - increased spacing
        ([1.5, 1.5, 2.2], 1),  # Type 1 - increased spacing
    ]
    
    for i, (pos, atom_type) in enumerate(positions[:n_atoms]):
        atom = MCAtom()
        atom.x, atom.y, atom.z = pos
        atom.charge = 0.0  # NO CHARGE
        atom.type = atom_type
        atoms.append(atom)
        
        res = MCResidue()
        res.active = True
        res.fixed = True
        res.atomStart = i
        res.atomCount = 1
        res.type = atom_type
        residues.append(res)
    
    state.atoms = atoms
    state.activeAtomCount = len(atoms)
    state.residues = residues
    state.activeResidueCount = len(residues)
    
    return state


def create_lj_only_system_with_molecules(n_molecules=2, atoms_per_mol=3):
    """Create a system with molecules containing multiple atoms (for intramolecular testing)"""
    
    box_size = 3.0  # nm
    cutoff = 1.2    # nm
    
    state = MCState()
    state.info.box = [box_size, box_size, box_size]
    state.info.cutoff = cutoff
    
    # Force field with single LJ type for simplicity
    ff = MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.ljEps = [1.0]    # Simple value
    ff.ljSigma = [0.35]  # Simple value
    state.forcefield = ff
    
    atoms = []
    residues = []
    
    # Create molecules with atoms arranged in a line
    atom_idx = 0
    for mol in range(n_molecules):
        res = MCResidue()
        res.active = True
        res.fixed = True
        res.atomStart = atom_idx
        res.atomCount = atoms_per_mol
        res.type = 0
        
        # Place atoms in a line with 0.15 nm spacing (within LJ minimum)
        # Keep molecules well within box boundaries
        base_x = 0.5 + mol * 0.7  # Molecules 0.7 nm apart, starting at 0.5
        for i in range(atoms_per_mol):
            atom = MCAtom()
            atom.x = base_x + i * 0.15
            atom.y = 1.5
            atom.z = 1.5
            atom.charge = 0.0
            atom.type = 0
            atoms.append(atom)
            atom_idx += 1
        
        residues.append(res)
    
    state.atoms = atoms
    state.activeAtomCount = len(atoms)
    state.residues = residues
    state.activeResidueCount = len(residues)
    
    return state


def calculate_openmm_lj_energy(state, use_pme=False):
    """Calculate LJ-only energy using OpenMM"""
    if not OPENMM_AVAILABLE:
        return None
    
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
    
    # Create NonbondedForce
    nonbonded = NonbondedForce()
    if use_pme:
        nonbonded.setNonbondedMethod(NonbondedForce.PME)
    else:
        nonbonded.setNonbondedMethod(NonbondedForce.CutoffPeriodic)
    nonbonded.setCutoffDistance(state.info.cutoff * nanometer)
    
    # Add particles with LJ parameters only
    ff = state.forcefield
    for atom in state.atoms:
        nonbonded.addParticle(
            0.0 * elementary_charge,  # No charge
            ff.ljSigma[atom.type] * nanometer,
            ff.ljEps[atom.type] * kilojoule_per_mole
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
    return energy_state.getPotentialEnergy().value_in_unit(kilojoule_per_mole)
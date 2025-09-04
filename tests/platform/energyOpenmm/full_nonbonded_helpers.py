"""
Comprehensive test comparing full nonbonded interactions (PME + LJ) between PyGCMC and OpenMM

This test fills the gap in existing tests by including both electrostatic and LJ interactions
"""

import pytest
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pygcmc
from pygcmc import MCState, MCAtom, MCResidue, MCForceField
from pygcmc import initializePMEParameters, computeSystemEnergyPME, computeSystemEnergyPMEFixed, computeSystemEnergyCutoffFixed

# Only import OpenMM if available
try:
    from openmm import *
    from openmm.app import *
    from openmm.unit import *
    OPENMM_AVAILABLE = True
except ImportError:
    OPENMM_AVAILABLE = False


def create_test_system_with_lj():
    """Create a test system with both charges and LJ parameters"""
    
    box_size = 4.0  # nm
    cutoff = 1.2    # nm
    
    state = MCState()
    state.info.box = [box_size, box_size, box_size]
    state.info.cutoff = cutoff
    
    # Force field with real LJ parameters
    ff = MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    # Argon-like LJ parameters (epsilon in kJ/mol, sigma in nm)
    ff.ljEps = [0.996]    # ~1 kJ/mol
    ff.ljSigma = [0.340]  # 3.4 Angstrom
    state.forcefield = ff
    
    atoms = []
    residues = []
    
    # Create a system with various atom types
    # Place atoms closer together but still > 0.4 nm apart
    test_atoms = [
        # Position, charge, LJ type
        ([2.0, 2.0, 2.0], 0.5, 0),    # Positive charge at center
        ([2.5, 2.0, 2.0], -0.5, 0),   # Negative charge - 0.5 nm away
        ([2.0, 2.5, 2.0], 0.2, 0),    # Small positive - 0.5 nm away
        ([2.0, 2.0, 2.5], -0.2, 0),   # Small negative - 0.5 nm away
        ([2.5, 2.5, 2.0], 0.0, 0),    # Neutral (LJ only) - 0.707 nm from center
        ([2.0, 2.5, 2.5], 0.0, 0),    # Another neutral - 0.707 nm from center
    ]
    
    for i, (pos, charge, lj_type) in enumerate(test_atoms):
        atom = MCAtom()
        atom.x, atom.y, atom.z = pos
        atom.charge = charge
        atom.type = lj_type
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


def calculate_openmm_full_energy(state, alpha):
    """Calculate full nonbonded energy using OpenMM (PME + LJ)"""
    if not OPENMM_AVAILABLE:
        return None, None, None
    
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
    
    # Create NonbondedForce with PME for electrostatics
    nonbonded = NonbondedForce()
    nonbonded.setNonbondedMethod(NonbondedForce.PME)
    nonbonded.setCutoffDistance(state.info.cutoff * nanometer)
    nonbonded.setEwaldErrorTolerance(1e-6)
    
    # Add particles with charges and LJ parameters
    ff = state.forcefield
    for atom in state.atoms:
        charge = atom.charge * elementary_charge
        sigma = ff.ljSigma[atom.type] * nanometer
        epsilon = ff.ljEps[atom.type] * kilojoule_per_mole
        nonbonded.addParticle(charge, sigma, epsilon)
    
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
    
    # Get total energy
    energy_state = context.getState(getEnergy=True)
    total_energy = energy_state.getPotentialEnergy().value_in_unit(kilojoule_per_mole)
    
    # Get energy breakdown by zeroing charges or LJ
    # 1. LJ only (zero charges)
    for i in range(len(state.atoms)):
        nonbonded.setParticleParameters(
            i, 0.0 * elementary_charge,
            ff.ljSigma[state.atoms[i].type] * nanometer,
            ff.ljEps[state.atoms[i].type] * kilojoule_per_mole
        )
    nonbonded.updateParametersInContext(context)
    
    lj_state = context.getState(getEnergy=True)
    lj_energy = lj_state.getPotentialEnergy().value_in_unit(kilojoule_per_mole)
    
    # 2. Electrostatics only (restore charges, zero LJ)
    for i in range(len(state.atoms)):
        nonbonded.setParticleParameters(
            i, state.atoms[i].charge * elementary_charge,
            ff.ljSigma[state.atoms[i].type] * nanometer,
            0.0 * kilojoule_per_mole
        )
    nonbonded.updateParametersInContext(context)
    
    elec_state = context.getState(getEnergy=True)
    elec_energy = elec_state.getPotentialEnergy().value_in_unit(kilojoule_per_mole)
    
    return total_energy, elec_energy, lj_energy



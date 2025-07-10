"""
Deep analysis of PME electrostatic energy differences between PyGCMC and OpenMM

This script analyzes various sources of PME discrepancies:
1. Real space cutoff handling
2. Reciprocal space (k-space) calculation
3. Self-energy correction
4. PME parameter effects (alpha, mesh size, spline order)
"""

import pytest
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pygcmc
from pygcmc import MCState, MCAtom, MCResidue, MCForceField
from pygcmc import initializePMEParameters, computeSystemEnergyPME, computeSystemEnergyPMEComplete
from pygcmc import initializeEwaldParameters, computeSystemEnergyEwald

try:
    from openmm import *
    from openmm.app import *
    from openmm.unit import *
    OPENMM_AVAILABLE = True
except ImportError:
    OPENMM_AVAILABLE = False


def create_simple_charged_system():
    """Create a simple system with two opposite charges"""
    state = MCState()
    state.info.box = [3.0, 3.0, 3.0]
    state.info.cutoff = 1.2
    
    # Simple force field
    ff = MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.ljEps = [0.0]    # No LJ
    ff.ljSigma = [0.1]
    state.forcefield = ff
    
    # Two atoms with opposite charges
    atoms = []
    
    # Positive charge
    atom1 = MCAtom()
    atom1.x, atom1.y, atom1.z = 1.5, 1.5, 1.2
    atom1.charge = 1.0
    atom1.type = 0
    atoms.append(atom1)
    
    # Negative charge
    atom2 = MCAtom()
    atom2.x, atom2.y, atom2.z = 1.5, 1.5, 1.8
    atom2.charge = -1.0
    atom2.type = 0
    atoms.append(atom2)
    
    # Create residues
    residues = []
    for i in range(2):
        res = MCResidue()
        res.active = True
        res.fixed = True
        res.atomStart = i
        res.atomCount = 1
        res.type = 0
        residues.append(res)
    
    state.atoms = atoms
    state.activeAtomCount = 2
    state.residues = residues
    state.activeResidueCount = 2
    
    return state


def calculate_openmm_pme_details(state, alpha, mesh_size):
    """Get detailed PME information from OpenMM"""
    if not OPENMM_AVAILABLE:
        return {}
    
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
    
    # Try to set PME parameters explicitly
    nonbonded.setEwaldErrorTolerance(1e-9)  # Very high precision
    
    for atom in state.atoms:
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
    
    # Get PME parameters that OpenMM actually uses
    actual_alpha = nonbonded.getPMEParametersInContext(context)[0]
    actual_mesh = nonbonded.getPMEParametersInContext(context)[1:4]
    
    state_obj = context.getState(getEnergy=True)
    energy = state_obj.getPotentialEnergy().value_in_unit(kilojoule_per_mole)
    
    return {
        'energy': energy,
        'alpha': actual_alpha,
        'mesh': actual_mesh
    }


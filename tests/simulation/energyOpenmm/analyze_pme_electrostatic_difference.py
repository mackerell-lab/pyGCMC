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


def analyze_pme_components(state):
    """Analyze PME energy components step by step"""
    
    print("\n" + "="*70)
    print("PME Component Analysis")
    print("="*70)
    
    # Test different alpha values
    alphas = [3.0, 4.0, 5.0, 5.6, 6.0, 7.0]
    mesh_sizes = [[24, 24, 24], [32, 32, 32], [48, 48, 48]]
    
    print("\n1. Alpha Parameter Sensitivity:")
    print("-"*50)
    print(f"{'Alpha':>6} | {'PyGCMC Energy':>15} | {'OpenMM Energy':>15} | {'Difference':>12} | {'Rel Error':>10}")
    print("-"*50)
    
    for alpha in alphas:
        mesh_size = [32, 32, 32]
        spline_order = 4
        
        # PyGCMC calculation
        pygcmc.setPMEParameters(alpha, mesh_size, spline_order)
        initializePMEParameters(state.info.cutoff, state.info.box, alpha, mesh_size, spline_order)
        elec, vdw, total = computeSystemEnergyPMEComplete(state)
        
        # OpenMM calculation
        openmm_result = calculate_openmm_pme_details(state, alpha, mesh_size)
        openmm_energy = openmm_result['energy']
        
        diff = elec - openmm_energy
        rel_error = abs(diff) / abs(openmm_energy) * 100 if openmm_energy != 0 else 0
        
        print(f"{alpha:6.1f} | {elec:15.6f} | {openmm_energy:15.6f} | {diff:12.6f} | {rel_error:9.3f}%")
    
    print("\n2. Mesh Size Sensitivity (alpha=5.6):")
    print("-"*50)
    print(f"{'Mesh':>12} | {'PyGCMC Energy':>15} | {'OpenMM Energy':>15} | {'Difference':>12} | {'Rel Error':>10}")
    print("-"*50)
    
    alpha = 5.6
    for mesh_size in mesh_sizes:
        spline_order = 4
        
        # PyGCMC calculation
        pygcmc.setPMEParameters(alpha, mesh_size, spline_order)
        initializePMEParameters(state.info.cutoff, state.info.box, alpha, mesh_size, spline_order)
        elec, vdw, total = computeSystemEnergyPMEComplete(state)
        
        # OpenMM calculation
        openmm_result = calculate_openmm_pme_details(state, alpha, mesh_size)
        openmm_energy = openmm_result['energy']
        
        diff = elec - openmm_energy
        rel_error = abs(diff) / abs(openmm_energy) * 100 if openmm_energy != 0 else 0
        
        mesh_str = f"{mesh_size[0]}x{mesh_size[1]}x{mesh_size[2]}"
        print(f"{mesh_str:>12} | {elec:15.6f} | {openmm_energy:15.6f} | {diff:12.6f} | {rel_error:9.3f}%")
    
    # Compare with Ewald summation
    print("\n3. PME vs Ewald Comparison:")
    print("-"*50)
    print("(Skipping Ewald comparison due to API differences)")


def analyze_error_sources():
    """Analyze specific sources of PME errors"""
    
    print("\n" + "="*70)
    print("PME Error Source Analysis")
    print("="*70)
    
    # Test with different charge separations
    separations = [0.3, 0.5, 0.6, 0.8, 1.0, 1.2, 1.5]
    
    print("\n4. Distance-Dependent Error Analysis:")
    print("-"*70)
    print(f"{'Distance':>8} | {'PyGCMC':>12} | {'OpenMM':>12} | {'Diff':>10} | {'Rel Error':>10} | {'Coulomb':>12}")
    print("-"*70)
    
    for sep in separations:
        state = MCState()
        state.info.box = [3.0, 3.0, 3.0]
        state.info.cutoff = 1.2
        
        ff = MCForceField()
        ff.numTotalTypes = 1
        ff.numMovementTypes = 1
        ff.ljEps = [0.0]
        ff.ljSigma = [0.1]
        state.forcefield = ff
        
        # Two charges separated by distance 'sep'
        atoms = []
        
        atom1 = MCAtom()
        atom1.x, atom1.y, atom1.z = 1.5, 1.5, 1.5 - sep/2
        atom1.charge = 1.0
        atom1.type = 0
        atoms.append(atom1)
        
        atom2 = MCAtom()
        atom2.x, atom2.y, atom2.z = 1.5, 1.5, 1.5 + sep/2
        atom2.charge = -1.0
        atom2.type = 0
        atoms.append(atom2)
        
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
        
        # Calculate energies
        alpha = 5.6
        mesh_size = [32, 32, 32]
        spline_order = 4
        
        pygcmc.setPMEParameters(alpha, mesh_size, spline_order)
        initializePMEParameters(state.info.cutoff, state.info.box, alpha, mesh_size, spline_order)
        elec, vdw, total = computeSystemEnergyPMEComplete(state)
        
        openmm_result = calculate_openmm_pme_details(state, alpha, mesh_size)
        openmm_energy = openmm_result['energy']
        
        # Simple Coulomb energy for reference
        ke = 138.935456  # Coulomb constant in kJ*nm/(mol*e^2)
        coulomb = -ke / sep
        
        diff = elec - openmm_energy
        rel_error = abs(diff) / abs(openmm_energy) * 100 if openmm_energy != 0 else 0
        
        print(f"{sep:8.2f} | {elec:12.4f} | {openmm_energy:12.4f} | {diff:10.4f} | {rel_error:9.3f}% | {coulomb:12.4f}")


def test_pme_self_energy():
    """Test self-energy calculation differences"""
    
    print("\n5. Self-Energy Calculation:")
    print("-"*50)
    
    # Single charge system
    state = MCState()
    state.info.box = [3.0, 3.0, 3.0]
    state.info.cutoff = 1.2
    
    ff = MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.ljEps = [0.0]
    ff.ljSigma = [0.1]
    state.forcefield = ff
    
    # Single charge
    atom = MCAtom()
    atom.x, atom.y, atom.z = 1.5, 1.5, 1.5
    atom.charge = 2.0  # Use charge 2 to make self-energy more significant
    atom.type = 0
    
    res = MCResidue()
    res.active = True
    res.fixed = True
    res.atomStart = 0
    res.atomCount = 1
    res.type = 0
    
    state.atoms = [atom]
    state.activeAtomCount = 1
    state.residues = [res]
    state.activeResidueCount = 1
    
    # Calculate with different alphas
    alphas = [4.0, 5.0, 5.6, 6.0, 7.0]
    
    print(f"{'Alpha':>6} | {'Self Energy (Theory)':>20} | {'PyGCMC Total':>15} | {'OpenMM Total':>15}")
    print("-"*70)
    
    for alpha in alphas:
        mesh_size = [32, 32, 32]
        spline_order = 4
        
        # Theoretical self-energy: -alpha * q^2 / sqrt(pi)
        ke = 138.935456  # kJ*nm/(mol*e^2)
        self_energy_theory = -alpha * atom.charge**2 * ke / math.sqrt(math.pi)
        
        # PyGCMC
        pygcmc.setPMEParameters(alpha, mesh_size, spline_order)
        initializePMEParameters(state.info.cutoff, state.info.box, alpha, mesh_size, spline_order)
        elec, vdw, total = computeSystemEnergyPMEComplete(state)
        
        # OpenMM
        openmm_result = calculate_openmm_pme_details(state, alpha, mesh_size)
        openmm_energy = openmm_result['energy']
        
        print(f"{alpha:6.1f} | {self_energy_theory:20.6f} | {elec:15.6f} | {openmm_energy:15.6f}")


if __name__ == "__main__":
    print("PyGCMC vs OpenMM PME Electrostatic Difference Analysis")
    
    # Basic component analysis
    state = create_simple_charged_system()
    analyze_pme_components(state)
    
    # Error source analysis
    analyze_error_sources()
    
    # Self-energy test
    test_pme_self_energy()
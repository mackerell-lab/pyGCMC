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


@pytest.mark.skipif(not OPENMM_AVAILABLE, reason="OpenMM not available")
def test_full_nonbonded_pme_lj():
    """Test full nonbonded energy (PME + LJ) comparison with OpenMM"""
    
    # Create test system
    state = create_test_system_with_lj()
    
    print(f"\nTest system with LJ and charges:")
    print(f"  Atoms: {state.activeAtomCount}")
    print(f"  Box: {state.info.box[0]} nm")
    print(f"  Cutoff: {state.info.cutoff} nm")
    print(f"  LJ parameters: eps={state.forcefield.ljEps[0]} kJ/mol, sigma={state.forcefield.ljSigma[0]} nm")
    
    # PME parameters
    alpha = 3.2
    mesh_size = [32, 32, 32]
    spline_order = 4
    
    # Initialize PyGCMC PME
    pygcmc.setPMEParameters(alpha, mesh_size, spline_order)
    initializePMEParameters(state.info.cutoff, state.info.box, alpha, mesh_size, spline_order)
    
    # Calculate PyGCMC energy using fixed function
    elec_total, vdw_total, pme_dict = computeSystemEnergyPMEFixed(state)
    pygcmc_elec = pme_dict['total']
    pygcmc_lj = vdw_total
    
    pygcmc_total = pygcmc_elec + pygcmc_lj
    
    # Calculate OpenMM energy
    openmm_total, openmm_elec, openmm_lj = calculate_openmm_full_energy(state, alpha)
    
    if openmm_total is not None:
        print(f"\nEnergy comparison:")
        print(f"  PyGCMC:")
        print(f"    Electrostatic: {pygcmc_elec:.6f} kJ/mol")
        print(f"    LJ: {pygcmc_lj:.6f} kJ/mol")
        print(f"    Total: {pygcmc_total:.6f} kJ/mol")
        print(f"  OpenMM:")
        print(f"    Electrostatic: {openmm_elec:.6f} kJ/mol")
        print(f"    LJ: {openmm_lj:.6f} kJ/mol")
        print(f"    Total: {openmm_total:.6f} kJ/mol")
        
        # Calculate differences
        elec_diff = abs(pygcmc_elec - openmm_elec)
        elec_rel = elec_diff / abs(openmm_elec) if openmm_elec != 0 else 0
        
        lj_diff = abs(pygcmc_lj - openmm_lj)
        lj_rel = lj_diff / abs(openmm_lj) if openmm_lj != 0 else 0
        
        total_diff = abs(pygcmc_total - openmm_total)
        total_rel = total_diff / abs(openmm_total) if openmm_total != 0 else 0
        
        print(f"\nRelative differences:")
        print(f"  Electrostatic: {elec_rel*100:.3f}%")
        print(f"  LJ: {lj_rel*100:.3f}%")
        print(f"  Total: {total_rel*100:.3f}%")
        
        # Assertions
        # PME implementations can differ slightly between PyGCMC and OpenMM
        assert elec_rel < 0.05, f"Electrostatic energies differ by {elec_rel*100:.3f}% (> 5%)"
        assert lj_rel < 0.01, f"LJ energies differ by {lj_rel*100:.3f}% (> 1%)"
        assert total_rel < 0.05, f"Total energies differ by {total_rel*100:.3f}% (> 5%)"
        
        print("\n✓ Full nonbonded energy calculation matches OpenMM!")
    else:
        print("\nOpenMM not available for comparison")


@pytest.mark.skipif(not OPENMM_AVAILABLE, reason="OpenMM not available")
def test_lj_only_comparison():
    """Test LJ-only energy comparison with OpenMM"""
    
    # Create system with no charges
    state = create_test_system_with_lj()
    
    # Zero all charges
    for atom in state.atoms:
        atom.charge = 0.0
    
    print(f"\nLJ-only test system:")
    print(f"  Atoms: {state.activeAtomCount}")
    print(f"  All charges set to 0")
    
    # No PME needed for LJ-only
    
    # Calculate PyGCMC LJ energy using fixed function
    elec, vdw, total = computeSystemEnergyCutoffFixed(state)
    pygcmc_lj = vdw
    
    # Calculate OpenMM LJ energy
    if OPENMM_AVAILABLE:
        # Use simple cutoff for LJ-only
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
        nonbonded.setNonbondedMethod(NonbondedForce.CutoffPeriodic)
        nonbonded.setCutoffDistance(state.info.cutoff * nanometer)
        
        ff = state.forcefield
        for atom in state.atoms:
            nonbonded.addParticle(
                0.0 * elementary_charge,
                ff.ljSigma[atom.type] * nanometer,
                ff.ljEps[atom.type] * kilojoule_per_mole
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
        openmm_lj = energy_state.getPotentialEnergy().value_in_unit(kilojoule_per_mole)
        
        print(f"\nLJ energy comparison:")
        print(f"  PyGCMC: {pygcmc_lj:.6f} kJ/mol")
        print(f"  OpenMM: {openmm_lj:.6f} kJ/mol")
        
        diff = abs(pygcmc_lj - openmm_lj)
        rel_diff = diff / abs(openmm_lj) if openmm_lj != 0 else 0
        print(f"  Relative difference: {rel_diff*100:.3f}%")
        
        assert rel_diff < 0.01, f"LJ energies differ by {rel_diff*100:.3f}% (> 1%)"
        
        print("\n✓ LJ-only energy calculation matches OpenMM!")


if __name__ == "__main__":
    test_full_nonbonded_pme_lj()
    test_lj_only_comparison()
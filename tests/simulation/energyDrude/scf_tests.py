# tests/simulation/energyDrude/scf_tests.py
"""SCF convergence and optimization tests for Drude oscillators."""

import pytest
import math
import pygcmc
from pygcmc import MCState, DrudeForce, DrudeSCFParams
from .helpers import (
    create_simple_drude_system,
    create_water_box,
    ONE_4PI_EPS0
)


def test_drude_scf_convergence():
    """Test SCF convergence for single Drude particle"""
    state = create_simple_drude_system()
    
    # Create Drude force
    drude_force = DrudeForce()
    
    # Add Drude particle
    charge = -1.0
    polarizability = 0.001
    drude_force.addParticle(
        drudeIndex=1, parentIndex=0,
        aniso1Index=-1, aniso2Index=-1,
        aniso3Index=-1, aniso4Index=-1,
        charge=charge, polarizability=polarizability,
        aniso12=1.0, aniso34=1.0
    )
    
    # Set very tight SCF parameters
    scf_params = DrudeSCFParams()
    scf_params.tolerance = 1e-8  # Very tight tolerance
    scf_params.maxIterations = 100
    scf_params.dampingFactor = 0.5
    scf_params.forceCutoff = 10.0
    drude_force.setSCFParameters(scf_params)
    
    # Place Drude far from equilibrium
    state.atoms[1].x = 6.0  # 1 nm away from parent
    
    # Calculate energy with SCF
    energy = drude_force.calculateEnergySCF(state)
    
    # After SCF, Drude should be very close to parent
    dx = state.atoms[1].x - state.atoms[0].x
    dy = state.atoms[1].y - state.atoms[0].y
    dz = state.atoms[1].z - state.atoms[0].z
    distance = math.sqrt(dx*dx + dy*dy + dz*dz)
    
    # Should converge to near zero separation
    assert distance < 1e-6


def test_drude_scf_tolerance():
    """Test effect of different SCF tolerances"""
    tolerances = [1e-4, 1e-6, 1e-8]
    energies = []
    
    for tol in tolerances:
        state = create_simple_drude_system()
        
        # Add external charge to create field
        external_charge = pygcmc.MCAtom()
        external_charge.x = 6.0
        external_charge.y = 5.0
        external_charge.z = 5.0
        external_charge.charge = 2.0
        external_charge.type = 2
        # external_charge.mass = 1.0  # Mass not needed for SCF
        state.atoms.append(external_charge)
        state.activeAtomCount = 3
        
        # Create Drude force
        drude_force = DrudeForce()
        
        charge = -1.0
        polarizability = 0.001
        drude_force.addParticle(
            drudeIndex=1, parentIndex=0,
            aniso1Index=-1, aniso2Index=-1,
            aniso3Index=-1, aniso4Index=-1,
            charge=charge, polarizability=polarizability,
            aniso12=1.0, aniso34=1.0
        )
        
        # Set SCF parameters
        scf_params = DrudeSCFParams()
        scf_params.tolerance = tol
        scf_params.maxIterations = 100
        drude_force.setSCFParameters(scf_params)
        
        # Place Drude away from equilibrium
        state.atoms[1].x = 5.5
        
        # Calculate energy
        energy = drude_force.calculateEnergySCF(state)
        energies.append(energy)
    
    # Energies should converge as tolerance gets tighter
    assert abs(energies[1] - energies[2]) < abs(energies[0] - energies[1])
    # Tightest tolerance should give most accurate energy
    assert abs(energies[1] - energies[2]) < 1e-6


def test_drude_scf_damping():
    """Test SCF damping for stability with strong forces"""
    state = create_simple_drude_system()
    
    # Add very strong external charge
    external_charge = pygcmc.MCAtom()
    external_charge.x = 5.2  # Very close
    external_charge.y = 5.0
    external_charge.z = 5.0
    external_charge.charge = 10.0  # Strong charge
    external_charge.type = 2
    # external_charge.mass = 1.0  # Mass not needed for SCF
    state.atoms.append(external_charge)
    state.activeAtomCount = 3
    
    # Create Drude force
    drude_force = DrudeForce()
    
    charge = -1.0
    polarizability = 0.001
    drude_force.addParticle(
        drudeIndex=1, parentIndex=0,
        aniso1Index=-1, aniso2Index=-1,
        aniso3Index=-1, aniso4Index=-1,
        charge=charge, polarizability=polarizability,
        aniso12=1.0, aniso34=1.0
    )
    
    # Test with different damping factors
    damping_factors = [0.3, 0.5, 0.7]
    converged_count = []
    
    for damping in damping_factors:
        # Reset Drude position
        state.atoms[1].x = 5.5
        state.atoms[1].y = 5.0
        state.atoms[1].z = 5.0
        
        scf_params = DrudeSCFParams()
        scf_params.tolerance = 1e-6
        scf_params.maxIterations = 50
        scf_params.dampingFactor = damping
        scf_params.forceCutoff = 10.0
        drude_force.setSCFParameters(scf_params)
        
        # Try to calculate energy
        try:
            energy = drude_force.calculateEnergySCF(state)
            converged_count.append(1)
        except:
            converged_count.append(0)
    
    # Lower damping should be more stable with strong forces
    assert sum(converged_count) >= 1  # At least one should converge


def test_drude_scf_multiple_particles():
    """Test SCF convergence with multiple Drude particles"""
    state = MCState()
    
    # Set box size
    state.info.box = [10.0, 10.0, 10.0]
    state.info.setTemperature(300.0)
    state.info.cutoff = 5.0
    state.info.pbc = True
    
    # Create chain of atoms with Drude particles
    atoms = []
    num_atoms = 4
    
    for i in range(num_atoms):
        # Parent atom
        atom = pygcmc.MCAtom()
        atom.x = 3.0 + i * 1.5
        atom.y = 5.0
        atom.z = 5.0
        atom.charge = 1.0
        atom.type = 0
        # atom.mass = 12.0  # Mass not needed for SCF
        atoms.append(atom)
        
        # Drude particle
        drude = pygcmc.MCAtom()
        drude.x = atom.x + 0.1 * (i % 2 * 2 - 1)  # Alternate displacement
        drude.y = atom.y + 0.1 * ((i+1) % 2 * 2 - 1)
        drude.z = atom.z
        drude.charge = -1.0
        drude.type = 1
        # drude.mass = 0.4  # Mass not needed for SCF
        atoms.append(drude)
    
    state.atoms = atoms
    state.activeAtomCount = len(atoms)
    
    # Create residues
    residues = []
    for i in range(num_atoms):
        res = pygcmc.MCResidue()
        res.atomStart = i * 2
        res.atomCount = 2
        res.active = True
        res.fixed = False
        res.type = 0
        residues.append(res)
    
    state.residues = residues
    state.activeResidueCount = num_atoms
    
    # Set force field
    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 2
    ff.numMovementTypes = 2
    ff.ljSigma = [0.3] * 4
    ff.ljEps = [0.0] * 4
    state.forcefield = ff
    
    # Create Drude force
    drude_force = DrudeForce()
    
    # Add all Drude particles
    charge = -1.0
    polarizability = 0.001
    
    for i in range(num_atoms):
        drude_force.addParticle(
            drudeIndex=i*2+1, parentIndex=i*2,
            aniso1Index=-1, aniso2Index=-1,
            aniso3Index=-1, aniso4Index=-1,
            charge=charge, polarizability=polarizability,
            aniso12=1.0, aniso34=1.0
        )
    
    # Set SCF parameters
    scf_params = DrudeSCFParams()
    scf_params.tolerance = 1e-6
    scf_params.maxIterations = 100
    drude_force.setSCFParameters(scf_params)
    
    # Calculate energy
    energy_initial = drude_force.calculateEnergySCF(state)
    
    # All Drude particles should find equilibrium
    total_displacement = 0.0
    for i in range(num_atoms):
        parent_idx = i * 2
        drude_idx = i * 2 + 1
        dx = state.atoms[drude_idx].x - state.atoms[parent_idx].x
        dy = state.atoms[drude_idx].y - state.atoms[parent_idx].y
        dz = state.atoms[drude_idx].z - state.atoms[parent_idx].z
        total_displacement += math.sqrt(dx*dx + dy*dy + dz*dz)
    
    # Average displacement should be small but non-zero due to interactions
    avg_displacement = total_displacement / num_atoms
    assert avg_displacement > 1e-4  # Non-zero due to interactions
    assert avg_displacement < 0.1   # But still small


def test_drude_scf_with_external_field():
    """Test SCF optimization in presence of external electric field"""
    state = create_simple_drude_system()
    
    # Create uniform external field using charges
    # Place positive charges on one side, negative on other
    for i in range(2):
        for sign, x in [(1, 2.0), (-1, 8.0)]:
            charge_atom = pygcmc.MCAtom()
            charge_atom.x = x
            charge_atom.y = 4.0 + i * 2.0
            charge_atom.z = 5.0
            charge_atom.charge = sign * 5.0
            charge_atom.type = 2
            # charge_atom.mass = 1.0  # Mass not needed for SCF
            state.atoms.append(charge_atom)
    
    state.activeAtomCount = len(state.atoms)
    
    # Create Drude force
    drude_force = DrudeForce()
    
    charge = -1.0
    polarizability = 0.001
    drude_force.addParticle(
        drudeIndex=1, parentIndex=0,
        aniso1Index=-1, aniso2Index=-1,
        aniso3Index=-1, aniso4Index=-1,
        charge=charge, polarizability=polarizability,
        aniso12=1.0, aniso34=1.0
    )
    
    # Calculate energy
    energy = drude_force.calculateEnergySCF(state)
    
    # Drude should be displaced in field direction (positive x)
    dx = state.atoms[1].x - state.atoms[0].x
    assert dx > 0.01  # Should be displaced in positive x direction
    
    # Calculate expected displacement from field
    # Simplified: displacement ≈ α * E / k
    # where E is electric field, k is spring constant
    k = ONE_4PI_EPS0 * charge * charge / polarizability
    
    # Energy should be non-zero due to field interaction
    assert energy > 0.0


def test_drude_scf_iteration_limit():
    """Test behavior when SCF hits iteration limit"""
    state = create_simple_drude_system()
    
    # Create impossible to converge scenario
    # Add oscillating external charges
    for i in range(10):
        charge_atom = pygcmc.MCAtom()
        charge_atom.x = 5.0 + 0.1 * math.cos(i * math.pi / 5)
        charge_atom.y = 5.0 + 0.1 * math.sin(i * math.pi / 5)
        charge_atom.z = 5.0
        charge_atom.charge = (-1)**i * 2.0
        charge_atom.type = 2
        # charge_atom.mass = 1.0  # Mass not needed for SCF
        state.atoms.append(charge_atom)
    
    state.activeAtomCount = len(state.atoms)
    
    # Create Drude force
    drude_force = DrudeForce()
    
    charge = -1.0
    polarizability = 0.001
    drude_force.addParticle(
        drudeIndex=1, parentIndex=0,
        aniso1Index=-1, aniso2Index=-1,
        aniso3Index=-1, aniso4Index=-1,
        charge=charge, polarizability=polarizability,
        aniso12=1.0, aniso34=1.0
    )
    
    # Set very few iterations
    scf_params = DrudeSCFParams()
    scf_params.tolerance = 1e-12  # Impossible tolerance
    scf_params.maxIterations = 5   # Very few iterations
    drude_force.setSCFParameters(scf_params)
    
    # Should still return an energy even if not converged
    energy = drude_force.calculateEnergySCF(state)
    assert energy is not None
    assert energy >= 0.0  # Energy should be non-negative
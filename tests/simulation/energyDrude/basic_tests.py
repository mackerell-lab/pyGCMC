# tests/simulation/energyDrude/basic_tests.py
"""Basic Drude oscillator tests."""

import pytest
import math
import pygcmc
from pygcmc import MCState, DrudeForce
from .helpers import (
    create_simple_drude_system,
    ONE_4PI_EPS0,
    calculate_drude_polarizability
)


def test_drude_initialization():
    """Test basic Drude force initialization"""
    drude_force = DrudeForce()
    
    # Test empty initialization
    assert drude_force.getNumParticles() == 0
    assert drude_force.getNumScreenedPairs() == 0
    
    # Add a Drude particle
    charge = -1.0
    polarizability = 0.001  # nm^3
    drude_idx = drude_force.addParticle(
        drudeIndex=1,
        parentIndex=0,
        aniso1Index=-1,
        aniso2Index=-1,
        aniso3Index=-1,
        aniso4Index=-1,
        charge=charge,
        polarizability=polarizability,
        aniso12=1.0,
        aniso34=1.0
    )
    
    assert drude_idx == 0
    assert drude_force.getNumParticles() == 1


def test_drude_harmonic_energy():
    """Test harmonic restraint energy calculation"""
    state = create_simple_drude_system()
    
    # Create Drude force
    drude_force = DrudeForce()
    
    # Add Drude particle with known parameters
    charge = -1.0
    polarizability = 0.001  # nm^3
    drude_force.addParticle(
        drudeIndex=1,  # Drude particle index
        parentIndex=0,  # Parent atom index
        aniso1Index=-1,
        aniso2Index=-1,
        aniso3Index=-1,
        aniso4Index=-1,
        charge=charge,
        polarizability=polarizability,
        aniso12=1.0,
        aniso34=1.0
    )
    
    # For an isolated Drude oscillator with no external field,
    # SCF should bring Drude to parent position, giving zero energy
    energy = drude_force.calculateEnergySCF(state)
    
    # After SCF, Drude should be at parent position
    dx = state.atoms[1].x - state.atoms[0].x
    dy = state.atoms[1].y - state.atoms[0].y
    dz = state.atoms[1].z - state.atoms[0].z
    distance = math.sqrt(dx*dx + dy*dy + dz*dz)
    
    # Check that Drude moved to parent position
    assert distance < 1e-6
    # And energy should be essentially zero
    assert abs(energy) < 1e-6
    
    # Test harmonic energy without SCF
    # Move Drude away from parent
    state.atoms[1].x = 5.1  # 0.1 nm displacement
    
    # Create new DrudeForce with no SCF iterations
    drude_force2 = DrudeForce()
    scf_params = pygcmc.DrudeSCFParams()
    scf_params.maxIterations = 0  # No optimization
    drude_force2.setSCFParameters(scf_params)
    
    drude_force2.addParticle(
        drudeIndex=1, parentIndex=0,
        aniso1Index=-1, aniso2Index=-1,
        aniso3Index=-1, aniso4Index=-1,
        charge=charge, polarizability=polarizability,
        aniso12=1.0, aniso34=1.0
    )
    
    # Calculate energy without optimization
    energy_displaced = drude_force2.calculateEnergySCF(state)
    
    # Expected harmonic energy
    k = ONE_4PI_EPS0 * charge * charge / polarizability
    expected = 0.5 * k * 0.1 * 0.1
    
    # Check harmonic energy is correct
    assert abs(energy_displaced - expected) < 1.0  # Within 1 kJ/mol


def test_drude_equilibrium_position():
    """Test that Drude particle finds equilibrium position"""
    state = create_simple_drude_system()
    
    # Apply external field by adding charge to parent
    state.atoms[0].charge = 2.0  # Parent has +2 charge
    
    # Create Drude force
    drude_force = DrudeForce()
    
    # Add Drude particle
    charge = -1.0
    polarizability = 0.001  # nm^3
    drude_force.addParticle(
        drudeIndex=1,
        parentIndex=0,
        aniso1Index=-1,
        aniso2Index=-1,
        aniso3Index=-1,
        aniso4Index=-1,
        charge=charge,
        polarizability=polarizability,
        aniso12=1.0,
        aniso34=1.0
    )
    
    # Set initial Drude position away from equilibrium
    state.atoms[1].x = 5.5  # 0.5 nm from parent
    state.atoms[1].y = 5.0
    state.atoms[1].z = 5.0
    
    # Calculate energy with SCF
    energy_initial = drude_force.calculateEnergySCF(state)
    
    # After SCF, Drude should be at equilibrium
    # At equilibrium, the force from the harmonic restraint balances external forces
    # For an isolated system, equilibrium is at parent position
    dx = state.atoms[1].x - state.atoms[0].x
    dy = state.atoms[1].y - state.atoms[0].y
    dz = state.atoms[1].z - state.atoms[0].z
    displacement = math.sqrt(dx*dx + dy*dy + dz*dz)
    
    # Should be very close to parent (within SCF tolerance)
    assert displacement < 1e-5


def test_drude_anisotropic_polarizability():
    """Test anisotropic polarizability"""
    state = MCState()
    
    # Set box size
    state.info.box = [10.0, 10.0, 10.0]
    state.info.setTemperature(300.0)
    state.info.cutoff = 5.0
    state.info.pbc = True
    
    # Create atoms for anisotropic system
    # Parent atom
    parent = pygcmc.MCAtom()
    parent.x = 5.0
    parent.y = 5.0
    parent.z = 5.0
    parent.charge = 1.0
    parent.type = 0
    parent.mass = 12.0
    
    # Drude particle
    drude = pygcmc.MCAtom()
    drude.x = 5.1
    drude.y = 5.0
    drude.z = 5.0
    drude.charge = -1.0
    drude.type = 1
    drude.mass = 0.4
    
    # Anisotropy axis atoms
    axis1 = pygcmc.MCAtom()
    axis1.x = 6.0
    axis1.y = 5.0
    axis1.z = 5.0
    axis1.charge = 0.0
    axis1.type = 2
    axis1.mass = 1.0
    
    axis2 = pygcmc.MCAtom()
    axis2.x = 4.0
    axis2.y = 5.0
    axis2.z = 5.0
    axis2.charge = 0.0
    axis2.type = 2
    axis2.mass = 1.0
    
    state.atoms = [parent, drude, axis1, axis2]
    state.activeAtomCount = 4
    
    # Create residue
    res = pygcmc.MCResidue()
    res.atomStart = 0
    res.atomCount = 4
    res.active = True
    res.fixed = False
    res.type = 0
    
    state.residues = [res]
    state.activeResidueCount = 1
    
    # Set force field
    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 3
    ff.numMovementTypes = 3
    ff.ljSigma = [0.3] * 9
    ff.ljEps = [0.0] * 9
    state.forcefield = ff
    
    # Create Drude force with anisotropic polarizability
    drude_force = DrudeForce()
    
    # Add anisotropic Drude particle
    # aniso12 < 1 means less polarizable along axis 1-2
    charge = -1.0
    polarizability = 0.001
    aniso12 = 0.5  # Half polarizability along x-axis
    drude_force.addParticle(
        drudeIndex=1,
        parentIndex=0,
        aniso1Index=2,  # axis1
        aniso2Index=3,  # axis2
        aniso3Index=-1,
        aniso4Index=-1,
        charge=charge,
        polarizability=polarizability,
        aniso12=aniso12,
        aniso34=1.0
    )
    
    # Set Drude displaced along x-axis
    state.atoms[1].x = 5.1
    
    # Calculate energy
    energy = drude_force.calculateEnergySCF(state)
    
    # Energy should be higher due to reduced polarizability along x
    # E = 0.5 * k_aniso * r^2 where k_aniso > k_iso
    k_iso = ONE_4PI_EPS0 * charge * charge / polarizability
    k_x = ONE_4PI_EPS0 * charge * charge / (polarizability * aniso12)
    
    dx = 0.1  # displacement along x
    expected_energy = 0.5 * k_x * dx * dx
    
    # Should be approximately equal (some difference due to SCF optimization)
    assert abs(energy - expected_energy) / expected_energy < 0.1


def test_drude_multiple_particles():
    """Test system with multiple Drude particles"""
    state = MCState()
    
    # Set box size
    state.info.box = [10.0, 10.0, 10.0]
    state.info.setTemperature(300.0)
    state.info.cutoff = 5.0
    state.info.pbc = True
    
    # Create two atoms with their Drude particles
    atoms = []
    
    # First atom and Drude
    atom1 = pygcmc.MCAtom()
    atom1.x = 3.0
    atom1.y = 5.0
    atom1.z = 5.0
    atom1.charge = 1.0
    atom1.type = 0
    atom1.mass = 12.0
    atoms.append(atom1)
    
    drude1 = pygcmc.MCAtom()
    drude1.x = 3.1
    drude1.y = 5.0
    drude1.z = 5.0
    drude1.charge = -1.0
    drude1.type = 1
    drude1.mass = 0.4
    atoms.append(drude1)
    
    # Second atom and Drude
    atom2 = pygcmc.MCAtom()
    atom2.x = 7.0
    atom2.y = 5.0
    atom2.z = 5.0
    atom2.charge = 1.0
    atom2.type = 0
    atom2.mass = 12.0
    atoms.append(atom2)
    
    drude2 = pygcmc.MCAtom()
    drude2.x = 7.0
    drude2.y = 5.1
    drude2.z = 5.0
    drude2.charge = -1.0
    drude2.type = 1
    drude2.mass = 0.4
    atoms.append(drude2)
    
    state.atoms = atoms
    state.activeAtomCount = 4
    
    # Create residues
    residues = []
    for i in range(2):
        res = pygcmc.MCResidue()
        res.atomStart = i * 2
        res.atomCount = 2
        res.active = True
        res.fixed = False
        res.type = 0
        residues.append(res)
    
    state.residues = residues
    state.activeResidueCount = 2
    
    # Set force field
    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 2
    ff.numMovementTypes = 2
    ff.ljSigma = [0.3] * 4
    ff.ljEps = [0.0] * 4
    state.forcefield = ff
    
    # Create Drude force
    drude_force = DrudeForce()
    
    # Add both Drude particles
    charge = -1.0
    polarizability = 0.001
    
    drude_force.addParticle(
        drudeIndex=1, parentIndex=0,
        aniso1Index=-1, aniso2Index=-1,
        aniso3Index=-1, aniso4Index=-1,
        charge=charge, polarizability=polarizability,
        aniso12=1.0, aniso34=1.0
    )
    
    drude_force.addParticle(
        drudeIndex=3, parentIndex=2,
        aniso1Index=-1, aniso2Index=-1,
        aniso3Index=-1, aniso4Index=-1,
        charge=charge, polarizability=polarizability,
        aniso12=1.0, aniso34=1.0
    )
    
    assert drude_force.getNumParticles() == 2
    
    # Calculate energy
    energy = drude_force.calculateEnergySCF(state)
    
    # Should have energy from both Drude oscillators
    assert energy > 0
    
    # Energy should be sum of individual harmonic energies
    # (since particles are far apart, no interaction)
    k = ONE_4PI_EPS0 * charge * charge / polarizability
    energy1 = 0.5 * k * 0.1 * 0.1  # First Drude displaced 0.1 nm in x
    energy2 = 0.5 * k * 0.1 * 0.1  # Second Drude displaced 0.1 nm in y
    expected_total = energy1 + energy2
    
    # Should be close (some difference due to SCF optimization)
    assert abs(energy - expected_total) / expected_total < 0.1
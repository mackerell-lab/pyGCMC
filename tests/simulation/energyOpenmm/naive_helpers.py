# tests/simulation/openmm/naive_comparison.py

import pytest
import math
import pygcmc
import os
import warnings

# Filter SWIG-related DeprecationWarning in advance
warnings.filterwarnings("ignore", category=DeprecationWarning,
                        message="builtin type SwigPyPacked has no __module__ attribute")
warnings.filterwarnings("ignore", category=DeprecationWarning,
                        message="builtin type SwigPyObject has no __module__ attribute")
warnings.filterwarnings("ignore", category=DeprecationWarning,
                        message="builtin type swigvarlink has no __module__ attribute")
from openmm import *
from openmm.app import *
from openmm.unit import *

def create_test_system():
    """Create a simple test system with a few molecules."""
    # Create a system with periodic boundary conditions
    system = System()
    box_size = 3.0 * nanometers
    system.setDefaultPeriodicBoxVectors(
        Vec3(box_size, 0, 0),
        Vec3(0, box_size, 0),
        Vec3(0, 0, box_size)
    )
    
    # Add particles
    # Movement molecule (benzene-like): 6 carbons in a ring
    # Each carbon has charge 0 (for testing LJ only)
    movement_atoms = []
    radius = 0.15  # nm
    for i in range(6):
        angle = i * 2 * math.pi / 6
        x = radius * math.cos(angle)
        y = radius * math.sin(angle)
        z = 0.0
        system.addParticle(12.0)  # mass in amu
        movement_atoms.append([x, y, z])
    
    # Fixed molecule (water-like): 3 atoms
    # All charges set to 0 (for testing LJ only)
    fixed_atoms = [
        [1.0, 0.0, 0.0],  # O
        [1.1, 0.1, 0.0],  # H1
        [1.1, -0.1, 0.0]  # H2
    ]
    system.addParticle(16.0)  # O
    system.addParticle(1.0)   # H1
    system.addParticle(1.0)   # H2
    
    # Create nonbonded force
    nb_force = NonbondedForce()
    
    # Add movement molecule atoms (carbons with OPLS-AA like parameters)
    for _ in range(6):
        nb_force.addParticle(-0.1, 0.34, 0.36)  # charge=-0.1e, sigma=0.34nm, epsilon=0.36kJ/mol
    
    # Add fixed molecule atoms (TIP3P-like parameters)
    nb_force.addParticle(-0.834, 0.3166, 0.650)  # O: charge=-0.834e
    nb_force.addParticle(0.417, 0.0, 0.0)       # H1: charge=0.417e
    nb_force.addParticle(0.417, 0.0, 0.0)       # H2: charge=0.417e
    
    # Set up periodic boundary conditions
    nb_force.setNonbondedMethod(NonbondedForce.NoCutoff)  # Changed to NoCutoff to match naive implementation
    system.addForce(nb_force)
    
    # Create topology
    topology = Topology()
    chain = topology.addChain()
    
    # Add movement molecule (benzene)
    res_movement = topology.addResidue('BEN', chain)
    element = Element.getBySymbol('C')
    for pos in movement_atoms:
        topology.addAtom('C', element, res_movement)
    
    # Add fixed molecule (water)
    res_fixed = topology.addResidue('HOH', chain)
    o_element = Element.getBySymbol('O')
    h_element = Element.getBySymbol('H')
    topology.addAtom('O', o_element, res_fixed)
    topology.addAtom('H1', h_element, res_fixed)
    topology.addAtom('H2', h_element, res_fixed)
    
    # Create positions
    positions = []
    positions.extend([Vec3(*pos) for pos in movement_atoms])
    positions.extend([Vec3(*pos) for pos in fixed_atoms])
    positions = [pos * nanometers for pos in positions]
    
    return system, topology, positions

def calculate_nonbonded_energy(system, positions, movement_atoms, fixed_atoms, use_pbc=False, custom_expression=None):
    """Calculate nonbonded energy between specified groups using direct Coulomb and LJ potentials.
    
    The energy expression includes both electrostatic and van der Waals terms:
    E = kC * q1 * q2 / r +  # direct Coulomb
        4 * sqrt(eps1*eps2) * ((sigma/r)^12 - (sigma/r)^6)  # LJ
    
    Args:
        system: OpenMM System object
        positions: List of Vec3 positions
        movement_atoms: Set of atom indices for movement group
        fixed_atoms: Set of atom indices for fixed group
        use_pbc: Whether to use periodic boundary conditions
        custom_expression: Optional custom energy expression to use instead of default
    """
    # Default nonbonded energy expression, including Coulomb and LJ
    if custom_expression is None:
        energy_expression = """
        kC * q1 * q2 / r + 
        4 * sqrt(eps1*eps2) * (
            (0.5*(sigma1+sigma2)/r)^12 - 
            (0.5*(sigma1+sigma2)/r)^6
        )"""
    else:
        energy_expression = custom_expression
    
    custom_force = CustomNonbondedForce(energy_expression)
    
    # Add parameters for each particle
    custom_force.addPerParticleParameter("q")      # Charge
    custom_force.addPerParticleParameter("sigma")  # LJ sigma
    custom_force.addPerParticleParameter("eps")    # LJ epsilon
    
    # Add global parameters
    custom_force.addGlobalParameter("kC", 138.935456)  # Coulomb constant (kJ·nm/mol/e^2)
    
    # Get parameters from the original NonbondedForce
    nb_force = None
    for force in system.getForces():
        if isinstance(force, NonbondedForce):
            nb_force = force
            break
    if nb_force is None:
        raise ValueError("No NonbondedForce found in system")
    
    # Add particle parameters
    num_particles = system.getNumParticles()
    for i in range(num_particles):
        charge, sigma, epsilon = nb_force.getParticleParameters(i)
        # OpenMM's sigma unit is nm, epsilon unit is kJ/mol
        custom_force.addParticle([charge, sigma, epsilon])
    
    # Only calculate interactions between specified groups
    custom_force.addInteractionGroup(movement_atoms, fixed_atoms)
    
    # Set nonbonded method and cutoff distance
    if use_pbc:
        custom_force.setNonbondedMethod(CustomNonbondedForce.CutoffPeriodic)
        custom_force.setCutoffDistance(1.0 * nanometers)  # Set the same cutoff distance as naive implementation
    else:
        custom_force.setNonbondedMethod(CustomNonbondedForce.CutoffNonPeriodic)
        custom_force.setCutoffDistance(1.0 * nanometers)  # Set the same cutoff distance as naive implementation
    
    # Create energy system containing only the custom force
    energy_system = System()
    for i in range(num_particles):
        energy_system.addParticle(system.getParticleMass(i))
    energy_system.setDefaultPeriodicBoxVectors(*system.getDefaultPeriodicBoxVectors())
    energy_system.addForce(custom_force)
    
    # Use Reference platform to calculate energy
    integrator = VerletIntegrator(0.001 * picoseconds)
    platform = Platform.getPlatformByName('Reference')
    context = Context(energy_system, integrator, platform)
    context.setPositions(positions)
    
    state = context.getState(getEnergy=True)
    energy = state.getPotentialEnergy()
    
    # Print detailed energy information for debugging
    print(f"\nDetailed energy calculation:")
    print(f"Number of movement atoms: {len(movement_atoms)}")
    print(f"Number of fixed atoms: {len(fixed_atoms)}")
    print(f"PBC: {use_pbc}")
    print(f"Energy: {energy.value_in_unit(kilojoules_per_mole):.6f} kJ/mol")
    
    del context, integrator, energy_system
    return energy

def convert_openmm_state_to_mcstate():
    """Convert OpenMM test system to MCState for naive implementation."""
    # Create OpenMM test system
    system, topology, positions = create_test_system()
    
    # Create MCState
    state = pygcmc.MCState()
    
    # Set cutoff distance to 1.0nm, consistent with OpenMM
    state.info.cutoff = 1.0
    
    # 1. Set up force field
    state.forcefield.numTotalTypes = 3  # Three types: C, O, and H
    state.forcefield.numMovementTypes = 1  # C is the movement type
    
    # Get force field parameters from OpenMM system
    nb_force = None
    for force in system.getForces():
        if isinstance(force, NonbondedForce):
            nb_force = force
            break
    
    # Get parameters
    c_params = nb_force.getParticleParameters(0)  # Carbon parameters
    o_params = nb_force.getParticleParameters(6)  # Oxygen parameters
    h_params = nb_force.getParticleParameters(7)  # Hydrogen parameters
    
    # Note: OpenMM's epsilon already includes the factor of 4, so no need to divide by 4
    c_c_eps = c_params[2].value_in_unit(kilojoules_per_mole)
    c_c_sigma = c_params[1].value_in_unit(nanometers)
    o_o_eps = o_params[2].value_in_unit(kilojoules_per_mole)
    o_o_sigma = o_params[1].value_in_unit(nanometers)
    h_h_eps = h_params[2].value_in_unit(kilojoules_per_mole)  # Should be 0
    h_h_sigma = h_params[1].value_in_unit(nanometers)         # Should be 0
    
    # Use OpenMM mixing rules
    def mix_params(eps1, sigma1, eps2, sigma2):
        if eps1 == 0 or eps2 == 0 or sigma1 == 0 or sigma2 == 0:
            return 0.0, 0.0
        # OpenMM mixing rules:
        # - sigma: arithmetic mean 0.5*(sigma1 + sigma2)
        # - epsilon: geometric mean sqrt(eps1 * eps2)
        mixed_sigma = 0.5 * (sigma1 + sigma2)  # Modified here, using 0.5*(sigma1 + sigma2)
        mixed_eps = math.sqrt(eps1 * eps2)
        return mixed_eps, mixed_sigma
    
    # Calculate mixed parameters
    c_o_eps, c_o_sigma = mix_params(c_c_eps, c_c_sigma, o_o_eps, o_o_sigma)
    c_h_eps, c_h_sigma = mix_params(c_c_eps, c_c_sigma, h_h_eps, h_h_sigma)
    o_h_eps, o_h_sigma = mix_params(o_o_eps, o_o_sigma, h_h_eps, h_h_sigma)
    
    # Set force field parameter matrices (numTotalTypes * numTotalTypes = 3 * 3)
    # Complete interaction matrix:
    # [C-C, C-O, C-H]
    # [O-C, O-O, O-H]
    # [H-C, H-O, H-H]
    state.forcefield.ljEps = [
        c_c_eps, c_o_eps, c_h_eps,    # C with (C,O,H) interactions
        c_o_eps, o_o_eps, o_h_eps,    # O with (C,O,H) interactions
        c_h_eps, o_h_eps, h_h_eps     # H with (C,O,H) interactions
    ]
    state.forcefield.ljSigma = [
        c_c_sigma, c_o_sigma, c_h_sigma,    # C with (C,O,H) interactions
        c_o_sigma, o_o_sigma, o_h_sigma,    # O with (C,O,H) interactions
        c_h_sigma, o_h_sigma, h_h_sigma     # H with (C,O,H) interactions
    ]
    
    # Set movement types
    state.movementAtomTypes = [0]  # Type 0 (C) is movement type
    state.numMovementAtomTypes = 1
    
    # 2. Set atoms
    atoms = []
    # Add benzene atoms
    for i in range(6):
        pos = positions[i].value_in_unit(nanometers)
        atom = pygcmc.MCAtom()
        atom.x = pos[0]
        atom.y = pos[1]
        atom.z = pos[2]
        atom.charge = c_params[0].value_in_unit(elementary_charge)
        atom.type = 0  # Carbon type
        atoms.append(atom)
    
    # Add water molecule atoms
    for i in range(6, 9):
        pos = positions[i].value_in_unit(nanometers)
        atom = pygcmc.MCAtom()
        atom.x = pos[0]
        atom.y = pos[1]
        atom.z = pos[2]
        if i == 6:  # Oxygen
            atom.charge = o_params[0].value_in_unit(elementary_charge)
            atom.type = 1  # Oxygen type
        else:  # Hydrogen
            atom.charge = h_params[0].value_in_unit(elementary_charge)
            atom.type = 2  # Hydrogen type (separate from oxygen)
        atoms.append(atom)
    
    state.atoms = atoms
    state.activeAtomCount = len(atoms)
    
    # 3. Set residues
    # Movement residue (benzene)
    movement_res = pygcmc.MCResidue()
    movement_res.active = True
    movement_res.type = 0  # Movement type
    movement_res.atomStart = 0
    movement_res.atomCount = 6
    
    # Fixed residue (water)
    fixed_res = pygcmc.MCResidue()
    fixed_res.active = True
    fixed_res.type = 1  # Fixed type
    fixed_res.atomStart = 6
    fixed_res.atomCount = 3
    
    state.residues = [movement_res, fixed_res]
    state.activeResidueCount = 2
    
    # 4. Set movement residue info
    movement_info = pygcmc.MCMovementResidueInfo()
    movement_info.startIndex = 0
    movement_info.activeCount = 1
    movement_info.totalCount = 1
    
    state.movementResidues = [movement_info]
    
    # 5. Set box and cutoff
    box_vectors = system.getDefaultPeriodicBoxVectors()
    state.info.box = [
        box_vectors[0][0].value_in_unit(nanometers),
        box_vectors[1][1].value_in_unit(nanometers),
        box_vectors[2][2].value_in_unit(nanometers)
    ]
    state.info.cutoff = nb_force.getCutoffDistance().value_in_unit(nanometers)
    
    return state, system, positions


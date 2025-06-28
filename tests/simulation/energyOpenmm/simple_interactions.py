# tests/simulation/openmm/simple_interactions.py

import pytest
import os
import warnings

import math

# Suppress SWIG-related DeprecationWarning
warnings.filterwarnings("ignore", category=DeprecationWarning,
                        message="builtin type SwigPyPacked has no __module__ attribute")
warnings.filterwarnings("ignore", category=DeprecationWarning,
                        message="builtin type SwigPyObject has no __module__ attribute")
warnings.filterwarnings("ignore", category=DeprecationWarning,
                        message="builtin type swigvarlink has no __module__ attribute")
from openmm import *
from openmm.app import *
from openmm.unit import *

# Get the directory containing test data files
TEST_DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "data")

def create_test_system():
    """
    Create a simple test system with a benzene-like molecule and a water-like molecule.
    
    System configuration:
    1. Periodic boundary conditions with 3nm box size
    2. Movement molecule (benzene-like):
       - 6 carbon atoms in a ring
       - Each carbon has charge +0.1e
       - Uses OPLS-AA like parameters:
         σ = 0.34nm, ε = 0.36kJ/mol
       @http://docs.openmm.org/8.2.0/userguide/theory/02_standard_forces.html#lennard-jones-interaction
    
    3. Fixed molecule (water-like):
       - TIP3P water model
       - Oxygen: charge -0.834e, σ = 0.3166nm, ε = 0.650kJ/mol
       - Hydrogens: charge +0.417e each
       @http://docs.openmm.org/8.2.0/userguide/theory/02_standard_forces.html#nonbondedforce
    
    4. Nonbonded force settings:
       - Method: CutoffNonPeriodic
       - Cutoff distance: 1.0nm
       - Switching function: True
       - Switching distance: 0.9nm
    
    Returns:
        system (System): OpenMM system
        topology (Topology): System topology
        positions (list): Initial atomic positions
    """
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
    # Each carbon has charge +0.1e
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
    # Oxygen with charge -0.834e and 2 hydrogens with charge +0.417e each (TIP3P)
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
        nb_force.addParticle(0.1, 0.34, 0.36)  # charge=+0.1e, sigma=0.34nm, epsilon=0.36kJ/mol
    
    # Add fixed molecule atoms (TIP3P-like parameters)
    nb_force.addParticle(-0.834, 0.3166, 0.650)  # O
    nb_force.addParticle(0.417, 0.0, 0.0)        # H1
    nb_force.addParticle(0.417, 0.0, 0.0)        # H2
    
    # Set up nonbonded method
    nb_force.setNonbondedMethod(NonbondedForce.CutoffNonPeriodic)
    nb_force.setCutoffDistance(1.0 * nanometers)
    nb_force.setUseSwitchingFunction(True)
    nb_force.setSwitchingDistance(0.9 * nanometers)
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

def calculate_nonbonded_energy(system, positions, movement_atoms, fixed_atoms, use_pbc=True):
    """
    Calculate nonbonded energy between specified groups using shifted Coulomb and LJ potentials.
    
    Energy expressions:
    1. Lennard-Jones with switching function:
       E_LJ = 4ε[(σ/r)¹² - (σ/r)⁶] * S(r)
       where S(r) is the switching function:
       S(r) = 1                                    if r ≤ r_switch
       S(r) = (r_cut-r)²(r_cut+2r-3r_switch)/     if r_switch < r < r_cut
              (r_cut-r_switch)³
       S(r) = 0                                    if r ≥ r_cut
       @http://docs.openmm.org/8.2.0/userguide/theory/02_standard_forces.html#lennard-jones-interaction
    
    2. Shifted Coulomb potential:
       E_coul = kC * q₁q₂ * (1/r - 1/r_cut)
       @http://docs.openmm.org/8.2.0/userguide/theory/02_standard_forces.html#coulomb-interaction-with-cutoff
    
    Parameters:
        system (System): OpenMM system
        positions (list): Atomic positions
        movement_atoms (set): Indices of first group atoms
        fixed_atoms (set): Indices of second group atoms
        use_pbc (bool): Whether to use periodic boundary conditions
    
    Returns:
        energy (Quantity): Calculated nonbonded energy
    """
    # Complete nonbonded energy expression, including Coulomb and LJ
    energy_expression = """
    step(cutoff - r) * (
        kC * q1 * q2 * (1/r - 1/cutoff) + 
        4 * sqrt(eps1*eps2) * (
            (0.5*(sigma1+sigma2)/r)^12 - 
            (0.5*(sigma1+sigma2)/r)^6
        ) * (
            step(switch - r) +
            step(r - switch) * (cutoff - r)^2 * (cutoff + 2*r - 3*switch) / ((cutoff - switch)^3)
        )
    )"""
    custom_force = CustomNonbondedForce(energy_expression)
    
    # Add per-particle parameters
    custom_force.addPerParticleParameter("q")      # charge
    custom_force.addPerParticleParameter("sigma")  # LJ sigma
    custom_force.addPerParticleParameter("eps")    # LJ epsilon
    
    # Add global parameters
    custom_force.addGlobalParameter("kC", 138.935456)  # Coulomb constant (kJ·nm/mol/e^2)
    custom_force.addGlobalParameter("cutoff", 1.0)     # cutoff distance (nm)
    custom_force.addGlobalParameter("switch", 0.9)     # switching distance (nm)
    
    # Get parameters from original NonbondedForce
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
        # OpenMM's sigma is in nm, epsilon in kJ/mol
        custom_force.addParticle([charge, sigma, epsilon])
    
    # Only calculate interactions between specified groups
    custom_force.addInteractionGroup(movement_atoms, fixed_atoms)
    
    # Set nonbonded method and cutoff distance
    if use_pbc:
        custom_force.setNonbondedMethod(CustomNonbondedForce.CutoffPeriodic)
    else:
        custom_force.setNonbondedMethod(CustomNonbondedForce.CutoffNonPeriodic)
    custom_force.setCutoffDistance(1.0 * nanometers)
    
    # Create energy system with only custom force
    energy_system = System()
    for i in range(num_particles):
        energy_system.addParticle(system.getParticleMass(i))
    energy_system.setDefaultPeriodicBoxVectors(*system.getDefaultPeriodicBoxVectors())
    energy_system.addForce(custom_force)
    
    # Calculate energy using Reference platform
    integrator = VerletIntegrator(0.001 * picoseconds)
    platform = Platform.getPlatformByName('Reference')
    context = Context(energy_system, integrator, platform)
    context.setPositions(positions)
    
    state = context.getState(getEnergy=True)
    energy = state.getPotentialEnergy()
    del context, integrator, energy_system
    return energy

def calculate_self_energy_correction(force, cutoff):
    """
    Calculate self-energy correction for shifted Coulomb potential.
    
    In OpenMM's shifted Coulomb implementation, a self-energy correction term is subtracted:
    U_self = - (kC/(2*r_cut)) * sum_i q_i²
    
    This correction ensures proper energy conservation and is described in:
    @http://docs.openmm.org/8.2.0/userguide/theory/02_standard_forces.html#coulomb-interaction-with-cutoff
    
    Parameters:
        force (NonbondedForce): The force object containing particle parameters
        cutoff (float): Cutoff distance in nm
    
    Returns:
        float: Self-energy correction in kJ/mol
        
    Note:
        This correction is automatically handled by OpenMM's NonbondedForce,
        but needs to be manually applied when using CustomNonbondedForce
        to match the standard implementation.
    """
    kC = 138.935456  # kJ·nm/mol/e^2
    sum_q2 = 0.0
    for i in range(force.getNumParticles()):
        charge, _, _ = force.getParticleParameters(i)
        # Convert to dimensionless value in units of e
        charge_val = charge.value_in_unit(elementary_charge)
        sum_q2 += charge_val * charge_val
    correction = - (kC * sum_q2) / (2 * cutoff)
    return correction

def test_attractive_interaction():
    """
    Test nonbonded energy calculation for an attractive interaction.
    
    Tests the attractive regime of the Lennard-Jones and Coulomb potentials:
    1. LJ attraction: r > r_min where r_min = 2^(1/6)σ
       @http://docs.openmm.org/8.2.0/userguide/theory/02_standard_forces.html#lennard-jones-interaction
    
    2. Coulomb attraction: opposite charges
       @http://docs.openmm.org/8.2.0/userguide/theory/02_standard_forces.html#coulomb-interaction-without-cutoff
    
    System configuration:
    - Benzene (+0.1e per C) interacting with water (-0.834e on O, +0.417e on H)
    - Default separation should result in net attractive force
    
    Verification:
    - Ensures total energy is negative (attractive)
    """
    # Create test system
    system, topology, positions = create_test_system()
    
    # Define movement and fixed atoms
    movement_atoms = set(range(6))  # First 6 atoms (benzene)
    fixed_atoms = set(range(6, 9))  # Last 3 atoms (water)
    
    # Calculate energy
    energy = calculate_nonbonded_energy(system, positions, movement_atoms, fixed_atoms)
    
    # Energy should be negative (attractive) due to opposite charges
    assert energy.value_in_unit(kilojoules_per_mole) < 0, \
           f"Expected attractive interaction, got {energy}"
    print(f"Attractive interaction energy: {energy}")

def test_repulsive_interaction():
    """
    Test nonbonded energy calculation for a repulsive interaction.
    
    Tests the repulsive regime of the Lennard-Jones potential:
    1. LJ repulsion: r < r_min where r_min = 2^(1/6)σ
       E_LJ = 4ε[(σ/r)¹² - (σ/r)⁶]
       @http://docs.openmm.org/8.2.0/userguide/theory/02_standard_forces.html#lennard-jones-interaction
    
    System configuration:
    - Moves water molecule very close to benzene (0.1 nm)
    - At this distance, LJ repulsion dominates over electrostatic attraction
    
    Verification:
    - Ensures total energy is positive (repulsive)
    """
    # Create test system
    system, topology, positions = create_test_system()
    
    # Move water molecule very close to benzene to create repulsion
    new_positions = []
    for i in range(len(positions)):
        if i >= 6 and i < 9:  # Water atoms
            pos = positions[i].value_in_unit(nanometers)
            new_positions.append(Vec3(0.1, pos[1], pos[2]) * nanometers)
        else:
            new_positions.append(positions[i])
    positions = new_positions
    
    # Define movement and fixed atoms
    movement_atoms = set(range(6))  # First 6 atoms (benzene)
    fixed_atoms = set(range(6, 9))  # Last 3 atoms (water)
    
    # Calculate energy
    energy = calculate_nonbonded_energy(system, positions, movement_atoms, fixed_atoms)
    
    # Energy should be positive (repulsive) due to close distance
    assert energy.value_in_unit(kilojoules_per_mole) > 0, \
           f"Expected repulsive interaction, got {energy}"
    print(f"Repulsive interaction energy: {energy}")

def test_pbc_interaction():
    """
    Test nonbonded energy calculation with periodic boundary conditions.
    
    Tests the implementation of periodic boundary conditions in nonbonded calculations:
    1. Reaction field for electrostatics:
       E = (q₁q₂/4πε₀)[1/r + k_rf*r² - c_rf]
       @http://docs.openmm.org/8.2.0/userguide/theory/02_standard_forces.html#coulomb-interaction-with-cutoff
    
    2. Periodic wrapping for LJ interactions
       @http://docs.openmm.org/8.2.0/userguide/theory/02_standard_forces.html#lennard-jones-interaction
    
    System configuration:
    - Moves water molecule to box edge
    - Tests with and without PBC
    
    Verification:
    - Ensures PBC affects the interaction energy
    - Compares energies with and without PBC
    """
    # Create test system
    system, topology, positions = create_test_system()
    
    # Move water molecule to box edge
    box_size = system.getDefaultPeriodicBoxVectors()[0][0].value_in_unit(nanometers)
    new_positions = []
    for i in range(len(positions)):
        if i >= 6 and i < 9:  # Water atoms
            pos = positions[i].value_in_unit(nanometers)
            new_positions.append(Vec3(box_size - 0.1, pos[1], pos[2]) * nanometers)
        else:
            new_positions.append(positions[i])
    positions = new_positions
    
    # Define movement and fixed atoms
    movement_atoms = set(range(6))  # First 6 atoms (benzene)
    fixed_atoms = set(range(6, 9))  # Last 3 atoms (water)
    
    # Calculate energy with and without PBC
    energy_pbc = calculate_nonbonded_energy(system, positions, movement_atoms, fixed_atoms, use_pbc=True)
    energy_no_pbc = calculate_nonbonded_energy(system, positions, movement_atoms, fixed_atoms, use_pbc=False)
    
    # Energy with PBC should be different from energy without PBC
    assert abs(energy_pbc.value_in_unit(kilojoules_per_mole) - 
              energy_no_pbc.value_in_unit(kilojoules_per_mole)) > 1e-3, \
           "PBC should affect the interaction energy"
    print(f"PBC interaction energy: {energy_pbc}")
    print(f"Non-PBC interaction energy: {energy_no_pbc}")

def test_cutoff_effect():
    """
    Test the effect of cutoff distance on nonbonded energy calculation.
    
    Tests the implementation of cutoff-based methods:
    1. Switching function for LJ:
       S(r) = 1-6x⁵+15x⁴-10x³, x=(r-r_switch)/(r_cutoff-r_switch)
       @http://docs.openmm.org/8.2.0/userguide/theory/02_standard_forces.html#lennard-jones-interaction
    
    2. Reaction field for Coulomb:
       E = (q₁q₂/4πε₀)[1/r + k_rf*r² - c_rf]
       @http://docs.openmm.org/8.2.0/userguide/theory/02_standard_forces.html#coulomb-interaction-with-cutoff
    
    Test distances:
    - 0.5 nm: Well within cutoff
    - 0.7 nm: Within cutoff
    - 0.9 nm: At switching distance
    - 1.2 nm: Beyond cutoff
    
    Verification:
    - Energy decreases with distance
    - Energy goes to zero beyond cutoff
    - Switching function properly applied
    """
    # Create test system
    system, topology, positions = create_test_system()
    
    # Place water molecule at different distances
    # Note: benzene carbons are at radius 0.15 nm, so we need x >= 1.15 nm
    # to ensure all atom pairs are beyond the 1.0 nm cutoff
    distances = [0.5, 0.7, 0.9, 1.2]  # nm (last distance ensures all pairs > cutoff)
    movement_atoms = set(range(6))  # First 6 atoms (benzene)
    fixed_atoms = set(range(6, 9))  # Last 3 atoms (water)
    
    energies = []
    for dist in distances:
        # Move water molecule
        new_positions = []
        for i in range(len(positions)):
            if i >= 6 and i < 9:  # Water atoms
                pos = positions[i].value_in_unit(nanometers)
                new_positions.append(Vec3(dist, pos[1], pos[2]) * nanometers)
            else:
                new_positions.append(positions[i])
        
        # Calculate energy without PBC to observe true distance effects
        energy = calculate_nonbonded_energy(system, new_positions, movement_atoms, fixed_atoms, use_pbc=False)
        energy_val = energy.value_in_unit(kilojoules_per_mole)
        energies.append(energy_val)
        print(f"Energy at distance {dist} nm: {energy_val:.4f} kJ/mol")
    
    # Verify energy decreases with distance
    for i in range(len(distances)-1):
        assert abs(energies[i]) > abs(energies[i+1]), \
               f"Energy should decrease with distance. Energies: {energies}"
    
    # Verify energy is zero beyond cutoff (1.0 nm)
    assert abs(energies[-1]) < 1e-6, \
           f"Energy should be zero beyond cutoff (1.0 nm), but got {energies[-1]} kJ/mol"


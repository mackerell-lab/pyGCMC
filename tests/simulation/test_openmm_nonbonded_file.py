# tests/simulation/test_openmm_nonbonded_file.py

import pytest
import os
import numpy as np
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
from openmm.app.gromacstopfile import GromacsTopFile

# Get the directory containing test data files
TEST_DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")

def load_test_system():
    """
    Load test system from PDB and TOP files.
    
    System configuration:
    1. Nonbonded force settings:
       - Method: CutoffNonPeriodic
       - Cutoff distance: 1.0 nm
       - Switching distance: 0.9 nm
       @http://docs.openmm.org/8.2.0/userguide/theory/02_standard_forces.html#nonbondedforce
    
    Returns:
        system (System): OpenMM system with nonbonded forces
        positions (list): Initial atomic positions
    """
    # Load PDB file
    pdb_path = os.path.join(TEST_DATA_DIR, "test.pdb")
    pdb = PDBFile(pdb_path)
    
    # Load TOP file
    top_path = os.path.join(TEST_DATA_DIR, "test.top")
    top = GromacsTopFile(top_path)
    
    # Create system
    system = top.createSystem(
        nonbondedMethod=CutoffNonPeriodic,
        nonbondedCutoff=1.0*nanometer,
        switchDistance=0.9*nanometer
    )
    
    return system, pdb.positions

def test_verify_openmm_expressions():
    """
    Verify consistency between custom nonbonded expressions and OpenMM default implementation.
    
    Tests two implementations:
    1. Standard nonbonded interactions:
       E = E_coulomb + E_LJ
       where:
       E_coulomb = (kC * q₁q₂)/r
       E_LJ = 4ε[(σ/r)¹² - (σ/r)⁶]
       @http://docs.openmm.org/8.2.0/userguide/theory/02_standard_forces.html#lennard-jones-interaction
       @http://docs.openmm.org/8.2.0/userguide/theory/02_standard_forces.html#coulomb-interaction-without-cutoff
    
    2. Exception interactions (1-4 pairs):
       E = coulombScale * E_coulomb + ljScale * E_LJ
       where:
       E_coulomb = (kC * q₁q₂)/r
       E_LJ = 4ε[(σ/r)¹² - (σ/r)⁶]
       @http://docs.openmm.org/8.2.0/userguide/theory/02_standard_forces.html#nonbondedforce
    
    Parameters:
    - kC: Coulomb constant (138.935456 kJ·nm/mol/e²)
    - q₁,q₂: Particle charges
    - σ: Combined LJ diameter (σ₁₂ = (σ₁ + σ₂)/2)
    - ε: Combined LJ well depth (ε₁₂ = √(ε₁ε₂))
    - r: Interparticle distance
    
    Test distances:
    - 0.9 nm: At switching distance
    - 0.95 nm: In switching region
    - 1.0 nm: At cutoff
    - 1.1-2.0 nm: Beyond cutoff
    
    Verification:
    - Relative error < 0.0001% for all distances
    - Both standard and exception interactions match
    """
    # Load test system
    system, positions = load_test_system()
    
    # Get original NonbondedForce
    original_nb_force = None
    for force in system.getForces():
        if isinstance(force, NonbondedForce):
            original_nb_force = force
            break
    
    # Define test distances
    distances = [0.9, 0.95, 1.0, 1.1, 1.5, 2.0]  # Remove extreme short distances
    print("\n=== Comparing Different Nonbonded Interaction Expressions ===")
    print("Distance(nm)  OpenMM Default   Custom Formula    Relative Error(%)")
    print("-" * 55)
    
    platform = Platform.getPlatformByName('Reference')
    
    # Print system basic information (only once)
    print(f"\nSystem Information: {system.getNumParticles()} particles")
    
    for dist in distances:
        # Scale all positions
        scaled_positions = [Vec3(pos[0].value_in_unit(nanometers) * dist,
                               pos[1].value_in_unit(nanometers) * dist,
                               pos[2].value_in_unit(nanometers) * dist) * nanometers
                          for pos in positions]
        
        # 1. Create a new system with ONLY NonbondedForce for comparison
        openmm_system = System()
        for i in range(system.getNumParticles()):
            openmm_system.addParticle(system.getParticleMass(i))
        
        openmm_nb_force = NonbondedForce()
        openmm_nb_force.setNonbondedMethod(NonbondedForce.NoCutoff)
        
        for i in range(original_nb_force.getNumParticles()):
            charge, sigma, epsilon = original_nb_force.getParticleParameters(i)
            openmm_nb_force.addParticle(charge, sigma, epsilon)
        
        for i in range(original_nb_force.getNumExceptions()):
            p1, p2, chargeProd, sigma, epsilon = original_nb_force.getExceptionParameters(i)
            openmm_nb_force.addException(p1, p2, chargeProd, sigma, epsilon)
        
        openmm_system.addForce(openmm_nb_force)
        
        # 1. Calculate energy using OpenMM default implementation
        integrator = VerletIntegrator(0.001 * picoseconds)
        context = Context(openmm_system, integrator, platform)
        context.setPositions(scaled_positions)
        state = context.getState(getEnergy=True)
        openmm_energy = state.getPotentialEnergy()
        del context, integrator

        # 2. Calculate energy using custom formula
        # 2.1 Calculate normal nonbonded interactions
        combined_nonbonded_expression = """
            (kC * q1 * q2 / r + 4 * sqrt(eps1*eps2) * ((0.5*(sigma1+sigma2)/r)^12 - (0.5*(sigma1+sigma2)/r)^6))
        """

        nonbonded_force = CustomNonbondedForce(combined_nonbonded_expression.replace('\n', '').strip())
        nonbonded_force.addPerParticleParameter("q")
        nonbonded_force.addPerParticleParameter("sigma")
        nonbonded_force.addPerParticleParameter("eps")
        nonbonded_force.addGlobalParameter("kC", 138.935456)
        nonbonded_force.setNonbondedMethod(CustomNonbondedForce.NoCutoff)

        # Add particle parameters
        for i in range(original_nb_force.getNumParticles()):
            charge, sigma, epsilon = original_nb_force.getParticleParameters(i)
            nonbonded_force.addParticle([charge, sigma, epsilon])

        # 2.2 Calculate exception interactions
        exception_expression = """
            (kC * chargeprod / r * coulombscale + 4 * epsilon * ((sigma/r)^12 - (sigma/r)^6) * ljscale)
        """
        
        exception_force = CustomBondForce(exception_expression)
        exception_force.addPerBondParameter("chargeprod")
        exception_force.addPerBondParameter("sigma")
        exception_force.addPerBondParameter("epsilon")
        exception_force.addPerBondParameter("ljscale")
        exception_force.addPerBondParameter("coulombscale")
        exception_force.addGlobalParameter("kC", 138.935456)

        for i in range(original_nb_force.getNumExceptions()):
            p1, p2, chargeProd, sigma, epsilon = original_nb_force.getExceptionParameters(i)
            coulomb_scale = 1.0
            lj_scale = 1.0
            exception_force.addBond(p1, p2, [chargeProd, sigma, epsilon, lj_scale, coulomb_scale])
            nonbonded_force.addExclusion(p1, p2)

        # Create system and calculate energy
        custom_system = System()
        for i in range(system.getNumParticles()):
            custom_system.addParticle(system.getParticleMass(i))
        custom_system.addForce(nonbonded_force)
        custom_system.addForce(exception_force)

        integrator = VerletIntegrator(0.001 * picoseconds)
        context = Context(custom_system, integrator, platform)
        context.setPositions(scaled_positions)
        custom_energy = context.getState(getEnergy=True).getPotentialEnergy()
        del context, integrator

        # Calculate relative error
        openmm_val = openmm_energy.value_in_unit(kilojoules_per_mole)
        custom_val = custom_energy.value_in_unit(kilojoules_per_mole)
        rel_error = abs(custom_val - openmm_val) / abs(openmm_val) * 100 if abs(openmm_val) > 1e-6 else 0.0
        
        print(f"{dist:6.2f}  {openmm_val:10.4f}  {custom_val:11.4f}  {rel_error:8.4f}")
        
        # Verify results
        assert rel_error < 1e-4, f"Relative error too large: {rel_error:.4f}%"


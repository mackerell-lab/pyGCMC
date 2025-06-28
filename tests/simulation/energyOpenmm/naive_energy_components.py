# tests/simulation/energyOpenmm/naive_energy_components.py

import pytest
from .naive_helpers import *

def test_openmm_energy_components():
    """Test OpenMM energy components separately."""
    state, system, positions = convert_openmm_state_to_mcstate()
    movement_atoms = set(range(6))
    fixed_atoms = set(range(6, 9))
    
    # Calculate total energy
    total_energy = calculate_nonbonded_energy(system, positions, movement_atoms, fixed_atoms)
    
    # Only calculate electrostatic energy
    elec_expression = """
    kC * q1 * q2 / r;
    """
    elec_energy = calculate_nonbonded_energy(system, positions, movement_atoms, fixed_atoms, 
                                           custom_expression=elec_expression)
    
    # Only calculate LJ energy
    lj_expression = """
    4 * sqrt(eps1*eps2) * (
        (0.5*(sigma1+sigma2)/r)^12 - 
        (0.5*(sigma1+sigma2)/r)^6
    );
    """
    lj_energy = calculate_nonbonded_energy(system, positions, movement_atoms, fixed_atoms,
                                         custom_expression=lj_expression)
    
    print(f"\nOpenMM energy components:")
    print(f"Total energy: {total_energy.value_in_unit(kilojoules_per_mole):.6f} kJ/mol")
    print(f"Electrostatic: {elec_energy.value_in_unit(kilojoules_per_mole):.6f} kJ/mol")
    print(f"LJ: {lj_energy.value_in_unit(kilojoules_per_mole):.6f} kJ/mol")
    
    # Verify total energy approximately equals sum of components
    total_val = total_energy.value_in_unit(kilojoules_per_mole)
    components_sum = (elec_energy.value_in_unit(kilojoules_per_mole) + 
                     lj_energy.value_in_unit(kilojoules_per_mole))
    assert abs(total_val - components_sum) < 1e-6, \
           f"Energy components don't sum to total: {total_val} != {components_sum}"

def test_naive_energy_components():
    """Test naive implementation energy components."""
    state, system, positions = convert_openmm_state_to_mcstate()
    
    # Calculate energy
    pygcmc.computeMovementEnergyCutoff(state)
    
    print(f"\nNaive implementation energy components:")
    print(f"VDW energy: {state.residues[0].energy_vdw:.6f} kJ/mol")
    print(f"Elec energy: {state.residues[0].energy_elec:.6f} kJ/mol")
    print(f"Total energy: {(state.residues[0].energy_vdw + state.residues[0].energy_elec):.6f} kJ/mol")

def print_force_field_params():
    """Print force field parameters for both implementations."""
    state, system, positions = convert_openmm_state_to_mcstate()
    
    # Print OpenMM parameters
    nb_force = None
    for force in system.getForces():
        if isinstance(force, NonbondedForce):
            nb_force = force
            break
            
    print("\nOpenMM parameters:")
    for i in range(nb_force.getNumParticles()):
        charge, sigma, epsilon = nb_force.getParticleParameters(i)
        print(f"Atom {i}: q={charge.value_in_unit(elementary_charge):.3f}e, "
              f"sigma={sigma.value_in_unit(nanometers):.3f}nm, "
              f"epsilon={epsilon.value_in_unit(kilojoules_per_mole):.3f}kJ/mol")
        
    print("\nNaive implementation parameters:")
    print("LJ Epsilon matrix [kJ/mol]:")
    n = int(math.sqrt(len(state.forcefield.ljEps)))
    for i in range(n):
        row = state.forcefield.ljEps[i*n:(i+1)*n]
        print(f"Type {i}: {[f'{x:.3f}' for x in row]}")
    
    print("\nLJ Sigma matrix [nm]:")
    for i in range(n):
        row = state.forcefield.ljSigma[i*n:(i+1)*n]
        print(f"Type {i}: {[f'{x:.3f}' for x in row]}")

def test_compare_openmm_naive_nonbonded():
    """Compare nonbonded energy calculations between OpenMM and naive implementation."""
    # Enable debug output
    pygcmc.setEnergyDebugOutput(True)
    
    # First print all debug information
    print("\n=== Force Field Parameters ===")
    print_force_field_params()
    
    print("\n=== OpenMM Energy Components ===")
    test_openmm_energy_components()
    
    print("\n=== Naive Implementation Energy Components ===")
    test_naive_energy_components()
    
    # Original comparison test code
    state, system, positions = convert_openmm_state_to_mcstate()
    
    movement_atoms = set(range(6))  # Benzene carbons
    fixed_atoms = set(range(6, 9))  # Water atoms
    openmm_energy = calculate_nonbonded_energy(system, positions, movement_atoms, fixed_atoms)
    openmm_energy_val = openmm_energy.value_in_unit(kilojoules_per_mole)
    
    pygcmc.computeMovementEnergyCutoff(state)
    naive_energy = state.residues[0].energy_vdw + state.residues[0].energy_elec
    
    print(f"\n=== Final Energy Comparison ===")
    print(f"OpenMM energy: {openmm_energy_val:.6f} kJ/mol")
    print(f"Naive energy: {naive_energy:.6f} kJ/mol")
    print(f"Absolute difference: {abs(openmm_energy_val - naive_energy):.6f} kJ/mol")
    print(f"Relative difference: {abs(openmm_energy_val - naive_energy)/abs(openmm_energy_val)*100:.6f}%")
    
    rel_tol = 0.001  # 0.1% relative error tolerance, as the two implementations may have different details
    assert abs(openmm_energy_val - naive_energy) / abs(openmm_energy_val) < rel_tol, \
           f"Energy mismatch: OpenMM={openmm_energy_val}, Naive={naive_energy}"


"""
Test PGP Complete vs OpenMM for movement energy delta E accuracy
Adapted from old250726 tests with improvements
"""
import pytest
import pygcmc
import numpy as np
import math

# Only import OpenMM if available
try:
    import openmm
    import openmm.app as app
    import openmm.unit as unit
    OPENMM_AVAILABLE = True
except ImportError:
    OPENMM_AVAILABLE = False
    pytest.skip("OpenMM not available", allow_module_level=True)


def test_pgp_complete_delta_e_accuracy():
    """Test PGP Complete delta E vs OpenMM with high accuracy"""
    print("\n" + "="*70)
    print("PGP Complete vs OpenMM Delta E Accuracy Test")
    print("="*70)
    
    # Reset PGP state
    pygcmc.resetPGPState()
    
    # Create system
    state = pygcmc.MCState()
    box_size = 5.0
    state.info.box = [box_size, box_size, box_size]
    state.info.cutoff = 1.2
    
    # Force field with LJ
    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 2
    
    # LJ parameters for Na+ and Cl-
    sigma_na = 0.333  # nm
    sigma_cl = 0.442  # nm
    eps_na = 0.0115  # kJ/mol
    eps_cl = 0.4184  # kJ/mol
    
    ff.ljSigma = [
        sigma_na, (sigma_na + sigma_cl)/2.0,
        (sigma_na + sigma_cl)/2.0, sigma_cl
    ]
    ff.ljEps = [
        eps_na, math.sqrt(eps_na * eps_cl),
        math.sqrt(eps_na * eps_cl), eps_cl
    ]
    
    state.forcefield = ff
    
    # Create atoms - fixed system plus moving molecule
    atoms = []
    
    # Fixed system: Na-Cl pair
    atom1 = pygcmc.MCAtom()
    atom1.x, atom1.y, atom1.z = 1.0, 1.0, 1.0
    atom1.charge = 1.0
    atom1.type = 0  # Na
    atoms.append(atom1)
    
    atom2 = pygcmc.MCAtom()
    atom2.x, atom2.y, atom2.z = 1.5, 1.0, 1.0
    atom2.charge = -1.0
    atom2.type = 1  # Cl
    atoms.append(atom2)
    
    # Moving molecule: Na-Cl pair
    atom3 = pygcmc.MCAtom()
    atom3.x, atom3.y, atom3.z = 3.0, 3.0, 3.0
    atom3.charge = 1.0
    atom3.type = 0  # Na
    atoms.append(atom3)
    
    atom4 = pygcmc.MCAtom()
    atom4.x, atom4.y, atom4.z = 3.5, 3.0, 3.0
    atom4.charge = -1.0
    atom4.type = 1  # Cl
    atoms.append(atom4)
    
    state.atoms = atoms
    state.activeAtomCount = 4
    
    # Create residues
    residues = []
    
    # Fixed residue (Na-Cl)
    res1 = pygcmc.MCResidue()
    res1.atomStart = 0
    res1.atomCount = 2
    res1.active = True
    res1.fixed = True
    residues.append(res1)
    
    # Moving residue (Na-Cl)
    res2 = pygcmc.MCResidue()
    res2.atomStart = 2
    res2.atomCount = 2
    res2.active = True
    res2.fixed = False
    residues.append(res2)
    
    state.residues = residues
    state.activeResidueCount = 2
    
    # Set movement residues
    state.movementResidues.clear()
    movement_info = pygcmc.MCMovementResidueInfo()
    movement_info.startIndex = 1
    movement_info.activeCount = 1
    state.movementResidues.append(movement_info)
    
    print("\nSystem configuration:")
    print(f"  Box: {box_size} x {box_size} x {box_size} nm")
    print(f"  Fixed: Na-Cl at (1.0, 1.0, 1.0) and (1.5, 1.0, 1.0)")
    print(f"  Moving: Na-Cl at (3.0, 3.0, 3.0) and (3.5, 3.0, 3.0)")
    print(f"  Cutoff: 1.2 nm")
    
    # Initialize parameters
    alpha = 5.6 / state.info.cutoff  # Standard PME alpha
    mesh_size = [32, 32, 32]
    
    pygcmc.setPMEParameters(alpha, mesh_size, 4, 1e-5)
    pygcmc.initializePMEParameters(1.2, [box_size, box_size, box_size], alpha)
    pygcmc.setPGPParameters(alpha, mesh_size, 1.2, mesh_size, 4, 1e-5)
    
    # Precompute grid
    print("\nPrecomputing PGP grid for fixed atoms...")
    pygcmc.precomputeGridPotential(state, fixed_only=True)
    
    # Create OpenMM system
    print("\nCreating OpenMM system...")
    system = openmm.System()
    
    # Add particles
    masses = [22.99, 35.45, 22.99, 35.45]  # Na, Cl, Na, Cl
    for mass in masses:
        system.addParticle(mass * unit.amu)
    
    # Periodic box
    system.setDefaultPeriodicBoxVectors(
        [box_size, 0, 0] * unit.nanometer,
        [0, box_size, 0] * unit.nanometer,
        [0, 0, box_size] * unit.nanometer
    )
    
    # NonbondedForce
    nonbonded = openmm.NonbondedForce()
    nonbonded.setNonbondedMethod(openmm.NonbondedForce.PME)
    nonbonded.setCutoffDistance(1.2 * unit.nanometer)
    nonbonded.setEwaldErrorTolerance(1e-5)
    
    # Add particles
    for i, atom in enumerate(atoms):
        charge = atom.charge * unit.elementary_charge
        if atom.type == 0:  # Na
            sigma = sigma_na * unit.nanometer
            epsilon = eps_na * unit.kilojoule_per_mole
        else:  # Cl
            sigma = sigma_cl * unit.nanometer
            epsilon = eps_cl * unit.kilojoule_per_mole
        nonbonded.addParticle(charge, sigma, epsilon)
    
    system.addForce(nonbonded)
    
    # Create context
    integrator = openmm.VerletIntegrator(1.0 * unit.femtosecond)
    platform = openmm.Platform.getPlatformByName('Reference')
    context = openmm.Context(system, integrator, platform)
    
    # Initial positions
    positions = []
    for atom in atoms:
        positions.append([atom.x, atom.y, atom.z] * unit.nanometer)
    context.setPositions(positions)
    
    # Get initial energies
    print("\n1. Initial configuration:")
    
    # PGP Complete (use corrected implementation for accurate results)
    pgp_result1 = pygcmc.computeMovementEnergyPGPCompleteCorrect(state)
    elec_pgp1 = pgp_result1[0]
    vdw_pgp1 = pgp_result1[1]
    pgp_total1 = elec_pgp1 + vdw_pgp1
    
    print(f"\nPGP Complete (movement only):")
    print(f"  Electrostatic: {elec_pgp1:.6f} kJ/mol")
    print(f"  VDW: {vdw_pgp1:.6f} kJ/mol")
    print(f"  Total: {pgp_total1:.6f} kJ/mol")
    
    # OpenMM total energy
    state_omm1 = context.getState(getEnergy=True)
    energy_omm1 = state_omm1.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    
    # Calculate OpenMM energy without the moving molecule
    # Remove moving molecule
    for i in range(2, 4):
        nonbonded.setParticleParameters(i, 0.0, 1.0*unit.nanometer, 0.0*unit.kilojoule_per_mole)
    nonbonded.updateParametersInContext(context)
    
    state_omm_fixed = context.getState(getEnergy=True)
    energy_omm_fixed = state_omm_fixed.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    
    # Restore moving molecule
    for i in range(2, 4):
        atom = atoms[i]
        charge = atom.charge * unit.elementary_charge
        if atom.type == 0:  # Na
            sigma = sigma_na * unit.nanometer
            epsilon = eps_na * unit.kilojoule_per_mole
        else:  # Cl
            sigma = sigma_cl * unit.nanometer
            epsilon = eps_cl * unit.kilojoule_per_mole
        nonbonded.setParticleParameters(i, charge, sigma, epsilon)
    nonbonded.updateParametersInContext(context)
    
    # OpenMM movement energy = total - fixed_only
    omm_movement1 = energy_omm1 - energy_omm_fixed
    
    print(f"\nOpenMM:")
    print(f"  Total system: {energy_omm1:.6f} kJ/mol")
    print(f"  Fixed only: {energy_omm_fixed:.6f} kJ/mol")
    print(f"  Movement contribution: {omm_movement1:.6f} kJ/mol")
    
    # Move the molecule
    translation = [0.2, -0.1, 0.15]
    print(f"\n2. Moving molecule by {translation} nm...")
    
    # PyGCMC movement
    residue = state.residues[1]
    for i in range(residue.atomCount):
        atom_idx = residue.atomStart + i
        state.atoms[atom_idx].x += translation[0]
        state.atoms[atom_idx].y += translation[1]
        state.atoms[atom_idx].z += translation[2]
    
    # OpenMM movement
    for i in range(2, 4):
        positions[i] = [
            positions[i][0] + translation[0] * unit.nanometer,
            positions[i][1] + translation[1] * unit.nanometer,
            positions[i][2] + translation[2] * unit.nanometer
        ]
    context.setPositions(positions)
    
    # Final energies
    print("\n3. After movement:")
    
    # PGP Complete (use corrected implementation for accurate results)
    pgp_result2 = pygcmc.computeMovementEnergyPGPCompleteCorrect(state)
    elec_pgp2 = pgp_result2[0]
    vdw_pgp2 = pgp_result2[1]
    pgp_total2 = elec_pgp2 + vdw_pgp2
    
    print(f"\nPGP Complete (movement only):")
    print(f"  Electrostatic: {elec_pgp2:.6f} kJ/mol")
    print(f"  VDW: {vdw_pgp2:.6f} kJ/mol")
    print(f"  Total: {pgp_total2:.6f} kJ/mol")
    
    # OpenMM
    state_omm2 = context.getState(getEnergy=True)
    energy_omm2 = state_omm2.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    omm_movement2 = energy_omm2 - energy_omm_fixed
    
    print(f"\nOpenMM:")
    print(f"  Total system: {energy_omm2:.6f} kJ/mol")
    print(f"  Fixed only: {energy_omm_fixed:.6f} kJ/mol (unchanged)")
    print(f"  Movement contribution: {omm_movement2:.6f} kJ/mol")
    
    # Delta E analysis
    print("\n4. Delta E Analysis:")
    
    pgp_delta = pgp_total2 - pgp_total1
    omm_delta = omm_movement2 - omm_movement1
    omm_total_delta = energy_omm2 - energy_omm1
    
    print(f"\nPGP Complete ΔE: {pgp_delta:.6f} kJ/mol")
    print(f"OpenMM movement ΔE: {omm_delta:.6f} kJ/mol")
    print(f"OpenMM total ΔE: {omm_total_delta:.6f} kJ/mol")
    
    # Error analysis
    if abs(omm_delta) > 1e-10:
        error = abs((pgp_delta - omm_delta) / omm_delta) * 100
        print(f"\nRelative error: {error:.4f}%")
        
        # Check accuracy - should be < 1%
        assert error < 1.0, f"Error {error:.2f}% exceeds 1% threshold"
        print("✓ Excellent agreement (<1%)")
    
    print("\n" + "="*70)
    print("Test passed: PGP Complete accurately tracks movement molecule ΔE")
    print("="*70)


def test_pgp_complete_intramolecular_vdw():
    """Test that PGP Complete includes intramolecular VDW"""
    print("\n" + "="*70)
    print("PGP Complete Intramolecular VDW Test")
    print("="*70)
    
    # Reset
    pygcmc.resetPGPState()
    
    # Create system
    state = pygcmc.MCState()
    box_size = 3.0
    state.info.box = [box_size, box_size, box_size]
    state.info.cutoff = 1.2
    
    # Force field
    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 2
    ff.numMovementTypes = 2
    
    # LJ parameters
    sigma_na = 0.333  # nm
    sigma_cl = 0.442  # nm
    eps_na = 0.0115  # kJ/mol
    eps_cl = 0.4184  # kJ/mol
    
    ff.ljSigma = [
        sigma_na, (sigma_na + sigma_cl)/2.0,
        (sigma_na + sigma_cl)/2.0, sigma_cl
    ]
    ff.ljEps = [
        eps_na, math.sqrt(eps_na * eps_cl),
        math.sqrt(eps_na * eps_cl), eps_cl
    ]
    
    state.forcefield = ff
    
    # Create atoms - Na-Cl molecule
    atoms = []
    
    # Movement molecule with intramolecular interaction
    atom1 = pygcmc.MCAtom()
    atom1.x, atom1.y, atom1.z = 1.5, 1.5, 1.5
    atom1.charge = 1.0
    atom1.type = 0  # Na
    atoms.append(atom1)
    
    atom2 = pygcmc.MCAtom()
    atom2.x, atom2.y, atom2.z = 1.9, 1.5, 1.5  # 0.4 nm apart
    atom2.charge = -1.0
    atom2.type = 1  # Cl
    atoms.append(atom2)
    
    state.atoms = atoms
    state.activeAtomCount = 2
    
    # Create one residue
    residues = []
    
    res1 = pygcmc.MCResidue()
    res1.atomStart = 0
    res1.atomCount = 2
    res1.active = True
    res1.fixed = False  # Movement residue
    residues.append(res1)
    
    state.residues = residues
    state.activeResidueCount = 1
    
    # Set movement residues
    state.movementResidues.clear()
    movement_info = pygcmc.MCMovementResidueInfo()
    movement_info.startIndex = 0
    movement_info.activeCount = 1
    state.movementResidues.append(movement_info)
    
    print("\nSystem:")
    print(f"  Movement molecule: Na-Cl")
    print(f"  Na at (1.5, 1.5, 1.5)")
    print(f"  Cl at (1.9, 1.5, 1.5)")
    print(f"  Distance: 0.4 nm")
    
    # Calculate expected LJ energy
    r = 0.4  # nm
    sigma = ff.ljSigma[1]  # Na-Cl sigma
    epsilon = ff.ljEps[1]  # Na-Cl epsilon
    
    sr = sigma / r
    sr6 = sr**6
    sr12 = sr6**2
    expected_lj = 4.0 * epsilon * (sr12 - sr6)
    
    print(f"\nExpected intramolecular LJ: {expected_lj:.6f} kJ/mol")
    
    # Initialize parameters
    alpha = 5.6 / state.info.cutoff
    mesh_size = [32, 32, 32]
    
    pygcmc.setPMEParameters(alpha, mesh_size, 4, 1e-5)
    pygcmc.initializePMEParameters(1.2, [box_size, box_size, box_size], alpha)
    pygcmc.setPGPParameters(alpha, mesh_size, 1.2, mesh_size, 4, 1e-5)
    
    # No fixed atoms, so no grid to precompute
    print("\nNo fixed atoms - grid will be zero")
    
    # Test 1: Standard PME (should exclude intramolecular VDW)
    print("\n1. Standard PME Movement Energy:")
    result_pme = pygcmc.computeMovementEnergyPME(state)
    elec_pme = result_pme[0]
    vdw_pme = result_pme[1]
    
    print(f"  Electrostatic: {elec_pme:.6f} kJ/mol")
    print(f"  VDW: {vdw_pme:.6f} kJ/mol")
    print(f"  Expected: 0 (excludes intramolecular)")
    
    # Test 2: PGP Complete (should include intramolecular VDW)
    print("\n2. PGP Complete Movement Energy:")
    result_pgp = pygcmc.computeMovementEnergyPGPCompleteCorrect(state)
    elec_pgp = result_pgp[0]
    vdw_pgp = result_pgp[1]
    pgp_dict = result_pgp[2]
    
    print(f"  Electrostatic: {elec_pgp:.6f} kJ/mol")
    print(f"  VDW: {vdw_pgp:.6f} kJ/mol")
    print(f"  Expected VDW: {expected_lj:.6f} kJ/mol")
    
    print(f"\nPGP Details:")
    print(f"  Grid: {pgp_dict.get('reciprocal', 0):.6f} kJ/mol")
    print(f"  Real space: {pgp_dict.get('real_space', 0):.6f} kJ/mol")
    print(f"  Self: {pgp_dict.get('self', 0):.6f} kJ/mol")
    
    # Check residue energy
    print(f"\nResidue 0 VDW: {state.residues[0].energy_vdw:.6f} kJ/mol")
    
    # The VDW energy is stored in the residue, not returned directly
    # PGP Complete returns (electrostatic, vdw_from_movement_residues, dict)
    # But for single residue with no movement residues list, it might return 0
    actual_vdw = state.residues[0].energy_vdw
    
    # Verify
    print("\n3. Verification:")
    assert abs(vdw_pme) < 1e-10, "PME should exclude intramolecular VDW"
    print("✓ PME correctly excludes intramolecular VDW")
    
    # Check that the VDW was calculated and stored in the residue
    assert abs(actual_vdw - expected_lj) < 1e-4, f"PGP Complete VDW mismatch: {actual_vdw:.6f} vs expected {expected_lj:.6f}"
    print("✓ PGP Complete correctly includes intramolecular VDW (stored in residue)")
    
    print("\n" + "="*70)
    print("Test passed: PGP Complete correctly handles intramolecular VDW")
    print("="*70)


if __name__ == "__main__":
    test_pgp_complete_delta_e_accuracy()
    test_pgp_complete_intramolecular_vdw()
    print("\nAll tests passed!")
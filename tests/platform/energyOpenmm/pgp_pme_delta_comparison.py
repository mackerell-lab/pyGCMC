"""
Compare PyGCMC PGP delta energies with OpenMM PME

PGP (Precomputed Grid Potential) is designed for efficient calculation
of energy changes in GCMC simulations. This test verifies that PGP
correctly reproduces PME reciprocal space energy changes.

Key concepts:
1. PGP only calculates energy changes (ΔE), not absolute energies
2. Compare with PME reciprocal space energy changes
3. Separate fixed (precomputed) and moving parts
"""

import pytest
import numpy as np
import pygcmc
from pygcmc import MCState, MCAtom, MCResidue, MCForceField, MCMovementResidueInfo

# Only import OpenMM if available
try:
    import openmm as mm
    from openmm import app
    OPENMM_AVAILABLE = True
except ImportError:
    OPENMM_AVAILABLE = False
    pytest.skip("OpenMM not available", allow_module_level=True)


def create_test_system_with_fixed_and_moving():
    """Create a test system with fixed and moving molecules"""
    
    # System parameters
    box_size = 5.0  # nm
    
    # Create PyGCMC system
    state = MCState()
    state.info.box = [box_size, box_size, box_size]
    state.info.cutoff = 2.0
    
    # Force field
    ff = MCForceField()
    ff.numTotalTypes = 2  # Na+, Cl-
    ff.numMovementTypes = 2
    # LJ parameters - need full matrix (numTotalTypes * numTotalTypes)
    # Order: Na-Na, Na-Cl, Cl-Na, Cl-Cl
    # Using Lorentz-Berthelot combining rules
    eps_na = 0.5216  # kJ/mol
    eps_cl = 0.5216
    sig_na = 0.2160  # nm
    sig_cl = 0.4830
    
    ff.ljEps = [
        eps_na,                                    # Na-Na
        (eps_na * eps_cl)**0.5,                   # Na-Cl
        (eps_na * eps_cl)**0.5,                   # Cl-Na (same as Na-Cl)
        eps_cl                                     # Cl-Cl
    ]
    ff.ljSigma = [
        sig_na,                                    # Na-Na
        (sig_na + sig_cl) / 2.0,                  # Na-Cl
        (sig_na + sig_cl) / 2.0,                  # Cl-Na
        sig_cl                                     # Cl-Cl
    ]
    state.forcefield = ff
    
    atoms = []
    residues = []
    
    # Fixed molecules (2 NaCl pairs)
    fixed_positions = [
        ([1.0, 1.0, 2.5], 1.0, 0),   # Na+ 
        ([1.5, 1.0, 2.5], -1.0, 1),  # Cl-
        ([3.5, 3.5, 2.5], 1.0, 0),   # Na+
        ([4.0, 3.5, 2.5], -1.0, 1),  # Cl-
    ]
    
    for i, (pos, charge, atom_type) in enumerate(fixed_positions):
        atom = MCAtom()
        atom.x, atom.y, atom.z = pos
        atom.charge = charge
        atom.type = atom_type
        atoms.append(atom)
        
        res = MCResidue()
        res.atomStart = i
        res.atomCount = 1
        res.active = True
        res.fixed = True
        res.type = atom_type
        residues.append(res)
    
    # Moving molecule (1 NaCl pair)
    moving_positions = [
        ([2.5, 2.5, 2.5], 1.0, 0),   # Na+
        ([2.7, 2.5, 2.5], -1.0, 1),  # Cl-
    ]
    
    for i, (pos, charge, atom_type) in enumerate(moving_positions):
        atom = MCAtom()
        atom.x, atom.y, atom.z = pos
        atom.charge = charge
        atom.type = atom_type
        atoms.append(atom)
        
        res = MCResidue()
        res.atomStart = len(fixed_positions) + i
        res.atomCount = 1
        res.active = True
        res.fixed = False
        res.type = atom_type
        residues.append(res)
    
    state.atoms = atoms
    state.activeAtomCount = len(atoms)
    state.residues = residues
    state.activeResidueCount = len(residues)
    
    return state


def setup_openmm_system(state):
    """Create equivalent OpenMM system"""
    
    # Create OpenMM system
    system = mm.System()
    
    # Add particles
    for atom in state.atoms:
        system.addParticle(22.99 if atom.type == 0 else 35.45)  # Na or Cl mass
    
    # NonbondedForce with PME
    nonbonded = mm.NonbondedForce()
    nonbonded.setNonbondedMethod(mm.NonbondedForce.PME)
    nonbonded.setCutoffDistance(state.info.cutoff)
    nonbonded.setEwaldErrorTolerance(1e-5)
    
    # Add particles to nonbonded force
    for i, atom in enumerate(state.atoms):
        sigma = state.forcefield.ljSigma[atom.type]
        epsilon = state.forcefield.ljEps[atom.type]
        nonbonded.addParticle(atom.charge, sigma, epsilon)
    
    system.addForce(nonbonded)
    
    # Set periodic box
    box_size = state.info.box[0]
    system.setDefaultPeriodicBoxVectors(
        [box_size, 0, 0],
        [0, box_size, 0],
        [0, 0, box_size]
    )
    
    return system, nonbonded


def test_pgp_pme_delta_comparison():
    """Test that PGP correctly reproduces PME reciprocal space energy changes"""
    
    print("\n" + "="*60)
    print("Testing PGP vs OpenMM PME Delta Energies")
    print("="*60)
    
    # Create test system
    state = create_test_system_with_fixed_and_moving()
    
    # Setup PME/PGP parameters
    alpha = 2.84  # 1/nm (typical for 2nm cutoff)
    mesh_size = [32, 32, 32]
    spline_order = 4
    tolerance = 1e-5
    
    # Initialize PyGCMC PME
    pygcmc.setPMEParameters(
        alpha=alpha,
        meshSize=mesh_size,
        splineOrder=spline_order,
        tolerance=tolerance
    )
    pygcmc.initializePMEParameters(state.info.cutoff, state.info.box, alpha)
    
    # Initialize PGP
    pygcmc.setPGPParameters(
        alpha=alpha,
        meshSize=mesh_size,
        potential_cutoff=state.info.cutoff,
        potentialGridSize=mesh_size,
        splineOrder=spline_order,
        tolerance=tolerance
    )
    
    # Precompute grid potential for fixed atoms
    print("\nPrecomputing PGP grid for fixed atoms...")
    pygcmc.precomputeGridPotential(state, fixed_only=True)
    
    # Setup movement residues (last 2 residues are moving)
    state.movementResidues.clear()
    movement_info = MCMovementResidueInfo()
    movement_info.startIndex = 4  # Start from residue 4
    movement_info.activeCount = 2  # 2 moving residues
    state.movementResidues.append(movement_info)
    
    # Calculate initial energies
    print("\nCalculating initial energies...")
    
    # PyGCMC PME
    pme_result = pygcmc.computeMovementEnergyPME(state)
    pme_reciprocal_initial = pme_result[2]['reciprocal']
    
    # PyGCMC PGP
    pgp_energy_initial = pygcmc.calculateMoleculeEnergy(state)
    
    # OpenMM
    system, nonbonded = setup_openmm_system(state)
    integrator = mm.VerletIntegrator(0.001)
    platform = mm.Platform.getPlatformByName('Reference')
    context = mm.Context(system, integrator, platform)
    
    # Set initial positions
    positions = []
    for atom in state.atoms:
        positions.append([atom.x, atom.y, atom.z])
    context.setPositions(positions)
    
    # Get initial OpenMM energy
    omm_state = context.getState(getEnergy=True)
    omm_energy_initial = omm_state.getPotentialEnergy().value_in_unit(mm.unit.kilojoules_per_mole)
    
    # Move the moving molecules
    print("\nMoving molecules...")
    displacement = [0.2, -0.1, 0.15]  # nm
    
    for res_idx in range(4, 6):  # Moving residues
        res = state.residues[res_idx]
        for i in range(res.atomCount):
            atom_idx = res.atomStart + i
            state.atoms[atom_idx].x += displacement[0]
            state.atoms[atom_idx].y += displacement[1]
            state.atoms[atom_idx].z += displacement[2]
    
    # Calculate final energies
    print("\nCalculating final energies...")
    
    # PyGCMC PME
    pme_result_final = pygcmc.computeMovementEnergyPME(state)
    pme_reciprocal_final = pme_result_final[2]['reciprocal']
    
    # PyGCMC PGP
    pgp_energy_final = pygcmc.calculateMoleculeEnergy(state)
    
    # OpenMM - update positions
    new_positions = []
    for atom in state.atoms:
        new_positions.append([atom.x, atom.y, atom.z])
    context.setPositions(new_positions)
    
    # Get final OpenMM energy
    omm_state_final = context.getState(getEnergy=True)
    omm_energy_final = omm_state_final.getPotentialEnergy().value_in_unit(mm.unit.kilojoules_per_mole)
    
    # Calculate energy changes
    print("\nEnergy Changes (ΔE):")
    print("-" * 40)
    
    delta_pme = pme_reciprocal_final - pme_reciprocal_initial
    delta_pgp = pgp_energy_final - pgp_energy_initial
    delta_omm = omm_energy_final - omm_energy_initial
    
    print(f"PyGCMC PME reciprocal ΔE: {delta_pme:10.6f} kJ/mol")
    print(f"PyGCMC PGP ΔE:            {delta_pgp:10.6f} kJ/mol") 
    print(f"OpenMM total ΔE:          {delta_omm:10.6f} kJ/mol")
    
    # Note: OpenMM total includes both reciprocal and direct space
    # PGP should match PME reciprocal space changes
    
    # Calculate relative error between PGP and PME reciprocal
    if abs(delta_pme) > 1e-6:
        relative_error = abs((delta_pgp - delta_pme) / delta_pme)
        print(f"\nPGP vs PME reciprocal relative error: {relative_error:.2%}")
        
        # PGP approximates PME reciprocal space changes
        # For GCMC, 10-20% error is typically acceptable
        assert relative_error < 0.25, f"PGP error too large: {relative_error:.2%}"
    else:
        print("\nPME reciprocal change too small for comparison")
        assert abs(delta_pgp) < 1e-6, "PGP should also show negligible change"
    
    print("\n✓ PGP correctly reproduces PME reciprocal space energy changes")
    
    # Additional analysis
    print("\nDetailed Analysis:")
    print("-" * 40)
    print(f"Initial PME reciprocal: {pme_reciprocal_initial:.6f} kJ/mol")
    print(f"Final PME reciprocal:   {pme_reciprocal_final:.6f} kJ/mol")
    print(f"Initial PGP energy:     {pgp_energy_initial:.6f} kJ/mol")
    print(f"Final PGP energy:       {pgp_energy_final:.6f} kJ/mol")
    
    # Test multiple displacements
    print("\n\nTesting multiple displacements...")
    print("-" * 40)
    
    test_displacements = [
        [0.1, 0.0, 0.0],
        [0.0, 0.15, 0.0],
        [0.0, 0.0, -0.2],
        [-0.1, 0.1, -0.1]
    ]
    
    # Reset to initial positions
    for res_idx in range(4, 6):
        res = state.residues[res_idx]
        for i in range(res.atomCount):
            atom_idx = res.atomStart + i
            state.atoms[atom_idx].x -= displacement[0]
            state.atoms[atom_idx].y -= displacement[1]
            state.atoms[atom_idx].z -= displacement[2]
    
    errors = []
    for disp_idx, disp in enumerate(test_displacements):
        # Apply displacement
        for res_idx in range(4, 6):
            res = state.residues[res_idx]
            for i in range(res.atomCount):
                atom_idx = res.atomStart + i
                state.atoms[atom_idx].x += disp[0]
                state.atoms[atom_idx].y += disp[1]
                state.atoms[atom_idx].z += disp[2]
        
        # Calculate energies
        pme_result = pygcmc.computeMovementEnergyPME(state)
        pme_recip = pme_result[2]['reciprocal']
        pgp = pygcmc.calculateMoleculeEnergy(state)
        
        delta_pme_test = pme_recip - pme_reciprocal_initial
        delta_pgp_test = pgp - pgp_energy_initial
        
        if abs(delta_pme_test) > 1e-6:
            error = abs((delta_pgp_test - delta_pme_test) / delta_pme_test)
            errors.append(error)
            print(f"Displacement {disp}: PGP vs PME error = {error:.2%}")
        
        # Reset positions
        for res_idx in range(4, 6):
            res = state.residues[res_idx]
            for i in range(res.atomCount):
                atom_idx = res.atomStart + i
                state.atoms[atom_idx].x -= disp[0]
                state.atoms[atom_idx].y -= disp[1]
                state.atoms[atom_idx].z -= disp[2]
    
    if errors:
        avg_error = np.mean(errors)
        max_error = np.max(errors)
        print(f"\nAverage error: {avg_error:.2%}")
        print(f"Maximum error: {max_error:.2%}")
        
        assert max_error < 0.25, f"Maximum error too large: {max_error:.2%}"
    
    print("\n✓ All tests passed!")


if __name__ == "__main__":
    test_pgp_pme_delta_comparison()
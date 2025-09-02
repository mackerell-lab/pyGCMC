"""
Test PGP Complete vs OpenMM PME ΔE ratio in reciprocal-dominated regime
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


def test_pgp_openmm_reciprocal_ratio():
    """Test if PGP Complete and OpenMM PME ΔE have a constant ratio"""
    
    print("\n" + "="*60)
    print("Testing PGP Complete vs OpenMM PME ΔE Ratio")
    print("(Reciprocal-dominated regime)")
    print("="*60)
    
    # Create system with large box
    state = MCState()
    box_size = 20.0  # nm - very large box
    state.info.box = [box_size, box_size, box_size]
    state.info.cutoff = 1.0  # nm - small cutoff
    
    # Force field - NO LJ to focus purely on electrostatics
    ff = MCForceField()
    ff.numTotalTypes = 2
    ff.numMovementTypes = 2
    ff.ljEps = [0.0, 0.0, 0.0, 0.0]  # No LJ
    ff.ljSigma = [1.0, 1.0, 1.0, 1.0]
    state.forcefield = ff
    
    atoms = []
    residues = []
    
    # Fixed charges at corners - very far from center
    fixed_positions = [
        ([2.0, 2.0, 2.0], 1.0, 0),      # Na+ corner 1
        ([18.0, 2.0, 2.0], -1.0, 1),    # Cl- corner 2
        ([2.0, 18.0, 2.0], -1.0, 1),    # Cl- corner 3
        ([18.0, 18.0, 18.0], 1.0, 0),   # Na+ corner 4
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
    
    # Moving NaCl in center - far from all fixed atoms
    moving_positions = [
        ([10.0, 10.0, 10.0], 1.0, 0),   # Na+
        ([10.3, 10.0, 10.0], -1.0, 1),  # Cl-
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
    
    # Setup movement residues
    state.movementResidues.clear()
    movement_info = MCMovementResidueInfo()
    movement_info.startIndex = 4
    movement_info.activeCount = 2
    state.movementResidues.append(movement_info)
    
    # Initialize with high alpha for strong reciprocal dominance
    alpha = 5.6 / state.info.cutoff  # = 5.6
    mesh_size = [32, 32, 32]  # Reduced from 64x64x64 for performance
    spline_order = 4
    
    pygcmc.setPMEParameters(alpha, mesh_size, spline_order, 1e-5)
    pygcmc.initializePMEParameters(state.info.cutoff, state.info.box, alpha)
    
    pygcmc.setPGPParameters(alpha, mesh_size, state.info.cutoff, mesh_size, spline_order, 1e-5)
    pygcmc.precomputeGridPotential(state, fixed_only=True)
    
    # Check minimum distances
    print("\nChecking distances (cutoff = 1.0 nm):")
    min_dist = float('inf')
    for i in range(4):  # Fixed atoms
        for j in range(4, 6):  # Moving atoms
            dx = state.atoms[i].x - state.atoms[j].x
            dy = state.atoms[i].y - state.atoms[j].y
            dz = state.atoms[i].z - state.atoms[j].z
            dist = np.sqrt(dx**2 + dy**2 + dz**2)
            min_dist = min(min_dist, dist)
    print(f"  Minimum fixed-moving distance: {min_dist:.2f} nm")
    print(f"  Real space contribution: erfc({alpha:.2f}*{min_dist:.2f}) ≈ {np.exp(-(alpha*min_dist)**2):.2e}")
    
    # Setup OpenMM
    system = mm.System()
    for atom in state.atoms:
        system.addParticle(22.99 if atom.type == 0 else 35.45)
    
    nonbonded = mm.NonbondedForce()
    nonbonded.setNonbondedMethod(mm.NonbondedForce.PME)
    nonbonded.setCutoffDistance(state.info.cutoff)
    nonbonded.setEwaldErrorTolerance(1e-5)
    
    for atom in state.atoms:
        nonbonded.addParticle(atom.charge, 1.0, 0.0)  # No LJ
    
    system.addForce(nonbonded)
    system.setDefaultPeriodicBoxVectors(
        [box_size, 0, 0], [0, box_size, 0], [0, 0, box_size]
    )
    
    integrator = mm.VerletIntegrator(0.001)
    platform = mm.Platform.getPlatformByName('Reference')
    context = mm.Context(system, integrator, platform)
    
    # Test multiple displacements
    print("\n\nTesting ΔE for various displacements:")
    print("-" * 70)
    print("Displacement    PGP Complete ΔE    OpenMM PME ΔE    Ratio    Error%")
    print("-" * 70)
    
    displacements = [
        [0.1, 0.0, 0.0],
        [0.0, 0.2, 0.0],
        [0.2, -0.1, 0.1],
        [0.5, 0.0, 0.0],
    ]  # Reduced from 6 to 4 displacements for performance
    
    ratios = []
    
    # Initial energies
    positions = [[atom.x, atom.y, atom.z] for atom in state.atoms]
    context.setPositions(positions)
    
    # PGP Complete initial
    pgp_result = pygcmc.computeMovementEnergyPGPCompleteCorrect(state)
    pgp_initial = pgp_result[0] + pgp_result[1]  # elec + vdw
    
    # OpenMM initial
    omm_state = context.getState(getEnergy=True)
    omm_initial = omm_state.getPotentialEnergy().value_in_unit(mm.unit.kilojoules_per_mole)
    
    # Note: We compare full-system ΔE directly
    # Setting state.residues[i].active = False does not affect OpenMM Context energies
    # PGP Complete already handles movement residues correctly
    
    for disp in displacements:
        # Apply displacement
        for i in range(4, 6):
            state.atoms[i].x += disp[0]
            state.atoms[i].y += disp[1]
            state.atoms[i].z += disp[2]
            positions[i] = [state.atoms[i].x, state.atoms[i].y, state.atoms[i].z]
        
        context.setPositions(positions)
        
        # PGP Complete final
        pgp_result_final = pygcmc.computeMovementEnergyPGPCompleteCorrect(state)
        pgp_final = pgp_result_final[0] + pgp_result_final[1]
        pgp_delta = pgp_final - pgp_initial
        
        # OpenMM final - compare full system ΔE
        omm_state_final = context.getState(getEnergy=True)
        omm_final = omm_state_final.getPotentialEnergy().value_in_unit(mm.unit.kilojoules_per_mole)
        omm_delta = omm_final - omm_initial
        
        # Calculate ratio
        if abs(omm_delta) > 0.001:
            ratio = pgp_delta / omm_delta
            ratios.append(ratio)
            error = abs(1.0 - ratio) * 100
            print(f"{str(disp):15s} {pgp_delta:15.6f} {omm_delta:15.6f} {ratio:8.4f} {error:8.2f}%")
        else:
            print(f"{str(disp):15s} {pgp_delta:15.6f} {omm_delta:15.6f}    ---      ---")
        
        # Reset positions
        for i in range(4, 6):
            state.atoms[i].x -= disp[0]
            state.atoms[i].y -= disp[1]
            state.atoms[i].z -= disp[2]
            positions[i] = [state.atoms[i].x, state.atoms[i].y, state.atoms[i].z]
    
    if ratios:
        avg_ratio = np.mean(ratios)
        std_ratio = np.std(ratios)
        print("\n" + "-" * 70)
        print(f"Average ratio: {avg_ratio:.4f} ± {std_ratio:.4f}")
        print(f"Relative std: {std_ratio/avg_ratio*100:.2f}%")
    else:
        avg_ratio = None
        
    if avg_ratio is not None:
        # Check if it's close to a simple ratio
        simple_ratios = [0.5, 1.0, 2.0, 3.0, 4.0]
        for sr in simple_ratios:
            if abs(avg_ratio - sr) < 0.1:
                print(f"\nRatio is close to {sr}!")
                break
        
        # Additional analysis
        print("\n\nDetailed analysis:")
        print(f"PGP Complete uses grid interpolation for reciprocal")
        print(f"If ratio ≈ 1: PGP Complete is accurate")
        print(f"If ratio ≠ 1: There may be a systematic factor")
    
    # Assert the test passes
    assert avg_ratio is not None, "Failed to calculate average ratio"
    assert abs(avg_ratio - 1.0) < 0.01, f"Ratio {avg_ratio:.6f} deviates too much from 1.0"
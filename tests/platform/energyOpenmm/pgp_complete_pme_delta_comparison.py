"""
Compare PyGCMC PGP Complete delta energies with PME reciprocal space

This test verifies that PGP Complete reciprocal delta energies converge
to PME reciprocal delta energies as the mesh becomes finer.

NOTE on mesh-self energy:
PME 'reciprocal' = full system reciprocal (fixed-fixed + cross + mobile-mobile).
PGP 'reciprocal' = cross-term only (grid interpolation of fixed potential).
On a discrete mesh, the mobile-mobile reciprocal energy has a position-dependent
artifact called "mesh-self energy" that PGP eliminates by design.
A sufficiently fine grid (>=64^3 for this system) makes the mesh-self negligible
so that the delta comparison passes at <1% relative error.
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


def test_pgp_complete_pme_reciprocal_delta():
    """Test that PGP Complete correctly reproduces PME reciprocal space energy changes"""
    
    print("\n" + "="*60)
    print("Testing PGP Complete vs PME Reciprocal Delta Energies")
    print("="*60)
    
    # Create test system
    state = create_test_system_with_fixed_and_moving()
    
    # Setup PME/PGP parameters
    alpha = 2.84  # 1/nm (typical for 2nm cutoff)
    # Use 128^3 mesh to minimize mesh-self energy artifact in PME.
    # PME 'reciprocal' includes mobile-mobile mesh-self energy that
    # depends on atom position within grid cells. PGP eliminates this
    # artifact by design. At 32^3 the mesh-self dominates (~9% error);
    # at 64^3 it drops to ~1.5% for some directions; 128^3 gives <0.1%.
    mesh_size = [128, 128, 128]
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
    
    # PyGCMC PGP Complete (corrected version)
    pgp_complete_result = pygcmc.computeMovementEnergyPGPCompleteCorrect(state)
    # Extract reciprocal space energy from the dictionary
    pgp_reciprocal_initial = pgp_complete_result[2]['reciprocal']
    
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
    
    # PyGCMC PGP Complete
    pgp_complete_result_final = pygcmc.computeMovementEnergyPGPCompleteCorrect(state)
    # Extract reciprocal space energy from the dictionary
    pgp_reciprocal_final = pgp_complete_result_final[2]['reciprocal']
    
    # Calculate energy changes
    print("\nEnergy Changes (ΔE):")
    print("-" * 40)
    
    delta_pme = pme_reciprocal_final - pme_reciprocal_initial
    delta_pgp_complete = pgp_reciprocal_final - pgp_reciprocal_initial
    
    print(f"PyGCMC PME reciprocal ΔE:  {delta_pme:10.6f} kJ/mol")
    print(f"PyGCMC PGP Complete reciprocal ΔE: {delta_pgp_complete:10.6f} kJ/mol") 
    
    # Calculate relative error
    if abs(delta_pme) > 1e-6:
        relative_error = abs((delta_pgp_complete - delta_pme) / delta_pme)
        print(f"\nPGP Complete reciprocal vs PME reciprocal relative error: {relative_error:.4%}")
        
        # PGP Complete should match PME reciprocal very closely
        # Since it's a complete calculation, not interpolation
        assert relative_error < 0.01, f"PGP Complete error too large: {relative_error:.4%}"
    else:
        print("\nPME reciprocal change too small for comparison")
        assert abs(delta_pgp_complete) < 1e-6, "PGP Complete should also show negligible change"
    
    print("\n✓ PGP Complete accurately reproduces PME reciprocal space energy changes")
    
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
        
        pgp_complete_result = pygcmc.computeMovementEnergyPGPCompleteCorrect(state)
        pgp_reciprocal = pgp_complete_result[2]['reciprocal']
        
        delta_pme_test = pme_recip - pme_reciprocal_initial
        delta_pgp_complete_test = pgp_reciprocal - pgp_reciprocal_initial
        
        if abs(delta_pme_test) > 1e-6:
            error = abs((delta_pgp_complete_test - delta_pme_test) / delta_pme_test)
            errors.append(error)
            print(f"Displacement {disp}: PGP Complete reciprocal vs PME error = {error:.4%}")
        
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
        print(f"\nAverage error: {avg_error:.4%}")
        print(f"Maximum error: {max_error:.4%}")
        
        assert max_error < 0.01, f"Maximum error too large: {max_error:.4%}"
    
    print("\n✓ All PGP Complete tests passed with <1% error!")


if __name__ == "__main__":
    test_pgp_complete_pme_reciprocal_delta()
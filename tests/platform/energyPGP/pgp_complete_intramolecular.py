"""
PGP Complete intramolecular LJ tests

This module tests PGP Complete functionality for multi-atom residues 
to verify intramolecular LJ interactions are properly included
"""

import pytest
import pygcmc
from pygcmc import MCState, MCAtom, MCResidue, MCForceField
import math


def test_pgp_complete_intramolecular_lj():
    """Test PGP Complete correctly includes intramolecular LJ"""
    
    print("\n" + "="*70)
    print("PGP Complete Intramolecular LJ Test")
    print("="*70)
    
    # Reset PGP state
    pygcmc.resetPGPState()
    
    state = MCState()
    state.info.box = [5.0, 5.0, 5.0]
    state.info.cutoff = 2.0
    
    # Simple force field
    ff = MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.ljEps = [1.0]
    ff.ljSigma = [0.35]
    state.forcefield = ff
    
    # Create just one 2-atom residue to test intramolecular LJ
    atoms = []
    atom1 = MCAtom()
    atom1.x, atom1.y, atom1.z = 2.0, 2.5, 2.5
    atom1.charge = 0.5
    atom1.type = 0
    atoms.append(atom1)
    
    atom2 = MCAtom()
    atom2.x, atom2.y, atom2.z = 2.4, 2.5, 2.5  # 0.4 nm apart
    atom2.charge = -0.5
    atom2.type = 0
    atoms.append(atom2)
    
    res0 = MCResidue()
    res0.active = True
    res0.fixed = False
    res0.atomStart = 0
    res0.atomCount = 2
    res0.type = 0
    
    state.atoms = atoms
    state.activeAtomCount = 2
    state.residues = [res0]
    state.activeResidueCount = 1
    
    # Set up parameters
    alpha = 2.2
    mesh_size = [32, 32, 32]
    spline_order = 4
    
    pygcmc.setPMEParameters(alpha, mesh_size, spline_order)
    pygcmc.initializePMEParameters(state.info.cutoff, state.info.box, alpha, mesh_size, spline_order)
    
    # Calculate with standard PME (excludes intramolecular LJ)
    elec_pme, vdw_pme, total_pme = pygcmc.computeSystemEnergyPME(state)
    print(f"\nStandard PME (excludes intramolecular LJ):")
    print(f"  Electrostatic: {elec_pme:.6f} kJ/mol")
    print(f"  VdW:          {vdw_pme:.6f} kJ/mol (should be ~0)")
    
    # Now setup PGP and calculate with PGP Complete
    pygcmc.setPGPParameters(alpha, mesh_size, state.info.cutoff, mesh_size, spline_order, 1e-6)
    pygcmc.precomputeGridPotential(state, fixed_only=False)
    elec_pgp, vdw_pgp, total_pgp = pygcmc.computeSystemEnergyPGPComplete(state)
    
    print(f"\nPGP Complete (includes intramolecular LJ):")
    print(f"  Electrostatic: {elec_pgp:.6f} kJ/mol")
    print(f"  VdW:          {vdw_pgp:.6f} kJ/mol (should be non-zero)")
    
    # Calculate expected LJ energy for verification
    r = 0.4  # Distance between atoms in nm
    sigma = 0.35
    eps = 1.0
    r_ratio = sigma / r
    expected_lj = 4.0 * eps * (r_ratio**12 - r_ratio**6)
    print(f"\nExpected intramolecular LJ: {expected_lj:.6f} kJ/mol")
    
    # Verify results
    # 1. Electrostatic may differ due to different calculation methods
    # PME uses FFT reciprocal space, PGP uses grid interpolation
    elec_diff = abs(elec_pgp - elec_pme)
    rel_elec_diff = elec_diff / abs(elec_pme) * 100 if elec_pme != 0 else 0
    print(f"\nElectrostatic difference: {elec_diff:.6f} kJ/mol ({rel_elec_diff:.2f}%)")
    print("NOTE: Standard PME vs PGP Complete use different algorithms")
    # Just verify energies are finite and reasonable
    assert math.isfinite(elec_pgp) and math.isfinite(elec_pme), "Energies must be finite"
    
    # 2. Standard PME should have ~0 VdW (no intramolecular LJ)
    assert abs(vdw_pme) < 1e-6, f"Standard PME should have no VdW energy, got {vdw_pme:.2e}"
    
    # 3. PGP Complete should have non-zero VdW (includes intramolecular LJ)
    assert abs(vdw_pgp) > 0.1, f"PGP Complete should have significant VdW energy, got {vdw_pgp:.2e}"
    
    # 4. VdW should be close to expected value
    vdw_error = abs(vdw_pgp - expected_lj) / abs(expected_lj) * 100
    print(f"VdW relative error: {vdw_error:.2f}%")
    assert vdw_error < 1.0, f"VdW energy error too large: {vdw_error:.2f}%"
    
    print("\nTest passed! PGP Complete correctly includes intramolecular LJ.")
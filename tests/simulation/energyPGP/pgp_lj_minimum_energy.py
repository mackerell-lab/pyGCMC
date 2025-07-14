# tests/simulation/energyPGP/pgp_lj_minimum_energy.py
"""
Test PGP with atoms at LJ minimum distance to verify energy calculations.

This test ensures that PGP correctly calculates Lennard-Jones interactions
at the potential minimum where energy should equal -epsilon.
"""

import pytest
import math
import pygcmc
from . import pgp_wrapper
from .pgp_wrapper import setPGPParameters, initializePMEParameters, precomputeGridPotential, computeSystemEnergyPGP
from .pgp_wrapper import computeMovementEnergyPGP, computeSystemVdwEnergyCutoff
from pygcmc import MCAtom, MCResidue, MCState

def test_pgp_lj_close_interaction():
    """Test PGP with very close LJ interactions."""
    
    print("\n=== Test: PGP with close LJ interactions ===")
    
    # Create system
    state = MCState()
    
    # System info
    state.info.box = [3.0, 3.0, 3.0]
    state.info.cutoff = 1.2
    state.info.use_switching = False
    
    # Force field
    state.forcefield.numTotalTypes = 1
    state.forcefield.numMovementTypes = 1
    
    # Strong LJ parameters for clear signal
    eps = 2.0  # kJ/mol
    sigma = 0.3  # nm
    state.forcefield.ljEps = [eps]
    state.forcefield.ljSigma = [sigma]
    
    # Create two atoms very close (near LJ minimum)
    atoms = []
    
    # Atom 1 (fixed)
    atom1 = MCAtom()
    atom1.x = 1.5
    atom1.y = 1.5
    atom1.z = 1.5
    atom1.charge = 0.0
    atom1.type = 0
    atoms.append(atom1)
    
    # Atom 2 (moving) - place at LJ minimum distance
    # LJ minimum is at r = 2^(1/6) * sigma ≈ 1.122 * sigma
    r_min = 2**(1/6) * sigma
    atom2 = MCAtom()
    atom2.x = 1.5 + r_min
    atom2.y = 1.5
    atom2.z = 1.5
    atom2.charge = 0.0
    atom2.type = 0
    atoms.append(atom2)
    
    state.atoms = atoms
    state.activeAtomCount = 2
    
    # Residues
    residues = []
    for i in range(2):
        res = MCResidue()
        res.active = True
        res.fixed = (i == 0)  # First is fixed
        res.atomStart = i
        res.atomCount = 1
        res.type = 0
        residues.append(res)
    
    state.residues = residues
    state.activeResidueCount = 2
    
    print(f"Two atoms at distance r = {r_min:.4f} nm (LJ minimum)")
    print(f"Expected LJ energy ≈ -epsilon = {-eps:.4f} kJ/mol")
    
    # Calculate with PGP
    setPGPParameters(
        alpha=0.3, 
        meshSize=[16, 16, 16], 
        potential_cutoff=state.info.cutoff,
        potentialGridSize=[16, 16, 16], 
        splineOrder=4, 
        tolerance=1e-5
    )
    
    initializePMEParameters(state.info.cutoff, state.info.box, 0.3)
    precomputeGridPotential(state)
    
    computeSystemEnergyPGP(state)
    # Verify energy is close to -epsilon
    # Note: Energy is stored in both residues
    lj_total = state.residues[0].energy_vdw + state.residues[1].energy_vdw
    expected_lj = -eps  # At minimum
    relative_error = abs((lj_total - expected_lj) / expected_lj) if expected_lj != 0 else abs(lj_total)
    print(f"Relative error: {relative_error:.2%}")
    
    assert relative_error < 0.05, \
        f"LJ energy {lj_total:.6f} deviates too much from expected {expected_lj:.6f}"
    
    print("\n✅ PGP correctly calculates LJ energy at minimum!")

if __name__ == "__main__":
    test_pgp_lj_only_system()

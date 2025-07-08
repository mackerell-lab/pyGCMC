"""
PME LJ-only analysis tests

Tests for analyzing LJ energy behavior including distance scans
and mixing rule verification.
"""

import pytest
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pygcmc
from pygcmc import MCState, MCAtom, MCResidue, MCForceField
from pygcmc import computeSystemEnergyCutoff

from .pme_lj_only_helpers import calculate_openmm_lj_energy, OPENMM_AVAILABLE


@pytest.mark.skipif(not OPENMM_AVAILABLE, reason="OpenMM not available")
def test_lj_distance_scan():
    """Test LJ energy as a function of distance"""
    
    print("\nLJ energy distance scan:")
    
    box_size = 4.0
    cutoff = 1.5
    
    # Test distances from 0.3 to 1.4 nm
    distances = [0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.4]
    
    for dist in distances:
        state = MCState()
        state.info.box = [box_size, box_size, box_size]
        state.info.cutoff = cutoff
        
        # Simple LJ parameters
        ff = MCForceField()
        ff.numTotalTypes = 1
        ff.numMovementTypes = 1
        ff.ljEps = [1.0]     # kJ/mol (1x1 matrix = 1 value)
        ff.ljSigma = [0.35]  # nm
        state.forcefield = ff
        
        # Two atoms at specified distance
        atoms = []
        residues = []
        
        # Atom 1 at center
        atom1 = MCAtom()
        atom1.x, atom1.y, atom1.z = 2.0, 2.0, 2.0
        atom1.charge = 0.0
        atom1.type = 0
        atoms.append(atom1)
        
        # Atom 2 at distance
        atom2 = MCAtom()
        atom2.x = 2.0 + dist
        atom2.y, atom2.z = 2.0, 2.0
        atom2.charge = 0.0
        atom2.type = 0
        atoms.append(atom2)
        
        for i in range(2):
            res = MCResidue()
            res.active = True
            res.fixed = True
            res.atomStart = i
            res.atomCount = 1
            res.type = 0
            residues.append(res)
        
        state.atoms = atoms
        state.activeAtomCount = 2
        state.residues = residues
        state.activeResidueCount = 2
        
        # Calculate with cutoff
        computeSystemEnergyCutoff(state)
        pygcmc_lj = 0.0
        for res in state.residues:
            if res.active:
                pygcmc_lj += res.energy_vdw
        
        # Calculate with OpenMM
        openmm_lj = calculate_openmm_lj_energy(state, use_pme=False)
        
        # Calculate expected LJ energy
        sigma = ff.ljSigma[0]
        epsilon = ff.ljEps[0]
        r = dist
        sr6 = (sigma/r)**6
        expected_lj = 4.0 * epsilon * (sr6*sr6 - sr6)
        
        print(f"  Distance {dist:.1f} nm: PyGCMC={pygcmc_lj:.6f}, OpenMM={openmm_lj:.6f}, Expected={expected_lj:.6f}")
        
        # Check differences
        pygcmc_vs_expected = abs(pygcmc_lj - expected_lj) / abs(expected_lj) if expected_lj != 0 else 0
        if pygcmc_vs_expected > 0.01:
            print(f"    WARNING: PyGCMC differs from expected by {pygcmc_vs_expected*100:.1f}%")


@pytest.mark.skipif(not OPENMM_AVAILABLE, reason="OpenMM not available")
def test_lj_mixing_rules():
    """Test LJ mixing rules between different atom types"""
    
    print("\nLJ mixing rules test:")
    
    state = MCState()
    state.info.box = [3.0, 3.0, 3.0]
    state.info.cutoff = 1.2
    
    # Two different atom types
    ff = MCForceField()
    ff.numTotalTypes = 2
    ff.numMovementTypes = 2
    # Full 2x2 matrix
    eps0, eps1 = 1.0, 2.0
    sigma0, sigma1 = 0.3, 0.4
    ff.ljEps = [
        eps0,                          # 0-0
        math.sqrt(eps0 * eps1),        # 0-1
        math.sqrt(eps0 * eps1),        # 1-0
        eps1                           # 1-1
    ]
    ff.ljSigma = [
        sigma0,                        # 0-0
        (sigma0 + sigma1) / 2,         # 0-1
        (sigma0 + sigma1) / 2,         # 1-0
        sigma1                         # 1-1
    ]
    state.forcefield = ff
    
    # Two atoms of different types
    atoms = []
    residues = []
    
    atom1 = MCAtom()
    atom1.x, atom1.y, atom1.z = 1.5, 1.5, 1.5
    atom1.charge = 0.0
    atom1.type = 0
    atoms.append(atom1)
    
    atom2 = MCAtom()
    atom2.x, atom2.y, atom2.z = 1.9, 1.5, 1.5  # 0.4 nm apart
    atom2.charge = 0.0
    atom2.type = 1
    atoms.append(atom2)
    
    for i in range(2):
        res = MCResidue()
        res.active = True
        res.fixed = True
        res.atomStart = i
        res.atomCount = 1
        res.type = atoms[i].type
        residues.append(res)
    
    state.atoms = atoms
    state.activeAtomCount = 2
    state.residues = residues
    state.activeResidueCount = 2
    
    # Calculate energies
    computeSystemEnergyCutoff(state)
    pygcmc_lj = 0.0
    for res in state.residues:
        if res.active:
            pygcmc_lj += res.energy_vdw
    
    openmm_lj = calculate_openmm_lj_energy(state, use_pme=False)
    
    # Expected with Lorentz-Berthelot mixing rules
    eps_mixed = math.sqrt(ff.ljEps[0] * ff.ljEps[1])  # Geometric mean
    sigma_mixed = (ff.ljSigma[0] + ff.ljSigma[1]) / 2  # Arithmetic mean
    r = 0.4
    sr6 = (sigma_mixed/r)**6
    expected_lj = 4.0 * eps_mixed * (sr6*sr6 - sr6)
    
    print(f"  Mixed LJ parameters: eps={eps_mixed:.3f}, sigma={sigma_mixed:.3f}")
    print(f"  PyGCMC: {pygcmc_lj:.6f} kJ/mol")
    print(f"  OpenMM: {openmm_lj:.6f} kJ/mol")
    print(f"  Expected (LB rules): {expected_lj:.6f} kJ/mol")
    
    # Check if PyGCMC uses correct mixing rules
    diff_expected = abs(pygcmc_lj - expected_lj) / abs(expected_lj) if expected_lj != 0 else 0
    if diff_expected > 0.01:
        print(f"  ⚠️  PyGCMC differs from Lorentz-Berthelot by {diff_expected*100:.1f}%")
        print("  This suggests PyGCMC may use different mixing rules")


if __name__ == "__main__":
    test_lj_distance_scan()
    test_lj_mixing_rules()
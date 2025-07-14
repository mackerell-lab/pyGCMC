# tests/simulation/energyPGP/pgp_lj_limits.py
"""
Test PGP Lennard-Jones energy at extreme distances.

Verifies correct LJ 12-6 potential behavior at:
- Very short distances (r → σ)
- Near cutoff (r → cutoff)
- Beyond cutoff (r > cutoff → 0)
"""

import pytest
import math
import pygcmc
from . import pgp_wrapper
from .pgp_wrapper import setPGPParameters, initializePMEParameters, precomputeGridPotential, computeSystemEnergyPGP
from pygcmc import MCState, MCAtom, MCResidue
from pygcmc import MCForceField

def create_lj_pair_system(distance, epsilon=1.0, sigma=0.34, box_size=5.0, cutoff=1.2):
    """Create a system with two LJ particles at specified distance."""
    state = MCState()
    state.info.box = [box_size, box_size, box_size]
    state.info.cutoff = cutoff
    
    # Force field with LJ parameters
    ff = MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.ljEps = [epsilon]  # kJ/mol
    ff.ljSigma = [sigma]  # nm
    state.forcefield = ff
    
    # Create atoms with no charge (pure LJ)
    atoms = []
    
    # Atom 1 at center
    atom1 = MCAtom()
    atom1.x = box_size / 2
    atom1.y = box_size / 2
    atom1.z = box_size / 2
    atom1.charge = 0.0  # No charge
    atom1.type = 0
    atoms.append(atom1)
    
    # Atom 2 at specified distance
    atom2 = MCAtom()
    atom2.x = box_size / 2 + distance
    atom2.y = box_size / 2
    atom2.z = box_size / 2
    atom2.charge = 0.0  # No charge
    atom2.type = 0
    atoms.append(atom2)
    
    state.atoms = atoms
    state.activeAtomCount = 2
    
    # Create residues
    residues = []
    for i in range(2):
        res = MCResidue()
        res.active = True
        res.fixed = False
        res.atomStart = i
        res.atomCount = 1
        res.type = 0
        residues.append(res)
    
    state.residues = residues
    state.activeResidueCount = 2
    
    return state

def calculate_lj_analytical(r, epsilon, sigma):
    """Calculate analytical LJ 12-6 energy."""
    if r <= 0:
        return float('inf')
    r_ratio = sigma / r
    return 4 * epsilon * (r_ratio**12 - r_ratio**6)

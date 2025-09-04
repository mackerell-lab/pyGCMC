# tests/simulation/energyOpenmm/nbfix_setup.py
"""
Setup functions and utilities for NBFIX tests

This module provides:
1. Force field creation functions
2. System creation functions
3. Force field application utilities
"""

import pytest
import math
import pygcmc
import os
import warnings

# Filter SWIG-related warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)

try:
    from openmm import *
    from openmm.app import *
    from openmm.unit import *
    HAS_OPENMM = True
except ImportError:
    HAS_OPENMM = False

# Constants
ANGSTROM_TO_NM = 0.1
KCAL_TO_KJ = 4.184
kC = 138.935456  # Coulomb constant in kJ·nm/mol/e^2

# Get test data directory
TEST_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(TEST_DIR)), "data")


def create_nbfix_forcefield():
    """Create a force field with NBFIX parameters from toppar_water_ions.str"""
    ff = pygcmc.ForceField()
    
    # Load the water/ions parameter file
    water_ions_file = os.path.join(DATA_DIR, "toppar_water_ions.str")
    if os.path.exists(water_ions_file):
        pygcmc.PRMParser.parse_file_to_forcefield(water_ions_file, ff)
    
    # Manually add SOD-CLA NBFIX if not present (from literature values)
    # SOD    CLA      -0.083875   3.731
    try:
        result = ff.get_nbfix("SOD", "CLA")
        if not (len(result) >= 3 and result[2]):
            ff.add_nbfix("SOD", "CLA", -0.083875, 3.731)
    except:
        ff.add_nbfix("SOD", "CLA", -0.083875, 3.731)
    
    return ff


def create_simple_ion_water_mcstate():
    """Create a simple Na+ and water MCState for NBFIX testing"""
    # Create MCState directly
    state = pygcmc.MCState()
    
    # Set box and cutoff
    state.info.box = [3.0, 3.0, 3.0]  # 3 nm box
    state.info.cutoff = 1.2  # 1.2 nm cutoff
    
    # Define atom types
    sod_idx = state.atomTypes.get_or_add_type("SOD")  # 0
    ot_idx = state.atomTypes.get_or_add_type("OT")   # 1
    ht_idx = state.atomTypes.get_or_add_type("HT")   # 2
    cla_idx = state.atomTypes.get_or_add_type("CLA")  # 3
    
    # Add atoms
    atoms = []
    
    # Sodium ion
    atom_na = pygcmc.MCAtom()
    atom_na.x = 1.5  # nm
    atom_na.y = 1.5
    atom_na.z = 1.5
    atom_na.charge = 1.0
    atom_na.type = sod_idx  # SOD
    atoms.append(atom_na)
    
    # Water oxygen
    atom_o = pygcmc.MCAtom()
    atom_o.x = 1.8  # 0.3 nm from Na+
    atom_o.y = 1.5
    atom_o.z = 1.5
    atom_o.charge = -0.834
    atom_o.type = ot_idx  # OT
    atoms.append(atom_o)
    
    # Water hydrogen 1
    atom_h1 = pygcmc.MCAtom()
    atom_h1.x = 1.8957
    atom_h1.y = 1.5
    atom_h1.z = 1.5
    atom_h1.charge = 0.417
    atom_h1.type = ht_idx  # HT
    atoms.append(atom_h1)
    
    # Water hydrogen 2
    atom_h2 = pygcmc.MCAtom()
    atom_h2.x = 1.8
    atom_h2.y = 1.5957
    atom_h2.z = 1.5
    atom_h2.charge = 0.417
    atom_h2.type = ht_idx  # HT
    atoms.append(atom_h2)
    
    # Chloride ion
    atom_cl = pygcmc.MCAtom()
    atom_cl.x = 1.5
    atom_cl.y = 1.5
    atom_cl.z = 2.0  # 0.5 nm from Na+
    atom_cl.charge = -1.0
    atom_cl.type = cla_idx  # CLA
    atoms.append(atom_cl)
    
    state.atoms = atoms
    state.activeAtomCount = len(atoms)
    
    # Add residues
    residues = []
    
    # Sodium residue
    res_na = pygcmc.MCResidue()
    res_na.active = True
    res_na.atomStart = 0
    res_na.atomCount = 1
    res_na.type = 0
    residues.append(res_na)
    
    # Water residue
    res_water = pygcmc.MCResidue()
    res_water.active = True
    res_water.atomStart = 1
    res_water.atomCount = 3
    res_water.type = 1
    residues.append(res_water)
    
    # Chloride residue
    res_cl = pygcmc.MCResidue()
    res_cl.active = True
    res_cl.atomStart = 4
    res_cl.atomCount = 1
    res_cl.type = 2
    residues.append(res_cl)
    
    state.residues = residues
    state.activeResidueCount = len(residues)
    
    return state


def setup_mcstate_with_forcefield(state, forcefield):
    """Apply force field parameters to MCState"""
    # Get number of types
    n_types = len(state.atomTypes.atomTypes)
    
    # Initialize forcefield arrays
    state.forcefield.numTotalTypes = n_types
    
    # Create new arrays
    ljSigma = [0.0] * (n_types * n_types)
    ljEps = [0.0] * (n_types * n_types)
    
    # Apply force field parameters with NBFIX
    for i in range(n_types):
        for j in range(n_types):
            type1 = state.atomTypes.atomTypes[i]
            type2 = state.atomTypes.atomTypes[j]
            idx = i * n_types + j
            
            # Check for NBFIX first
            try:
                result = forcefield.get_nbfix(type1, type2)
                if len(result) >= 3 and result[2]:  # has_nbfix
                    # NBFIX found - use those parameters
                    # IMPORTANT: NBFIX parameters in CHARMM files:
                    # - epsilon is negative by convention (we use abs() to make it positive)
                    # - rmin is the actual Rmin value (NOT Rmin/2 like in standard LJ params)
                    epsilon = abs(result[0]) * KCAL_TO_KJ
                    rmin = result[1]  # This is Rmin (not Rmin/2)
                    # Convert Rmin to sigma for standard LJ form: sigma = Rmin / 2^(1/6)
                    sigma = rmin / math.pow(2.0, 1.0/6.0) * ANGSTROM_TO_NM
                else:
                    raise ValueError("No NBFIX")
            except Exception as e:
                # No NBFIX - use standard Lorentz-Berthelot mixing rules
                try:
                    lj1 = forcefield.get_lj_params(type1)
                    lj2 = forcefield.get_lj_params(type2)
                    
                    # Convert CHARMM parameters to standard LJ form
                    # IMPORTANT: Standard LJ params in CHARMM use rmin_half (Rmin/2)
                    # So we multiply by 2 to get Rmin, then convert to sigma
                    sigma1 = 2.0 * lj1.rmin_half / math.pow(2.0, 1.0/6.0) * ANGSTROM_TO_NM
                    sigma2 = 2.0 * lj2.rmin_half / math.pow(2.0, 1.0/6.0) * ANGSTROM_TO_NM
                    eps1 = abs(lj1.epsilon) * KCAL_TO_KJ
                    eps2 = abs(lj2.epsilon) * KCAL_TO_KJ
                    
                    # Lorentz-Berthelot mixing rules:
                    # sigma_ij = (sigma_i + sigma_j) / 2  (arithmetic mean)
                    # epsilon_ij = sqrt(epsilon_i * epsilon_j)  (geometric mean)
                    sigma = (sigma1 + sigma2) / 2.0
                    epsilon = math.sqrt(eps1 * eps2)
                except:
                    sigma = 0.0
                    epsilon = 0.0
            
            ljSigma[idx] = sigma
            ljEps[idx] = epsilon
            
    
    # Assign the complete arrays
    state.forcefield.ljSigma = ljSigma
    state.forcefield.ljEps = ljEps
    
    
    return state

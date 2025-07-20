#!/usr/bin/env python
"""Simple test to verify LJ calculation between water molecules"""

import sys
import numpy as np
sys.path.append('/home/zhaomt/gcmc/test107/pygcmc_dev/build/modules/bindings')

import pygcmc

def test_water_lj_simple():
    print("=== Simple LJ Test for Water ===\n")
    
    # Create state
    state = pygcmc.MCState()
    state.info.cutoff = 2.0
    state.info.box = [10.0, 10.0, 10.0]
    
    # Force field first - avoid segfault
    state.forcefield.numTotalTypes = 1
    state.forcefield.numMovementTypes = 1
    
    # SWM4-NDP LJ parameters
    sigma_O = 0.318395  # nm
    eps_O = 0.21094 * 4.184  # kJ/mol
    
    # Initialize LJ parameters (1x1 matrix for one type)
    state.forcefield.ljSigma = [sigma_O]
    state.forcefield.ljEps = [eps_O]
    
    # Test different distances
    distances = [0.25, 0.3, 0.35, 0.4, 0.5]  # nm
    
    for dist in distances:
        # Clear previous data
        atoms = []
        residues = []
        
        # Create two oxygen atoms in different residues
        # Residue 1 - just oxygen
        o1 = pygcmc.MCAtom()
        o1.x = 0.0
        o1.y = 0.0
        o1.z = 0.0
        o1.charge = 0.0  # No charge for pure LJ test
        o1.type = 0
        atoms.append(o1)
        
        res1 = pygcmc.MCResidue()
        res1.atomStart = 0
        res1.atomCount = 1
        res1.active = True
        res1.fixed = False
        res1.type = 0
        residues.append(res1)
        
        # Residue 2 - just oxygen
        o2 = pygcmc.MCAtom()
        o2.x = dist
        o2.y = 0.0
        o2.z = 0.0
        o2.charge = 0.0
        o2.type = 0
        atoms.append(o2)
        
        res2 = pygcmc.MCResidue()
        res2.atomStart = 1
        res2.atomCount = 1
        res2.active = True
        res2.fixed = False
        res2.type = 0
        residues.append(res2)
        
        # Update state
        state.atoms = atoms
        state.residues = residues
        state.activeAtomCount = 2
        state.activeResidueCount = 2
        
        # Calculate energy
        pygcmc.computeSystemEnergyCutoff(state)
        
        # Get total energy
        total_vdw = sum(res.energy_vdw for res in state.residues)
        total_elec = sum(res.energy_elec for res in state.residues)
        
        # Calculate expected LJ
        sigma_r = sigma_O / dist
        lj_expected = 4 * eps_O * (sigma_r**12 - sigma_r**6)
        
        print(f"Distance: {dist*10:.1f} Å")
        print(f"  Calculated VDW: {total_vdw:.3f} kJ/mol")
        print(f"  Expected LJ: {lj_expected:.3f} kJ/mol")
        print(f"  Difference: {abs(total_vdw - lj_expected):.3e} kJ/mol")
        print()

if __name__ == "__main__":
    test_water_lj_simple()
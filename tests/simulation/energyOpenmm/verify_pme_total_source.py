"""
Verify PME Total value sources

Check differences between state.ewald_energy and function return values
"""

import numpy as np
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pygcmc
from pygcmc import MCState, MCAtom, MCResidue, MCForceField
from pygcmc import initializePMEParameters, computeSystemEnergyPME


def test_pme_total_sources():
    """Test PME Total value consistency (regression test)
    
    Verify whether state.ewald_energy and function return values are consistent
    """
    
    # Create simple system
    state = MCState()
    state.info.box = [5.0, 5.0, 5.0]
    state.info.cutoff = 2.0
    
    # ------------------------------------------------------------
    # Initialize PME only once at the beginning of the function
    # Reuse the same set of global parameters to avoid memory issues from repeated release/reconstruction of grids
    # ------------------------------------------------------------
    alpha = 2.5
    mesh_size = [32, 32, 32]
    spline_order = 4
    try:
        pygcmc.setPMEParameters(alpha, mesh_size, spline_order)
        pygcmc.initializePMEParameters(
            state.info.cutoff,
            state.info.box,
            alpha,
            mesh_size,
            spline_order
        )
    except Exception as e:
        print(f"PME initialization warning: {e}")
        pass
    
    ff = MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.ljEps = [0.5]
    ff.ljSigma = [0.35]
    state.forcefield = ff
    
    # 2 atoms
    positions = [[2.0, 2.5, 2.5], [3.0, 2.5, 2.5]]
    charges = [1.0, -1.0]
    
    atoms = []
    for i in range(2):
        atom = MCAtom()
        atom.x, atom.y, atom.z = positions[i]
        atom.charge = charges[i]
        atom.type = 0
        atoms.append(atom)
    
    state.atoms = atoms
    state.activeAtomCount = 2
    
    # Test different residue configurations
    configs = [
        (1, "1 residue"),
        (2, "2 residues")
    ]
    
    for n_residues, desc in configs:
        print(f"\nTest: {desc}")
        print("-" * 50)
        
        # Set residues
        residues = []
        if n_residues == 1:
            res = MCResidue()
            res.active = True
            res.fixed = False
            res.atomStart = 0
            res.atomCount = 2
            res.type = 0
            residues.append(res)
        else:  # 2 residues
            for i in range(2):
                res = MCResidue()
                res.active = True
                res.fixed = False
                res.atomStart = i
                res.atomCount = 1
                res.type = 0
                residues.append(res)
        
        state.residues = residues
        state.activeResidueCount = n_residues
        
        # Initialize PME
        initializePMEParameters(state.info.cutoff, state.info.box, 2.5)
        
        # Call computeSystemEnergyPME and get return value
        result = computeSystemEnergyPME(state)
        
        # result is a tuple: (electrostatic_total, vdw, pme_dict)
        if isinstance(result, tuple) and len(result) == 3:
            elec_total, vdw, pme_dict = result
            
            print(f"\nFrom function return value:")
            print(f"  Electrostatic total energy: {elec_total:.2f}")
            print(f"  VDW energy: {vdw:.6f}")
            print(f"  Dictionary total: {pme_dict.get('total', 'NOT SET')}")
            print(f"  Dictionary real_space: {pme_dict.get('real_space', 'NOT SET')}")
            print(f"  Dictionary reciprocal: {pme_dict.get('reciprocal', 'NOT SET')}")
            print(f"  Dictionary self: {pme_dict.get('self', 'NOT SET')}")
        
        # From state.ewald_energy
        print(f"\nFrom state.ewald_energy:")
        print(f"  total: {state.ewald_energy.get('total', 'NOT SET')}")
        print(f"  real_space: {state.ewald_energy.get('real_space', 'NOT SET')}")
        print(f"  reciprocal: {state.ewald_energy.get('reciprocal', 'NOT SET')}")
        print(f"  self: {state.ewald_energy.get('self', 'NOT SET')}")
        
        # Manual calculation
        manual_total = (state.ewald_energy.get('real_space', 0) + 
                       state.ewald_energy.get('reciprocal', 0) + 
                       state.ewald_energy.get('self', 0))
        
        print(f"\nManually calculated electrostatic total energy: {manual_total:.2f}")
        
        # Check differences
        if isinstance(result, tuple) and len(result) == 3:
            _, _, pme_dict = result
            dict_total = pme_dict.get('total', 0)
            state_total = state.ewald_energy.get('total', 0)
            
            print(f"\nDifference analysis:")
            print(f"  Return dictionary total: {dict_total:.2f}")
            print(f"  State total: {state_total:.2f}")
            print(f"  Difference: {abs(dict_total - state_total):.2f}")
            
            if abs(dict_total - state_total) > 1e-6:
                print("  ⚠️  Return value and state total are inconsistent!")
    
    print("\n" + "="*80)


if __name__ == "__main__":
    test_pme_total_sources()
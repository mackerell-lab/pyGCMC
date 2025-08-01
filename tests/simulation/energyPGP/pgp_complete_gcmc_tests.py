"""
PGP Complete GCMC-specific tests

This module contains tests for PGP Complete under typical GCMC conditions.
"""

import pytest
import pygcmc
from pygcmc import MCState, MCAtom, MCResidue, MCForceField, MCMovementResidueInfo
import math


def test_pgp_complete_typical_gcmc_conditions():
    """Test PGP Complete under typical GCMC conditions"""
    
    # Reset PGP state
    pygcmc.resetPGPState()
    
    print("\n" + "="*70)
    print("PGP Complete Typical GCMC Conditions Test")
    print("="*70)
    
    # Test different box shapes (all reasonable for GCMC)
    test_configs = [
        ([5.0, 5.0, 5.0], "Cubic", [32, 32, 32]),
        ([4.0, 4.0, 6.0], "Slightly elongated", [32, 32, 64]),  # 48->64 for power of 2
        ([6.0, 4.0, 5.0], "Rectangular", [64, 32, 32]),  # 48,40->64,32 for power of 2
    ]
    
    for box, desc, mesh_size in test_configs:
        print(f"\nTesting {desc} box: {box}")
        
        # Reset
        pygcmc.resetPGPState()
        
        # Create state
        state = MCState()
        state.info.box = box
        state.info.cutoff = min(1.2, min(box) * 0.4)
        
        # Force field
        ff = MCForceField()
        ff.numTotalTypes = 2
        
        sigma_na = 0.333
        sigma_cl = 0.442
        eps_na = 0.0115
        eps_cl = 0.4184
        
        ff.ljSigma = [
            sigma_na, (sigma_na + sigma_cl)/2.0,
            (sigma_na + sigma_cl)/2.0, sigma_cl
        ]
        ff.ljEps = [
            eps_na, math.sqrt(eps_na * eps_cl),
            math.sqrt(eps_na * eps_cl), eps_cl
        ]
        
        state.forcefield = ff
        
        # Place atoms at reasonable positions
        atoms = []
        
        # Fixed Na-Cl near center
        fixed_pos = [box[0]*0.4, box[1]*0.4, box[2]*0.4]
        
        atom = MCAtom()
        atom.x, atom.y, atom.z = fixed_pos[0], fixed_pos[1], fixed_pos[2]
        atom.charge = 1.0
        atom.type = 0
        atoms.append(atom)
        
        atom = MCAtom()
        atom.x, atom.y, atom.z = fixed_pos[0] + 0.3, fixed_pos[1], fixed_pos[2]
        atom.charge = -1.0
        atom.type = 1
        atoms.append(atom)
        
        # Moving Na-Cl at reasonable distance
        distance = 0.8 * state.info.cutoff
        moving_pos = [fixed_pos[0] + distance, fixed_pos[1], fixed_pos[2]]
        
        atom = MCAtom()
        atom.x, atom.y, atom.z = moving_pos[0], moving_pos[1], moving_pos[2]
        atom.charge = 1.0
        atom.type = 0
        atoms.append(atom)
        
        atom = MCAtom()
        atom.x, atom.y, atom.z = moving_pos[0] + 0.3, moving_pos[1], moving_pos[2]
        atom.charge = -1.0
        atom.type = 1
        atoms.append(atom)
        
        state.atoms = atoms
        state.activeAtomCount = 4
        
        # Residues
        residues = []
        
        res = MCResidue()
        res.atomStart = 0
        res.atomCount = 2
        res.active = True
        res.fixed = True
        residues.append(res)
        
        res = MCResidue()
        res.atomStart = 2
        res.atomCount = 2
        res.active = True
        res.fixed = False
        residues.append(res)
        
        state.residues = residues
        state.activeResidueCount = 2
        
        state.movementResidues.clear()
        movement_info = MCMovementResidueInfo()
        movement_info.startIndex = 1
        movement_info.activeCount = 1
        state.movementResidues.append(movement_info)
        
        # Initialize PGP
        alpha = 5.6 / state.info.cutoff
        
        try:
            pygcmc.setPMEParameters(alpha, mesh_size, 4, 1e-5)
            pygcmc.initializePMEParameters(state.info.cutoff, state.info.box, alpha)
            pygcmc.setPGPParameters(alpha, mesh_size, state.info.cutoff, mesh_size, 4, 1e-5)
            pygcmc.precomputeGridPotential(state, fixed_only=True)
            
            # Test small displacement
            orig_pos = [(atom.x, atom.y, atom.z) for atom in state.atoms]
            
            # PGP before
            pgp_result1 = pygcmc.computeMovementEnergyPGPCompleteCorrect(state)
            pgp_e1 = pgp_result1[0] + pgp_result1[1]
            
            # Move atoms
            for i in [2, 3]:  # Moving atoms
                state.atoms[i].x += 0.05
            
            # PGP after
            pgp_result2 = pygcmc.computeMovementEnergyPGPCompleteCorrect(state)
            pgp_e2 = pgp_result2[0] + pgp_result2[1]
            
            pgp_delta = pgp_e2 - pgp_e1
            
            print(f"  Initial energy: {pgp_e1:.3f} kJ/mol")
            print(f"  Final energy: {pgp_e2:.3f} kJ/mol")
            print(f"  Delta E: {pgp_delta:.3f} kJ/mol")
            
            # Check results are reasonable
            assert math.isfinite(pgp_e1), "Initial energy must be finite"
            assert math.isfinite(pgp_e2), "Final energy must be finite"
            assert abs(pgp_delta) < 100, "Delta E should be reasonable"
            
            print(f"  ✅ {desc} box test PASSED")
            
        except Exception as e:
            print(f"  ❌ Error: {str(e)}")
            raise
    
    print("\n✅ All typical GCMC condition tests PASSED")
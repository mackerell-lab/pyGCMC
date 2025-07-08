"""
Movement residue PME energy calculation tests

Tests the computeMovementEnergyPME function for GCMC applications
"""

import pytest
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pygcmc
from pygcmc import MCState, MCAtom, MCResidue, MCForceField, MCMovementResidueInfo
from pygcmc import initializePMEParameters, computeSystemEnergyPME, computeMovementEnergyPME


def test_movement_residue_pme():
    """Test movement residue PME energy calculation"""
    
    # Create a simple system with one fixed charge and one moveable charge
    box_size = 3.0
    cutoff = 1.0
    
    state = MCState()
    state.info.box = [box_size, box_size, box_size]
    state.info.cutoff = cutoff
    
    # Force field - pure electrostatics
    ff = MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.ljEps = [0.0]
    ff.ljSigma = [0.3]
    state.forcefield = ff
    
    atoms = []
    residues = []
    
    # Fixed positive charge
    atom = MCAtom()
    atom.x, atom.y, atom.z = 1.0, 1.5, 1.5
    atom.charge = 1.0
    atom.type = 0
    atoms.append(atom)
    
    res = MCResidue()
    res.active = True
    res.fixed = True
    res.atomStart = 0
    res.atomCount = 1
    res.type = 0
    residues.append(res)
    
    # Movement residue - negative charge
    atom = MCAtom()
    atom.x, atom.y, atom.z = 2.0, 1.5, 1.5
    atom.charge = -1.0
    atom.type = 0
    atoms.append(atom)
    
    res = MCResidue()
    res.active = True
    res.fixed = False
    res.atomStart = 1
    res.atomCount = 1
    res.type = 0
    residues.append(res)
    
    state.atoms = atoms
    state.activeAtomCount = 2
    state.residues = residues
    state.activeResidueCount = 2
    
    # Set movement residue
    movement_info = MCMovementResidueInfo()
    movement_info.startIndex = 1
    movement_info.activeCount = 1
    state.movementResidues = [movement_info]
    
    # Initialize PME
    alpha = 5.0
    mesh_size = [24, 24, 24]
    spline_order = 4
    
    pygcmc.setPMEParameters(alpha, mesh_size, spline_order)
    initializePMEParameters(state.info.cutoff, state.info.box, alpha, mesh_size, spline_order)
    
    # Calculate system energy
    computeSystemEnergyPME(state)
    system_energy = state.ewald_energy.get('total', 0.0)
    
    # Calculate movement energy
    movement_result = computeMovementEnergyPME(state)
    
    # Extract electrostatic component
    if isinstance(movement_result, tuple) and len(movement_result) > 0:
        movement_elec = movement_result[0]
    else:
        movement_elec = 0.0
    
    print(f"\nMovement residue PME test:")
    print(f"  System: +1 charge (fixed) and -1 charge (movement)")
    print(f"  Distance: 1.0 nm")
    print(f"  System energy: {system_energy:.6f} kJ/mol")
    print(f"  Movement energy: {movement_elec:.6f} kJ/mol")
    
    # Basic checks
    assert movement_elec != 0.0, "Movement energy should be non-zero"
    
    # Note: computeMovementEnergyPME returns the reciprocal space contribution
    # of the movement residue, not the total interaction energy
    # This can be positive even for attractive interactions
    print(f"\nNote: Movement energy is the reciprocal space contribution")
    print(f"      of the movement residue, not the total interaction energy")
    
    print(f"\n✓ Movement residue PME calculation successful")


def test_movement_residue_removal():
    """Test energy change when removing movement residue"""
    
    # Create system with two charges
    box_size = 3.0
    cutoff = 1.0
    
    state = MCState()
    state.info.box = [box_size, box_size, box_size]
    state.info.cutoff = cutoff
    
    ff = MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.ljEps = [0.0]
    ff.ljSigma = [0.3]
    state.forcefield = ff
    
    # Two opposite charges
    atoms = []
    residues = []
    
    # Fixed charge
    atom1 = MCAtom()
    atom1.x, atom1.y, atom1.z = 1.5, 1.5, 1.5
    atom1.charge = 1.0
    atom1.type = 0
    atoms.append(atom1)
    
    res1 = MCResidue()
    res1.active = True
    res1.fixed = True
    res1.atomStart = 0
    res1.atomCount = 1
    res1.type = 0
    residues.append(res1)
    
    # Movement charge
    atom2 = MCAtom()
    atom2.x, atom2.y, atom2.z = 1.8, 1.5, 1.5
    atom2.charge = -1.0
    atom2.type = 0
    atoms.append(atom2)
    
    res2 = MCResidue()
    res2.active = True
    res2.fixed = False
    res2.atomStart = 1
    res2.atomCount = 1
    res2.type = 0
    residues.append(res2)
    
    state.atoms = atoms
    state.activeAtomCount = 2
    state.residues = residues
    state.activeResidueCount = 2
    
    # Set movement residue
    movement_info = MCMovementResidueInfo()
    movement_info.startIndex = 1
    movement_info.activeCount = 1
    state.movementResidues = [movement_info]
    
    # PME setup
    alpha = 5.0
    mesh_size = [24, 24, 24]
    spline_order = 4
    
    pygcmc.setPMEParameters(alpha, mesh_size, spline_order)
    initializePMEParameters(state.info.cutoff, state.info.box, alpha, mesh_size, spline_order)
    
    # Energy with both charges
    computeSystemEnergyPME(state)
    energy_with_both = state.ewald_energy.get('total', 0.0)
    
    # Deactivate movement residue
    state.residues[1].active = False
    state.activeResidueCount = 1
    state.activeAtomCount = 1
    
    # Energy with only fixed charge
    computeSystemEnergyPME(state)
    energy_fixed_only = state.ewald_energy.get('total', 0.0)
    
    # Energy difference should be the interaction energy
    energy_diff = energy_with_both - energy_fixed_only
    
    print(f"\nMovement residue removal test:")
    print(f"  Energy with both charges: {energy_with_both:.6f} kJ/mol")
    print(f"  Energy with fixed only: {energy_fixed_only:.6f} kJ/mol")
    print(f"  Energy difference: {energy_diff:.6f} kJ/mol")
    
    # Verify results
    assert energy_with_both < energy_fixed_only, "System with opposite charges should have lower energy"
    assert energy_diff < 0.0, "Removing attractive interaction should increase energy"
    
    print(f"\n✓ Movement residue removal test successful")


if __name__ == "__main__":
    test_movement_residue_pme()
    test_movement_residue_removal()
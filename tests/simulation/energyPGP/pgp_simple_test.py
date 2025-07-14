# tests/simulation/energyPGP/pgp_simple_test.py
"""
Simple test to verify PGP is working correctly after the fix.
"""

import pygcmc
from . import pgp_wrapper
from .pgp_wrapper import (
    initializePMEParameters, setPGPParameters, precomputeGridPotential,
    computeSystemEnergyPGP, computeSystemEnergyPME, setPMEParameters
)
from pygcmc import MCState, MCAtom, MCResidue
from pygcmc import MCForceField

def test_pgp_simple():
    """Simple test with charged particles - PGP vs PME movement energy."""
    print("\n=== Simple PGP Test ===")
    
    # Create system with more fixed atoms for better PGP performance
    state = MCState()
    state.info.box = [5.0, 5.0, 5.0]
    state.info.cutoff = 1.2
    
    # Force field
    ff = MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.ljEps = [0.0]
    ff.ljSigma = [0.3]
    state.forcefield = ff
    
    atoms = []
    residues = []
    
    # Create 4 fixed atoms in a square
    fixed_positions = [
        (2.0, 2.0, 2.5, 0.5),   # +0.5 charge
        (3.0, 2.0, 2.5, -0.5),  # -0.5 charge
        (2.0, 3.0, 2.5, -0.5),  # -0.5 charge
        (3.0, 3.0, 2.5, 0.5),   # +0.5 charge
    ]
    
    for i, (x, y, z, charge) in enumerate(fixed_positions):
        atom = MCAtom()
        atom.x = x
        atom.y = y
        atom.z = z
        atom.charge = charge
        atom.type = 0
        atoms.append(atom)
        
        res = MCResidue()
        res.active = True
        res.fixed = True
        res.atomStart = i
        res.atomCount = 1
        res.type = 0
        residues.append(res)
    
    # Add one moveable atom
    moveable = MCAtom()
    moveable.x = 2.5
    moveable.y = 2.5
    moveable.z = 2.5
    moveable.charge = 1.0
    moveable.type = 0
    atoms.append(moveable)
    
    moveable_res = MCResidue()
    moveable_res.active = True
    moveable_res.fixed = False
    moveable_res.atomStart = 4
    moveable_res.atomCount = 1
    moveable_res.type = 0
    residues.append(moveable_res)
    
    state.atoms = atoms
    state.activeAtomCount = 5
    state.residues = residues
    state.activeResidueCount = 5
    
    # Initialize
    alpha = 2.5
    mesh_size = [32, 32, 32]
    
    # Set up movement residues
    from pygcmc import MCMovementResidueInfo
    state.movementResidues = []
    movement_info = MCMovementResidueInfo()
    movement_info.startIndex = 4  # The moveable residue
    movement_info.activeCount = 1
    state.movementResidues.append(movement_info)
    
    # Initialize parameters
    initializePMEParameters(state.info.cutoff, state.info.box, alpha)
    setPMEParameters(alpha, mesh_size, 4, 1e-6)
    setPGPParameters(alpha, mesh_size, state.info.cutoff, mesh_size, 4, 1e-6)
    
    # PGP calculation
    print("\nPGP calculation:")
    precomputeGridPotential(state, fixed_only=True)
    
    # Calculate initial energy
    pgp_initial = pgp_wrapper.calculateMoleculeEnergy(state)
    print(f"  Initial PGP energy: {pgp_initial:.4f} kJ/mol")
    
    # PME movement energy calculation for comparison
    print("\nPME movement energy calculation:")
    pme_move = pgp_wrapper.computeMovementEnergyPME(state)
    pme_move_elec = pme_move[0]  # Electrostatic component
    print(f"  PME movement energy: {pme_move_elec:.4f} kJ/mol")
    
    # Move the particle slightly
    print("\nMoving particle by 0.1 nm in x direction...")
    state.atoms[4].x += 0.1
    
    # Calculate new energies
    pgp_moved = pgp_wrapper.calculateMoleculeEnergy(state)
    pme_move_new = pgp_wrapper.computeMovementEnergyPME(state)
    pme_move_elec_new = pme_move_new[0]
    
    # Calculate energy changes
    pgp_delta = pgp_moved - pgp_initial
    pme_delta = pme_move_elec_new - pme_move_elec
    
    print(f"\nEnergy changes:")
    print(f"  PGP delta:       {pgp_delta:.4f} kJ/mol")
    print(f"  PME move delta:  {pme_delta:.4f} kJ/mol")
    print(f"  Difference:      {abs(pgp_delta - pme_delta):.4f} kJ/mol")
    
    # For pure electrostatic systems, PGP and PME movement energy changes should be similar
    # Allow some tolerance due to different algorithms
    assert abs(pgp_delta - pme_delta) < 5.0, f"Energy change mismatch: PGP={pgp_delta}, PME={pme_delta}"
    print("\n✅ PGP calculation is working correctly!")

if __name__ == "__main__":
    test_pgp_simple()

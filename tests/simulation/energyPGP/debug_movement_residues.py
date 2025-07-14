"""
Debug movement residues calculation
"""

import pygcmc
from . import pgp_wrapper
from .pgp_wrapper import setPMEParameters, initializePMEParameters, computeMovementEnergyPME
from pygcmc import MCState, MCAtom, MCResidue
from pygcmc import MCForceField, MCMovementResidueInfo

from .test_decorators import pgp_unstable_test

@pgp_unstable_test
def test_movement_residues():
    """Debug movement residues"""
    
    state = MCState()
    state.info.box = [5.0, 5.0, 5.0]
    state.info.cutoff = 2.5
    
    # Force field with LJ
    ff = MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.ljEps = [1.0]
    ff.ljSigma = [0.4]
    state.forcefield = ff
    
    atoms = []
    residues = []
    
    # Three particles
    positions = [[2.0, 2.5, 2.5], [2.5, 2.5, 2.5], [3.0, 2.5, 2.5]]
    
    for i, pos in enumerate(positions):
        atom = MCAtom()
        atom.x, atom.y, atom.z = pos
        atom.charge = 0.0  # No charge
        atom.type = 0
        atoms.append(atom)
        
        res = MCResidue()
        res.active = True
        res.fixed = (i == 0)  # First is fixed, others are moveable
        res.atomStart = i
        res.atomCount = 1
        res.type = 0
        residues.append(res)
    
    state.atoms = atoms
    state.activeAtomCount = 3
    state.residues = residues
    state.activeResidueCount = 3
    
    # Set up movement residues
    movement_info = MCMovementResidueInfo()
    movement_info.startIndex = 1  # Start from second residue
    movement_info.activeCount = 2  # Two moveable residues
    state.movementResidues = [movement_info]
    
    print("\nInitial state:")
    for i, res in enumerate(state.residues):
        print(f"Residue {i}: fixed={res.fixed}, active={res.active}, atomStart={res.atomStart}")
    
    print(f"\nMovement residues: {len(state.movementResidues)}")
    for i, mr in enumerate(state.movementResidues):
        print(f"  Movement group {i}: startIndex={mr.startIndex}, activeCount={mr.activeCount}")
    
    # Check again before calculation
    print(f"\nMovement residues before computeSystemEnergyCutoff: {len(state.movementResidues)}")
    
    # Calculate cutoff energy to verify basic calculation works
    pygcmc.computeSystemEnergyCutoff(state)
    elec_total, vdw_total = pygcmc.getTotalEnergyComponents(state)
    print(f"\nTotal system energy (cutoff):")
    print(f"  VdW: {vdw_total:.6f}")
    
    print(f"\nResidue energies after system calculation:")
    for i, res in enumerate(state.residues):
        print(f"  Residue {i}: vdw={res.energy_vdw:.6f}")
    
    # Reset energies
    for res in state.residues:
        res.energy_vdw = 0.0
        res.energy_elec = 0.0
    
    # Calculate movement energy
    pygcmc.computeMovementEnergyCutoff(state)
    
    print(f"\nResidue energies after movement calculation:")
    for i, res in enumerate(state.residues):
        print(f"  Residue {i}: vdw={res.energy_vdw:.6f}")
    
    # Try PME movement
    alpha = 2.2
    mesh_size = [32, 32, 32]
    spline_order = 4
    
    pgp_wrapper.setPMEParameters(alpha, mesh_size, spline_order)
    pgp_wrapper.initializePMEParameters(state.info.cutoff, state.info.box, alpha)
    
    # Reset energies
    for res in state.residues:
        res.energy_vdw = 0.0
        res.energy_elec = 0.0
    
    # Enable debug output
    pygcmc.setEnergyDebugOutput(True)
    
    print("\n" + "="*50)
    print("Calling computeMovementEnergyPME...")
    print("="*50)
    
    pme_result = pygcmc.computeMovementEnergyPME(state)
    pme_elec = pme_result[0]
    pme_vdw = pme_result[1]
    
    print(f"\nPME movement energy:")
    print(f"  Electrostatic: {pme_elec:.6f}")
    print(f"  VdW: {pme_vdw:.6f}")
    
    print(f"\nResidue energies after PME movement:")
    for i, res in enumerate(state.residues):
        if res.active:
            print(f"  Residue {i}: vdw={res.energy_vdw:.6f}")

if __name__ == "__main__":
    test_movement_residues()

"""
Debug VdW calculation in movement energy
"""

import pygcmc
from pygcmc import MCState, MCAtom, MCResidue, MCForceField, MCMovementResidueInfo


def test_vdw_movement_debug():
    """Debug VdW in movement energy calculation"""
    
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
    
    # Two particles close together
    positions = [[2.0, 2.5, 2.5], [2.5, 2.5, 2.5]]  # 0.5 nm apart
    
    for i, pos in enumerate(positions):
        atom = MCAtom()
        atom.x, atom.y, atom.z = pos
        atom.charge = 0.0  # No charge
        atom.type = 0
        atoms.append(atom)
        
        res = MCResidue()
        res.active = True
        res.fixed = (i == 0)  # First is fixed, second is moveable
        res.atomStart = i
        res.atomCount = 1
        res.type = 0
        residues.append(res)
    
    state.atoms = atoms
    state.activeAtomCount = 2
    state.residues = residues
    state.activeResidueCount = 2
    
    # Set up movement residues
    state.movementResidues = []
    movement_info = MCMovementResidueInfo()
    movement_info.startIndex = 1  # Second residue
    movement_info.activeCount = 1
    state.movementResidues.append(movement_info)
    
    print("\nInitial state:")
    print(f"Residue 0 (fixed): pos=({atoms[0].x}, {atoms[0].y}, {atoms[0].z})")
    print(f"Residue 1 (moveable): pos=({atoms[1].x}, {atoms[1].y}, {atoms[1].z})")
    print(f"Distance: 0.5 nm")
    
    # Calculate cutoff energy
    pygcmc.computeSystemEnergyCutoff(state)
    elec_cutoff, vdw_cutoff = pygcmc.getTotalEnergyComponents(state)
    print(f"\nCutoff energy:")
    print(f"  Electrostatic: {elec_cutoff:.6f} (should be 0)")
    print(f"  VdW: {vdw_cutoff:.6f}")
    
    # Check residue energies after cutoff calculation
    print(f"\nResidue energies after cutoff:")
    for i, res in enumerate(state.residues):
        print(f"  Residue {i}: elec={res.energy_elec:.6f}, vdw={res.energy_vdw:.6f}")
    
    # Reset residue energies
    for res in state.residues:
        res.energy_elec = 0.0
        res.energy_vdw = 0.0
    
    # Calculate movement energy with cutoff
    pygcmc.computeMovementEnergyCutoff(state)
    
    # Check residue energies after movement calculation
    print(f"\nResidue energies after movement energy cutoff:")
    for i, res in enumerate(state.residues):
        print(f"  Residue {i}: elec={res.energy_elec:.6f}, vdw={res.energy_vdw:.6f}")
    
    # Set up PME
    alpha = 2.2
    mesh_size = [32, 32, 32]
    spline_order = 4
    
    pygcmc.setPMEParameters(alpha, mesh_size, spline_order)
    pygcmc.initializePMEParameters(state.info.cutoff, state.info.box, alpha, mesh_size, spline_order)
    
    # Calculate PME movement energy
    pme_result = pygcmc.computeMovementEnergyPME(state)
    pme_elec = pme_result[0]
    pme_vdw = pme_result[1]
    pme_dict = pme_result[2]
    
    print(f"\nPME movement energy:")
    print(f"  Electrostatic: {pme_elec:.6f} (should be 0)")
    print(f"  VdW: {pme_vdw:.6f}")
    print(f"  Dict: {pme_dict}")
    
    # Check residue energies after PME
    print(f"\nResidue energies after PME movement:")
    for i, res in enumerate(state.residues):
        print(f"  Residue {i}: elec={res.energy_elec:.6f}, vdw={res.energy_vdw:.6f}")
    
    # Calculate with PME Complete
    elec_complete, vdw_complete, total_complete = pygcmc.computeSystemEnergyPMEComplete(state)
    print(f"\nPME Complete:")
    print(f"  Electrostatic: {elec_complete:.6f}")
    print(f"  VdW: {vdw_complete:.6f}")
    print(f"  Total: {total_complete:.6f}")


if __name__ == "__main__":
    test_vdw_movement_debug()
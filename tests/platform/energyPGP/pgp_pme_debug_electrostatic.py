"""
Debug PGP vs PME Complete for simple electrostatic systems
"""

import pytest
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pygcmc
from pygcmc import MCState, MCAtom, MCResidue, MCForceField
from pygcmc import initializePMEParameters, computeSystemEnergyPMEComplete
from pygcmc import setPMEParameters, setPGPParameters
from pygcmc import precomputeGridPotential, calculateMoleculeEnergy
from pygcmc import computeMovementEnergyPME


def create_two_particle_system():
    """Create the simplest possible system: one fixed, one moveable particle"""
    state = MCState()
    state.info.box = [10.0, 10.0, 10.0]  # Large box to minimize periodic effects
    state.info.cutoff = 4.0
    
    # Force field with zero LJ
    ff = MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.ljEps = [0.0]
    ff.ljSigma = [0.0]
    state.forcefield = ff
    
    atoms = []
    residues = []
    
    # Fixed particle at origin
    fixed_atom = MCAtom()
    fixed_atom.x, fixed_atom.y, fixed_atom.z = 5.0, 5.0, 5.0
    fixed_atom.charge = 1.0
    fixed_atom.type = 0
    atoms.append(fixed_atom)
    
    fixed_res = MCResidue()
    fixed_res.active = True
    fixed_res.fixed = True
    fixed_res.atomStart = 0
    fixed_res.atomCount = 1
    fixed_res.type = 0
    residues.append(fixed_res)
    
    # Moveable particle
    moveable_atom = MCAtom()
    moveable_atom.x, moveable_atom.y, moveable_atom.z = 6.0, 5.0, 5.0
    moveable_atom.charge = -1.0
    moveable_atom.type = 0
    atoms.append(moveable_atom)
    
    moveable_res = MCResidue()
    moveable_res.active = True
    moveable_res.fixed = False
    moveable_res.atomStart = 1
    moveable_res.atomCount = 1
    moveable_res.type = 0
    residues.append(moveable_res)
    
    state.atoms = atoms
    state.activeAtomCount = 2
    state.residues = residues
    state.activeResidueCount = 2
    
    # Set up movement residues for PME movement energy
    state.movementResidues = []
    movement_info = pygcmc.MCMovementResidueInfo()
    movement_info.startIndex = 1
    movement_info.activeCount = 1
    state.movementResidues.append(movement_info)
    
    return state


def test_simple_two_particle_system():
    """Test the simplest possible case"""
    
    print("\n" + "="*70)
    print("Two Particle System Test")
    print("="*70)
    
    state = create_two_particle_system()
    
    # Use high precision parameters
    alpha = 2.5
    mesh_size = [128, 128, 128]
    spline_order = 6
    
    # Initialize PME
    setPMEParameters(alpha, mesh_size, spline_order)
    initializePMEParameters(state.info.cutoff, state.info.box, alpha, mesh_size, spline_order)
    
    # Initialize PGP with same parameters
    setPGPParameters(alpha, mesh_size, state.info.cutoff, mesh_size, spline_order, 1e-8)
    
    print("\nSystem setup:")
    print(f"  Fixed particle: pos=({state.atoms[0].x},{state.atoms[0].y},{state.atoms[0].z}), charge={state.atoms[0].charge}")
    print(f"  Moveable particle: pos=({state.atoms[1].x},{state.atoms[1].y},{state.atoms[1].z}), charge={state.atoms[1].charge}")
    
    # Calculate initial energies
    elec_init, vdw_init, total_init = computeSystemEnergyPMEComplete(state)
    print(f"\nInitial PME Complete: {total_init:.10f} kJ/mol")
    
    # Calculate PME movement energy
    pme_move_result = computeMovementEnergyPME(state)
    pme_move_elec = pme_move_result[0]
    pme_move_dict = pme_move_result[2]
    print(f"PME movement energy: {pme_move_elec:.10f} kJ/mol")
    print(f"  Components: {pme_move_dict}")
    
    # Precompute PGP grid
    precomputeGridPotential(state, fixed_only=True)
    pgp_init = calculateMoleculeEnergy(state)
    print(f"Initial PGP: {pgp_init:.10f} kJ/mol")
    
    # Compare initial state
    print(f"\nInitial comparison:")
    print(f"  PME movement reciprocal: {pme_move_dict.get('reciprocal', 0):.10f}")
    print(f"  PGP energy:             {pgp_init:.10f}")
    print(f"  Difference:             {abs(pme_move_dict.get('reciprocal', 0) - pgp_init):.10f}")
    
    # Move particle by small amount
    dx = 0.1
    state.atoms[1].x += dx
    
    print(f"\nMoved particle by {dx} nm in x direction")
    
    # Calculate new energies
    elec_final, vdw_final, total_final = computeSystemEnergyPMEComplete(state)
    pgp_final = calculateMoleculeEnergy(state)
    
    # Calculate PME movement energy after move
    pme_move_result_final = computeMovementEnergyPME(state)
    pme_move_final = pme_move_result_final[0]
    pme_move_dict_final = pme_move_result_final[2]
    
    # Calculate dE
    pme_de = total_final - total_init
    pgp_de = pgp_final - pgp_init
    pme_move_de = pme_move_final - pme_move_elec
    
    print(f"\nEnergy changes:")
    print(f"  PME Complete dE:  {pme_de:.10f} kJ/mol")
    print(f"  PGP dE:          {pgp_de:.10f} kJ/mol")
    print(f"  PME movement dE: {pme_move_de:.10f} kJ/mol")
    
    print(f"\nRelative errors:")
    rel_err_pgp = abs(pgp_de - pme_de) / abs(pme_de) * 100
    rel_err_move = abs(pme_move_de - pme_de) / abs(pme_de) * 100
    print(f"  PGP vs PME Complete:      {rel_err_pgp:.6f}%")
    print(f"  PME move vs PME Complete: {rel_err_move:.6f}%")
    
    # Check what PGP is actually calculating
    print(f"\nDiagnostic:")
    print(f"  PME reciprocal change: {pme_move_dict_final.get('reciprocal', 0) - pme_move_dict.get('reciprocal', 0):.10f}")
    print(f"  PGP change:           {pgp_de:.10f}")
    
    # Assert that the relative error is reasonable
    assert rel_err_pgp < 50.0, f"Relative error too large: {rel_err_pgp:.2f}%"


def test_pgp_self_consistency():
    """Test PGP self-consistency by comparing with known analytical result"""
    
    print("\n" + "="*70)
    print("PGP Self-Consistency Test")
    print("="*70)
    
    # Create a very simple system
    state = MCState()
    state.info.box = [20.0, 20.0, 20.0]  # Very large box
    state.info.cutoff = 8.0
    
    # Force field with zero LJ
    ff = MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.ljEps = [0.0]
    ff.ljSigma = [0.0]
    state.forcefield = ff
    
    atoms = []
    residues = []
    
    # Just two particles far apart
    # Fixed particle
    atom1 = MCAtom()
    atom1.x, atom1.y, atom1.z = 10.0, 10.0, 10.0
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
    
    # Moveable particle far away
    atom2 = MCAtom()
    atom2.x, atom2.y, atom2.z = 15.0, 10.0, 10.0  # 5 nm away
    atom2.charge = 1.0
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
    
    # Use parameters that should give good convergence
    alpha = 3.0
    mesh_size = [128, 128, 128]
    spline_order = 6
    
    # Initialize
    setPMEParameters(alpha, mesh_size, spline_order)
    initializePMEParameters(state.info.cutoff, state.info.box, alpha, mesh_size, spline_order)
    setPGPParameters(alpha, mesh_size, state.info.cutoff, mesh_size, spline_order, 1e-8)
    
    # Precompute grid
    precomputeGridPotential(state, fixed_only=True)
    
    # Calculate energy at different distances
    distances = [5.0, 4.0, 3.0, 2.0]
    
    print(f"\n{'Distance':>10} | {'PGP Energy':>15} | {'PME Energy':>15} | {'Analytical':>15}")
    print("-" * 60)
    
    for d in distances:
        # Set position
        state.atoms[1].x = 10.0 + d
        
        # Calculate energies
        pgp_energy = calculateMoleculeEnergy(state)
        elec, _, total = computeSystemEnergyPMEComplete(state)
        
        # Analytical coulomb energy (with unit conversion)
        # E = k * q1 * q2 / r where k = 138.935456 kJ*nm/(mol*e^2)
        k_coulomb = 138.935456
        analytical = k_coulomb * 1.0 * 1.0 / d
        
        print(f"{d:10.1f} | {pgp_energy:15.8f} | {total:15.8f} | {analytical:15.8f}")
    
    # The trend should be consistent even if absolute values differ


if __name__ == "__main__":
    test_simple_two_particle_system()
    test_pgp_self_consistency()
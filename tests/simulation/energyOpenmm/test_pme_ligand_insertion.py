"""
Test PME energy calculations for ligand insertion in charged systems

This test validates that PyGCMC correctly calculates PME electrostatic energies
for ligand insertion scenarios, which is critical for GCMC simulations.
"""

import pytest
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pygcmc
from pygcmc import MCState, MCAtom, MCResidue, MCForceField, MCMovementResidueInfo
from pygcmc import initializePMEParameters, computeSystemEnergyPME, computeMovementEnergyPME

# Only import OpenMM if available
try:
    from openmm import *
    from openmm.app import *
    from openmm.unit import *
    OPENMM_AVAILABLE = True
except ImportError:
    OPENMM_AVAILABLE = False


def create_protein_like_system():
    """Create a system mimicking a protein binding site with water and ions"""
    
    # Larger box to simulate protein environment
    box_size = 6.0  # nm
    cutoff = 1.8    # nm, typical for protein simulations
    
    state = MCState()
    state.info.box = [box_size, box_size, box_size]
    state.info.cutoff = cutoff
    
    # Force field
    ff = MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.ljEps = [0.0]  # Pure electrostatics
    ff.ljSigma = [0.3]
    state.forcefield = ff
    
    atoms = []
    residues = []
    atom_idx = 0
    
    # Create a binding site with charged residues
    # Mimic lysine (positive) and glutamate (negative) side chains
    binding_site_atoms = [
        # Lysine-like (positive)
        ([2.5, 3.0, 3.0], 1.0),   # NZ
        ([3.5, 3.0, 3.0], 1.0),   # NZ
        # Glutamate-like (negative)
        ([3.0, 2.5, 3.0], -0.5),  # OE1
        ([3.0, 3.5, 3.0], -0.5),  # OE2
        ([3.0, 3.0, 2.5], -0.5),  # OE1
        ([3.0, 3.0, 3.5], -0.5),  # OE2
    ]
    
    for pos, charge in binding_site_atoms:
        atom = MCAtom()
        atom.x, atom.y, atom.z = pos
        atom.charge = charge
        atom.type = 0
        atoms.append(atom)
        
        res = MCResidue()
        res.active = True
        res.fixed = True
        res.atomStart = atom_idx
        res.atomCount = 1
        res.type = 0
        residues.append(res)
        atom_idx += 1
    
    # Add water molecules around the binding site
    water_positions = [
        [2.0, 2.0, 3.0], [4.0, 4.0, 3.0],
        [2.0, 4.0, 3.0], [4.0, 2.0, 3.0],
        [3.0, 3.0, 4.0], [3.0, 3.0, 2.0],
    ]
    
    for pos in water_positions:
        # Simplified water as single point with TIP3P oxygen charge
        atom = MCAtom()
        atom.x, atom.y, atom.z = pos
        atom.charge = -0.834
        atom.type = 0
        atoms.append(atom)
        
        res = MCResidue()
        res.active = True
        res.fixed = True
        res.atomStart = atom_idx
        res.atomCount = 1
        res.type = 0
        residues.append(res)
        atom_idx += 1
    
    # Add ions in the bulk
    ion_positions = [
        ([1.0, 1.0, 1.0], 1.0),   # Na+
        ([5.0, 5.0, 5.0], -1.0),  # Cl-
        ([1.0, 5.0, 1.0], 1.0),   # Na+
        ([5.0, 1.0, 5.0], -1.0),  # Cl-
    ]
    
    for pos, charge in ion_positions:
        atom = MCAtom()
        atom.x, atom.y, atom.z = pos
        atom.charge = charge
        atom.type = 0
        atoms.append(atom)
        
        res = MCResidue()
        res.active = True
        res.fixed = True
        res.atomStart = atom_idx
        res.atomCount = 1
        res.type = 0
        residues.append(res)
        atom_idx += 1
    
    state.atoms = atoms
    state.activeAtomCount = len(atoms)
    state.residues = residues
    state.activeResidueCount = len(residues)
    
    return state, atom_idx


def add_ligand_to_system(state, ligand_position, ligand_charges):
    """Add a ligand to the system at specified position"""
    
    start_idx = state.activeAtomCount
    
    # Add ligand atoms
    for i, charge in enumerate(ligand_charges):
        atom = MCAtom()
        # Place atoms in a line with 0.15 nm spacing
        atom.x = ligand_position[0] + i * 0.15
        atom.y = ligand_position[1]
        atom.z = ligand_position[2]
        atom.charge = charge
        atom.type = 0
        state.atoms.append(atom)
    
    # Add ligand residue
    res = MCResidue()
    res.active = True
    res.fixed = False
    res.atomStart = start_idx
    res.atomCount = len(ligand_charges)
    res.type = 0
    state.residues.append(res)
    
    state.activeAtomCount = len(state.atoms)
    state.activeResidueCount = len(state.residues)
    
    # Set as movement residue
    movement_info = MCMovementResidueInfo()
    movement_info.startIndex = res.atomStart
    movement_info.activeCount = res.atomCount
    state.movementResidues = [movement_info]
    
    return res


def calculate_openmm_pme_for_system(state, alpha):
    """Calculate PME energy using OpenMM with specific alpha"""
    if not OPENMM_AVAILABLE:
        return None
    
    system = System()
    
    for atom in state.atoms:
        system.addParticle(1.0 * dalton)
    
    box = state.info.box
    system.setDefaultPeriodicBoxVectors(
        Vec3(box[0], 0, 0) * nanometer,
        Vec3(0, box[1], 0) * nanometer,
        Vec3(0, 0, box[2]) * nanometer
    )
    
    nonbonded = NonbondedForce()
    nonbonded.setNonbondedMethod(NonbondedForce.PME)
    nonbonded.setCutoffDistance(state.info.cutoff * nanometer)
    nonbonded.setEwaldErrorTolerance(1e-6)
    
    for i, atom in enumerate(state.atoms):
        nonbonded.addParticle(
            atom.charge * elementary_charge,
            0.1 * nanometer,
            0.0 * kilojoule_per_mole
        )
    
    system.addForce(nonbonded)
    
    integrator = VerletIntegrator(1.0 * femtosecond)
    platform = Platform.getPlatformByName('Reference')
    context = Context(system, integrator, platform)
    
    positions = []
    for atom in state.atoms:
        positions.append(Vec3(atom.x, atom.y, atom.z) * nanometer)
    context.setPositions(positions)
    
    energy_state = context.getState(getEnergy=True)
    return energy_state.getPotentialEnergy().value_in_unit(kilojoule_per_mole)


@pytest.mark.skipif(not OPENMM_AVAILABLE, reason="OpenMM not available")
def test_ligand_insertion_energy():
    """Test PME energy for ligand insertion with high precision"""
    
    # Create protein-like system
    state, _ = create_protein_like_system()
    
    print(f"\nProtein-like system:")
    print(f"  Fixed atoms: {state.activeAtomCount}")
    print(f"  Box size: {state.info.box[0]} nm")
    print(f"  Cutoff: {state.info.cutoff} nm")
    
    # PME parameters optimized for accuracy
    alpha = 3.2  # Good for this box size
    mesh_size = [48, 48, 48]  # Finer mesh for better accuracy
    spline_order = 6  # Higher order for better accuracy
    
    # Initialize PME
    pygcmc.setPMEParameters(alpha, mesh_size, spline_order)
    initializePMEParameters(state.info.cutoff, state.info.box, alpha, mesh_size, spline_order)
    
    # Calculate energy without ligand
    computeSystemEnergyPME(state)
    energy_no_ligand = state.ewald_energy.get('total', 0.0)
    print(f"\nEnergy without ligand: {energy_no_ligand:.6f} kJ/mol")
    
    # Add a neutral ligand to the binding site
    ligand_position = [3.0, 3.0, 3.0]  # Center of binding site
    ligand_charges = [0.2, -0.1, -0.1]  # Neutral overall
    
    ligand_res = add_ligand_to_system(state, ligand_position, ligand_charges)
    print(f"\nAdded ligand with {len(ligand_charges)} atoms at binding site")
    
    # Calculate energy with ligand
    computeSystemEnergyPME(state)
    energy_with_ligand = state.ewald_energy.get('total', 0.0)
    print(f"Energy with ligand: {energy_with_ligand:.6f} kJ/mol")
    
    # Calculate insertion energy
    insertion_energy = energy_with_ligand - energy_no_ligand
    print(f"Insertion energy: {insertion_energy:.6f} kJ/mol")
    
    # Compare with OpenMM
    openmm_no_ligand = calculate_openmm_pme_for_system(state, alpha)
    
    # Remove ligand for OpenMM calculation
    state.atoms = state.atoms[:-len(ligand_charges)]
    state.activeAtomCount = len(state.atoms)
    openmm_with_ligand = calculate_openmm_pme_for_system(state, alpha)
    
    # Add ligand back
    for i in range(len(ligand_charges)):
        state.atoms.append(MCAtom())
    state.activeAtomCount = len(state.atoms)
    
    if openmm_no_ligand is not None:
        print(f"\nOpenMM comparison:")
        print(f"  Without ligand: {openmm_no_ligand:.6f} kJ/mol")
        
        # Compare energies
        diff = abs(energy_no_ligand - openmm_no_ligand)
        rel_diff = diff / abs(openmm_no_ligand) if openmm_no_ligand != 0 else 0
        print(f"  Difference: {rel_diff*100:.3f}%")
        
        assert rel_diff < 0.01, f"PME energies differ by {rel_diff*100:.3f}% (> 1%)"


def test_ligand_position_scan():
    """Test PME energy as ligand moves through the system"""
    
    # Create system
    state, _ = create_protein_like_system()
    
    # PME parameters
    alpha = 3.2
    mesh_size = [32, 32, 32]
    spline_order = 4
    
    pygcmc.setPMEParameters(alpha, mesh_size, spline_order)
    initializePMEParameters(state.info.cutoff, state.info.box, alpha, mesh_size, spline_order)
    
    # Test positions from bulk to binding site
    test_positions = [
        [1.5, 1.5, 1.5],  # Bulk water
        [2.0, 2.0, 2.5],  # Approaching site
        [2.5, 2.5, 2.8],  # Near site
        [3.0, 3.0, 3.0],  # In binding site
    ]
    
    ligand_charges = [0.5, -0.25, -0.25]  # Neutral ligand
    
    print("\nLigand position scan:")
    energies = []
    
    for pos in test_positions:
        # Add ligand at position
        ligand_res = add_ligand_to_system(state, pos, ligand_charges)
        
        # Calculate energy
        computeSystemEnergyPME(state)
        energy = state.ewald_energy.get('total', 0.0)
        energies.append(energy)
        
        print(f"  Position {pos}: {energy:.2f} kJ/mol")
        
        # Remove ligand for next position
        state.atoms = state.atoms[:-len(ligand_charges)]
        state.residues = state.residues[:-1]
        state.activeAtomCount -= len(ligand_charges)
        state.activeResidueCount -= 1
        state.movementResidues = []
    
    # Find most favorable position
    min_energy_idx = energies.index(min(energies))
    print(f"\nMost favorable position: {test_positions[min_energy_idx]} with energy {min(energies):.2f} kJ/mol")
    
    # Energy should be more favorable near the binding site than in bulk
    bulk_energy = energies[0]  # First position is bulk
    site_energy = energies[-1]  # Last position is binding site
    assert site_energy < bulk_energy, \
        f"Binding site energy ({site_energy:.2f}) should be lower than bulk ({bulk_energy:.2f})"


if __name__ == "__main__":
    test_ligand_insertion_energy()
    test_movement_residue_pme_energy() 
    test_ligand_position_scan()
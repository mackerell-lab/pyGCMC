"""
Analyze LJ differences in multi-particle systems
"""

import pytest
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pygcmc
from pygcmc import MCState, MCAtom, MCResidue, MCForceField
from pygcmc import computeSystemEnergyCutoffFixed

# Only import OpenMM if available
try:
    from openmm import *
    from openmm.app import *
    from openmm.unit import *
    OPENMM_AVAILABLE = True
except ImportError:
    OPENMM_AVAILABLE = False


def create_test_system_progressive(n_atoms):
    """Create a system with n atoms in a grid pattern"""
    
    state = MCState()
    state.info.box = [5.0, 5.0, 5.0]  # Smaller box to have more interactions
    state.info.cutoff = 2.0  # Larger cutoff
    
    ff = MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.ljEps = [1.0]
    ff.ljSigma = [0.35]
    state.forcefield = ff
    
    atoms = []
    residues = []
    
    # Place atoms in a grid
    spacing = 0.6  # nm
    for i in range(n_atoms):
        atom = MCAtom()
        # Simple linear arrangement
        atom.x = 2.5 + (i % 3) * spacing
        atom.y = 2.5 + (i // 3) * spacing
        atom.z = 2.5
        atom.charge = 0.0
        atom.type = 0
        atoms.append(atom)
        
        res = MCResidue()
        res.active = True
        res.fixed = True
        res.atomStart = i
        res.atomCount = 1
        res.type = 0
        residues.append(res)
    
    state.atoms = atoms
    state.activeAtomCount = n_atoms
    state.residues = residues
    state.activeResidueCount = n_atoms
    
    return state


def calculate_openmm_energy(state):
    """Calculate energy using OpenMM"""
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
    nonbonded.setNonbondedMethod(NonbondedForce.CutoffPeriodic)
    nonbonded.setCutoffDistance(state.info.cutoff * nanometer)
    
    ff = state.forcefield
    for atom in state.atoms:
        nonbonded.addParticle(
            0.0 * elementary_charge,
            ff.ljSigma[0] * nanometer,
            ff.ljEps[0] * kilojoule_per_mole
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


def count_interactions(atoms, box, cutoff):
    """Count number of interactions within cutoff"""
    n_atoms = len(atoms)
    count = 0
    
    for i in range(n_atoms):
        for j in range(i+1, n_atoms):
            dx = atoms[i].x - atoms[j].x
            dy = atoms[i].y - atoms[j].y
            dz = atoms[i].z - atoms[j].z
            
            # Apply minimum image convention
            dx -= box[0] * round(dx / box[0])
            dy -= box[1] * round(dy / box[1])
            dz -= box[2] * round(dz / box[2])
            
            r = math.sqrt(dx*dx + dy*dy + dz*dz)
            
            if r < cutoff:
                count += 1
                
    return count


@pytest.mark.skipif(not OPENMM_AVAILABLE, reason="OpenMM not available")
def test_progressive_atoms():
    """Test systems with increasing number of atoms"""
    
    print("\nProgressive atom test:")
    print("N_atoms | PyGCMC    | OpenMM    | Diff(%)  | N_pairs")
    print("-" * 50)
    
    for n in [2, 3, 4, 5, 6]:
        state = create_test_system_progressive(n)
        
        # PyGCMC calculation
        elec, vdw, total = computeSystemEnergyCutoffFixed(state)
        
        # OpenMM calculation
        openmm_energy = calculate_openmm_energy(state)
        
        # Count interactions
        n_pairs = count_interactions(state.atoms, state.info.box, state.info.cutoff)
        
        if openmm_energy is not None:
            diff = abs(vdw - openmm_energy)
            rel_diff = diff / abs(openmm_energy) * 100 if openmm_energy != 0 else 0
            
            print(f"{n:7d} | {vdw:9.4f} | {openmm_energy:9.4f} | {rel_diff:8.2f} | {n_pairs:7d}")


@pytest.mark.skipif(not OPENMM_AVAILABLE, reason="OpenMM not available")
def test_cutoff_effects():
    """Test effect of different cutoffs"""
    
    print("\n\nCutoff effect test (6 atoms):")
    print("Cutoff  | PyGCMC    | OpenMM    | Diff(%)  | N_pairs")
    print("-" * 50)
    
    state = create_test_system_progressive(6)
    
    for cutoff in [0.8, 1.0, 1.2, 1.5, 2.0, 2.5]:
        state.info.cutoff = cutoff
        
        # PyGCMC calculation
        elec, vdw, total = computeSystemEnergyCutoffFixed(state)
        
        # OpenMM calculation
        openmm_energy = calculate_openmm_energy(state)
        
        # Count interactions
        n_pairs = count_interactions(state.atoms, state.info.box, state.info.cutoff)
        
        if openmm_energy is not None:
            diff = abs(vdw - openmm_energy)
            rel_diff = diff / abs(openmm_energy) * 100 if openmm_energy != 0 else 0
            
            print(f"{cutoff:7.1f} | {vdw:9.4f} | {openmm_energy:9.4f} | {rel_diff:8.2f} | {n_pairs:7d}")


@pytest.mark.skipif(not OPENMM_AVAILABLE, reason="OpenMM not available")
def test_pbc_boundary():
    """Test atoms near PBC boundary"""
    
    print("\n\nPBC boundary test:")
    
    state = MCState()
    state.info.box = [3.0, 3.0, 3.0]
    state.info.cutoff = 1.5
    
    ff = MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.ljEps = [1.0]
    ff.ljSigma = [0.35]
    state.forcefield = ff
    
    # Place atoms near box boundaries
    positions = [
        [0.1, 1.5, 1.5],   # Near x=0 boundary
        [2.9, 1.5, 1.5],   # Near x=box boundary
        [1.5, 0.1, 1.5],   # Near y=0 boundary
        [1.5, 2.9, 1.5],   # Near y=box boundary
    ]
    
    atoms = []
    residues = []
    
    for i, pos in enumerate(positions):
        atom = MCAtom()
        atom.x, atom.y, atom.z = pos
        atom.charge = 0.0
        atom.type = 0
        atoms.append(atom)
        
        res = MCResidue()
        res.active = True
        res.fixed = True
        res.atomStart = i
        res.atomCount = 1
        res.type = 0
        residues.append(res)
    
    state.atoms = atoms
    state.activeAtomCount = len(atoms)
    state.residues = residues
    state.activeResidueCount = len(residues)
    
    # Calculate energies
    elec, vdw, total = computeSystemEnergyCutoffFixed(state)
    openmm_energy = calculate_openmm_energy(state)
    
    print(f"  PyGCMC: {vdw:.6f} kJ/mol")
    print(f"  OpenMM: {openmm_energy:.6f} kJ/mol")
    
    if openmm_energy is not None:
        diff = abs(vdw - openmm_energy)
        rel_diff = diff / abs(openmm_energy) * 100 if openmm_energy != 0 else 0
        print(f"  Difference: {rel_diff:.2f}%")
        
    # Check distances with PBC
    print("\n  Pairwise distances (with PBC):")
    for i in range(len(atoms)):
        for j in range(i+1, len(atoms)):
            dx = atoms[i].x - atoms[j].x
            dy = atoms[i].y - atoms[j].y
            dz = atoms[i].z - atoms[j].z
            
            # Apply minimum image convention
            dx -= state.info.box[0] * round(dx / state.info.box[0])
            dy -= state.info.box[1] * round(dy / state.info.box[1])
            dz -= state.info.box[2] * round(dz / state.info.box[2])
            
            r = math.sqrt(dx*dx + dy*dy + dz*dz)
            print(f"    {i}-{j}: {r:.3f} nm")


if __name__ == "__main__":
    test_progressive_atoms()
    test_cutoff_effects()
    test_pbc_boundary()
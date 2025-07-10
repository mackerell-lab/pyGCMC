"""
Test to verify PyGCMC's handling of intramolecular LJ interactions

This test creates molecules with multiple atoms to check if PyGCMC
includes intramolecular (within-residue) LJ interactions.
"""

import math
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pygcmc
from pygcmc import MCState, MCAtom, MCResidue, MCForceField
from pygcmc import computeSystemEnergyCutoff

# Only import OpenMM if available
try:
    from openmm import *
    from openmm.app import *
    from openmm.unit import *
    OPENMM_AVAILABLE = True
except ImportError:
    OPENMM_AVAILABLE = False


def calculate_manual_lj(r, epsilon, sigma):
    """Calculate LJ energy manually"""
    sr6 = (sigma/r)**6
    return 4.0 * epsilon * (sr6*sr6 - sr6)


def test_single_molecule_lj():
    """Test LJ energy for a single molecule with multiple atoms"""
    
    print("\nSingle molecule LJ test:")
    print("="*60)
    
    state = MCState()
    state.info.box = [10.0, 10.0, 10.0]  # Large box to avoid PBC
    state.info.cutoff = 5.0
    
    # Simple force field
    ff = MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.ljEps = [1.0]     # kJ/mol
    ff.ljSigma = [0.35]  # nm
    state.forcefield = ff
    
    # Create a single residue with 3 atoms in a triangle
    atoms = []
    positions = [
        [5.0, 5.0, 5.0],
        [5.5, 5.0, 5.0],    # 0.5 nm from atom 0
        [5.25, 5.433, 5.0]  # Forms equilateral triangle
    ]
    
    for pos in positions:
        atom = MCAtom()
        atom.x, atom.y, atom.z = pos
        atom.charge = 0.0
        atom.type = 0
        atoms.append(atom)
    
    # Single residue containing all 3 atoms
    res = MCResidue()
    res.active = True
    res.fixed = True
    res.atomStart = 0
    res.atomCount = 3
    res.type = 0
    
    state.atoms = atoms
    state.activeAtomCount = 3
    state.residues = [res]
    state.activeResidueCount = 1
    
    # Calculate energy with PyGCMC
    computeSystemEnergyCutoff(state)
    pygcmc_lj = state.residues[0].energy_vdw
    
    # Calculate expected energy manually
    # If PyGCMC includes intramolecular interactions:
    # - Pair 0-1: r=0.5 nm
    # - Pair 0-2: r=0.5 nm (equilateral triangle)
    # - Pair 1-2: r=0.5 nm
    expected_with_intra = 0.0
    for i in range(3):
        for j in range(i+1, 3):
            dx = atoms[i].x - atoms[j].x
            dy = atoms[i].y - atoms[j].y
            dz = atoms[i].z - atoms[j].z
            r = math.sqrt(dx*dx + dy*dy + dz*dz)
            expected_with_intra += calculate_manual_lj(r, 1.0, 0.35)
    
    print(f"  Single residue with 3 atoms")
    print(f"  PyGCMC LJ energy: {pygcmc_lj:.6f} kJ/mol")
    print(f"  Expected with intramolecular: {expected_with_intra:.6f} kJ/mol")
    print(f"  Expected without intramolecular: 0.0 kJ/mol")
    
    if abs(pygcmc_lj) < 1e-6:
        print("\n  ✓ PyGCMC DOES NOT include intramolecular LJ interactions!")
        print("    Energy is 0 because all atoms are in the same residue")
    else:
        print("\n  ⚠️  PyGCMC includes intramolecular LJ interactions")
    
    # Test with OpenMM if available
    if OPENMM_AVAILABLE:
        system = System()
        
        for atom in atoms:
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
        
        # Add particles
        for atom in atoms:
            nonbonded.addParticle(
                0.0 * elementary_charge,
                0.35 * nanometer,
                1.0 * kilojoule_per_mole
            )
        
        system.addForce(nonbonded)
        
        integrator = VerletIntegrator(1.0 * femtosecond)
        platform = Platform.getPlatformByName('Reference')
        context = Context(system, integrator, platform)
        
        positions_nm = []
        for atom in atoms:
            positions_nm.append(Vec3(atom.x, atom.y, atom.z) * nanometer)
        context.setPositions(positions_nm)
        
        energy_state = context.getState(getEnergy=True)
        openmm_lj = energy_state.getPotentialEnergy().value_in_unit(kilojoule_per_mole)
        
        print(f"\n  OpenMM LJ energy: {openmm_lj:.6f} kJ/mol")
        print(f"  OpenMM includes all non-bonded interactions by default")


def test_two_molecules_lj():
    """Test LJ energy for two separate molecules"""
    
    print("\n\nTwo molecule LJ test:")
    print("="*60)
    
    state = MCState()
    state.info.box = [10.0, 10.0, 10.0]
    state.info.cutoff = 5.0
    
    ff = MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.ljEps = [1.0]
    ff.ljSigma = [0.35]
    state.forcefield = ff
    
    # Create two residues, each with 2 atoms
    atoms = []
    
    # Residue 1: atoms at (5,5,5) and (5.4,5,5)
    # Residue 2: atoms at (6,5,5) and (6.4,5,5)
    positions = [
        [5.0, 5.0, 5.0],    # Res 1, atom 0
        [5.4, 5.0, 5.0],    # Res 1, atom 1 (0.4 nm intra)
        [6.0, 5.0, 5.0],    # Res 2, atom 0 (0.6 nm from res1-atom1)
        [6.4, 5.0, 5.0],    # Res 2, atom 1 (0.4 nm intra)
    ]
    
    for pos in positions:
        atom = MCAtom()
        atom.x, atom.y, atom.z = pos
        atom.charge = 0.0
        atom.type = 0
        atoms.append(atom)
    
    # Two residues
    res1 = MCResidue()
    res1.active = True
    res1.fixed = True
    res1.atomStart = 0
    res1.atomCount = 2
    res1.type = 0
    
    res2 = MCResidue()
    res2.active = True
    res2.fixed = True
    res2.atomStart = 2
    res2.atomCount = 2
    res2.type = 0
    
    state.atoms = atoms
    state.activeAtomCount = 4
    state.residues = [res1, res2]
    state.activeResidueCount = 2
    
    # Calculate with PyGCMC
    computeSystemEnergyCutoff(state)
    
    res1_energy = state.residues[0].energy_vdw
    res2_energy = state.residues[1].energy_vdw
    total_pygcmc = res1_energy + res2_energy
    
    print(f"  Residue 1 energy: {res1_energy:.6f} kJ/mol")
    print(f"  Residue 2 energy: {res2_energy:.6f} kJ/mol")
    print(f"  Total PyGCMC: {total_pygcmc:.6f} kJ/mol")
    
    # Calculate expected energies
    print("\n  Expected interactions:")
    
    # Intramolecular (within residues) - PyGCMC skips these
    intra1 = calculate_manual_lj(0.4, 1.0, 0.35)  # Res1 internal
    intra2 = calculate_manual_lj(0.4, 1.0, 0.35)  # Res2 internal
    print(f"    Intramolecular Res1 (0-1): {intra1:.6f} kJ/mol")
    print(f"    Intramolecular Res2 (2-3): {intra2:.6f} kJ/mol")
    
    # Intermolecular (between residues) - PyGCMC includes these
    inter_pairs = [
        (0, 2, 1.0),   # Res1-atom0 to Res2-atom0
        (0, 3, 1.4),   # Res1-atom0 to Res2-atom1
        (1, 2, 0.6),   # Res1-atom1 to Res2-atom0
        (1, 3, 1.0),   # Res1-atom1 to Res2-atom1
    ]
    
    total_inter = 0.0
    print("\n  Intermolecular interactions:")
    for i, j, r in inter_pairs:
        energy = calculate_manual_lj(r, 1.0, 0.35)
        total_inter += energy
        print(f"    Atom {i}-{j} (r={r:.1f} nm): {energy:.6f} kJ/mol")
    
    print(f"\n  Total intermolecular only: {total_inter:.6f} kJ/mol")
    print(f"  Total with intramolecular: {total_inter + intra1 + intra2:.6f} kJ/mol")
    
    # Check if PyGCMC matches intermolecular only
    diff = abs(total_pygcmc - total_inter)
    if diff < 1e-3:
        print("\n  ✓ PyGCMC matches intermolecular-only calculation!")
        print("    Confirms: PyGCMC skips intramolecular LJ interactions")
    else:
        print(f"\n  ⚠️  Difference: {diff:.6f} kJ/mol")


if __name__ == "__main__":
    test_single_molecule_lj()
    test_two_molecules_lj()
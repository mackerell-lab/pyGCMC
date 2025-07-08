"""
Test LJ at very close distances to understand the issue
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


def create_two_atom_system_at_distance(distance):
    """Create two atoms at specified distance"""
    state = MCState()
    state.info.box = [10.0, 10.0, 10.0]
    state.info.cutoff = 5.0
    
    ff = MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.ljEps = [1.0]     # kJ/mol
    ff.ljSigma = [0.35]  # nm
    state.forcefield = ff
    
    atoms = []
    residues = []
    
    # Atom 1 at center
    atom1 = MCAtom()
    atom1.x, atom1.y, atom1.z = 5.0, 5.0, 5.0
    atom1.charge = 0.0
    atom1.type = 0
    atoms.append(atom1)
    
    # Atom 2 at specified distance
    atom2 = MCAtom()
    atom2.x = 5.0 + distance
    atom2.y, atom2.z = 5.0, 5.0
    atom2.charge = 0.0
    atom2.type = 0
    atoms.append(atom2)
    
    for i in range(2):
        res = MCResidue()
        res.active = True
        res.fixed = True
        res.atomStart = i
        res.atomCount = 1
        res.type = 0
        residues.append(res)
    
    state.atoms = atoms
    state.activeAtomCount = 2
    state.residues = residues
    state.activeResidueCount = 2
    
    return state


def calculate_lj_energy(r, sigma, epsilon):
    """Calculate LJ energy at distance r"""
    sr6 = (sigma/r)**6
    return 4.0 * epsilon * (sr6*sr6 - sr6)


@pytest.mark.skipif(not OPENMM_AVAILABLE, reason="OpenMM not available")
def test_close_distances():
    """Test LJ at various close distances"""
    
    print("\nLJ energy at close distances:")
    print("Distance | Expected   | PyGCMC    | OpenMM    | PyG-Exp(%) | PyG-OMM(%)")
    print("-" * 75)
    
    sigma = 0.35
    epsilon = 1.0
    
    # Test distances from very close to beyond sigma
    distances = [0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50]
    
    for dist in distances:
        state = create_two_atom_system_at_distance(dist)
        
        # Expected energy
        expected = calculate_lj_energy(dist, sigma, epsilon)
        
        # PyGCMC calculation
        elec, vdw, total = computeSystemEnergyCutoffFixed(state)
        
        # OpenMM calculation
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
        
        for atom in state.atoms:
            nonbonded.addParticle(
                0.0 * elementary_charge,
                sigma * nanometer,
                epsilon * kilojoule_per_mole
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
        openmm_energy = energy_state.getPotentialEnergy().value_in_unit(kilojoule_per_mole)
        
        # Calculate differences
        diff_exp = abs(vdw - expected) / abs(expected) * 100 if expected != 0 else 0
        diff_omm = abs(vdw - openmm_energy) / abs(openmm_energy) * 100 if openmm_energy != 0 else 0
        
        print(f"{dist:8.2f} | {expected:10.2f} | {vdw:10.2f} | {openmm_energy:10.2f} | {diff_exp:10.2f} | {diff_omm:10.2f}")
        
        # Check if PyGCMC caps the energy
        if dist < 0.3 and abs(vdw) < abs(expected) * 0.1:
            print(f"         WARNING: PyGCMC energy much smaller than expected at r={dist}")


def test_energy_capping():
    """Test if PyGCMC caps LJ energy at close distances"""
    
    print("\n\nChecking for energy capping:")
    
    state = MCState()
    state.info.box = [10.0, 10.0, 10.0]
    state.info.cutoff = 5.0
    
    ff = MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.ljEps = [1.0]
    ff.ljSigma = [0.35]
    state.forcefield = ff
    
    # Very close atoms
    atoms = []
    residues = []
    
    # Multiple pairs at different close distances
    pairs = [
        ([5.0, 5.0, 5.0], [5.1, 5.0, 5.0]),   # 0.1 nm
        ([6.0, 5.0, 5.0], [6.15, 5.0, 5.0]),  # 0.15 nm
        ([7.0, 5.0, 5.0], [7.2, 5.0, 5.0]),   # 0.2 nm
    ]
    
    atom_idx = 0
    for pair in pairs:
        for pos in pair:
            atom = MCAtom()
            atom.x, atom.y, atom.z = pos
            atom.charge = 0.0
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
    
    # Calculate energy
    elec, vdw, total = computeSystemEnergyCutoffFixed(state)
    
    print(f"  Total VDW energy for 3 close pairs: {vdw:.6f} kJ/mol")
    
    # Calculate expected
    total_expected = 0.0
    for i in [0, 2, 4]:
        dx = atoms[i].x - atoms[i+1].x
        r = abs(dx)
        expected = calculate_lj_energy(r, 0.35, 1.0)
        total_expected += expected
        print(f"  Pair {i//2+1} at r={r:.2f}: expected {expected:.2f} kJ/mol")
    
    print(f"  Total expected: {total_expected:.2f} kJ/mol")
    
    if abs(vdw) < abs(total_expected) * 0.1:
        print("  ⚠️  PyGCMC appears to cap LJ energy at close distances!")


if __name__ == "__main__":
    test_close_distances()
    test_energy_capping()
"""
Detailed comparison to understand remaining LJ energy differences
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


def calculate_manual_lj(atoms, forcefield, box, cutoff):
    """Manually calculate LJ energy to understand differences"""
    total_energy = 0.0
    n_atoms = len(atoms)
    
    for i in range(n_atoms):
        for j in range(i+1, n_atoms):
            # Calculate distance with PBC
            dx = atoms[i].x - atoms[j].x
            dy = atoms[i].y - atoms[j].y
            dz = atoms[i].z - atoms[j].z
            
            # Apply minimum image convention
            dx -= box[0] * round(dx / box[0])
            dy -= box[1] * round(dy / box[1])
            dz -= box[2] * round(dz / box[2])
            
            r = math.sqrt(dx*dx + dy*dy + dz*dz)
            
            if r < cutoff:
                # Get LJ parameters
                type_i = atoms[i].type
                type_j = atoms[j].type
                
                # For matrix indexing: index = type_i * numTypes + type_j
                num_types = forcefield.numTotalTypes
                idx = type_i * num_types + type_j
                
                epsilon = forcefield.ljEps[idx]
                sigma = forcefield.ljSigma[idx]
                
                # Calculate LJ energy
                sr6 = (sigma/r)**6
                energy = 4.0 * epsilon * (sr6*sr6 - sr6)
                total_energy += energy
                
                print(f"  Pair {i}-{j}: r={r:.3f}, eps={epsilon:.3f}, sigma={sigma:.3f}, E={energy:.6f}")
    
    return total_energy


@pytest.mark.skipif(not OPENMM_AVAILABLE, reason="OpenMM not available")
def test_simple_two_atom_system():
    """Test the simplest possible system - two atoms"""
    
    print("\nTwo-atom system test:")
    
    state = MCState()
    state.info.box = [10.0, 10.0, 10.0]  # Large box
    state.info.cutoff = 5.0
    
    ff = MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.ljEps = [1.0]
    ff.ljSigma = [0.35]
    state.forcefield = ff
    
    # Two atoms at exactly 0.5 nm
    atoms = []
    residues = []
    
    atom1 = MCAtom()
    atom1.x, atom1.y, atom1.z = 5.0, 5.0, 5.0
    atom1.charge = 0.0
    atom1.type = 0
    atoms.append(atom1)
    
    atom2 = MCAtom()
    atom2.x, atom2.y, atom2.z = 5.5, 5.0, 5.0
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
    
    # Calculate with PyGCMC
    elec, vdw, total = computeSystemEnergyCutoffFixed(state)
    
    # Calculate manually
    manual_lj = calculate_manual_lj(atoms, ff, state.info.box, state.info.cutoff)
    
    # Calculate with OpenMM
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
    openmm_lj = energy_state.getPotentialEnergy().value_in_unit(kilojoule_per_mole)
    
    print(f"\nResults:")
    print(f"  PyGCMC Fixed: {vdw:.6f} kJ/mol")
    print(f"  Manual calc:  {manual_lj:.6f} kJ/mol")
    print(f"  OpenMM:       {openmm_lj:.6f} kJ/mol")
    
    # Check residue energies
    print(f"\nResidue energies:")
    for i, res in enumerate(state.residues):
        print(f"  Residue {i}: {res.energy_vdw:.6f} kJ/mol")


@pytest.mark.skipif(not OPENMM_AVAILABLE, reason="OpenMM not available")
def test_mixed_types():
    """Test system with different atom types"""
    
    print("\n\nMixed atom types test:")
    
    state = MCState()
    state.info.box = [10.0, 10.0, 10.0]
    state.info.cutoff = 5.0
    
    ff = MCForceField()
    ff.numTotalTypes = 2
    ff.numMovementTypes = 2
    
    # Set up parameters
    eps0, eps1 = 1.0, 2.0
    sigma0, sigma1 = 0.3, 0.4
    
    # Full matrix
    ff.ljEps = [
        eps0,                          # 0-0
        math.sqrt(eps0 * eps1),        # 0-1
        math.sqrt(eps0 * eps1),        # 1-0
        eps1                           # 1-1
    ]
    ff.ljSigma = [
        sigma0,                        # 0-0
        (sigma0 + sigma1) / 2,         # 0-1
        (sigma0 + sigma1) / 2,         # 1-0
        sigma1                         # 1-1
    ]
    state.forcefield = ff
    
    # Two atoms of different types at 0.5 nm
    atoms = []
    residues = []
    
    atom1 = MCAtom()
    atom1.x, atom1.y, atom1.z = 5.0, 5.0, 5.0
    atom1.charge = 0.0
    atom1.type = 0
    atoms.append(atom1)
    
    atom2 = MCAtom()
    atom2.x, atom2.y, atom2.z = 5.5, 5.0, 5.0
    atom2.charge = 0.0
    atom2.type = 1
    atoms.append(atom2)
    
    for i in range(2):
        res = MCResidue()
        res.active = True
        res.fixed = True
        res.atomStart = i
        res.atomCount = 1
        res.type = atoms[i].type
        residues.append(res)
    
    state.atoms = atoms
    state.activeAtomCount = 2
    state.residues = residues
    state.activeResidueCount = 2
    
    # Calculate with PyGCMC
    elec, vdw, total = computeSystemEnergyCutoffFixed(state)
    
    # Calculate manually
    manual_lj = calculate_manual_lj(atoms, ff, state.info.box, state.info.cutoff)
    
    print(f"\nResults:")
    print(f"  PyGCMC Fixed: {vdw:.6f} kJ/mol")
    print(f"  Manual calc:  {manual_lj:.6f} kJ/mol")
    
    # Show parameters used
    print(f"\nParameters:")
    print(f"  Type 0: eps={eps0}, sigma={sigma0}")
    print(f"  Type 1: eps={eps1}, sigma={sigma1}")
    print(f"  Mixed: eps={ff.ljEps[1]:.3f}, sigma={ff.ljSigma[1]:.3f}")


if __name__ == "__main__":
    test_simple_two_atom_system()
    test_mixed_types()
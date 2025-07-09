"""
Analyze PME error components for medium complexity system
Focus on separating real space and reciprocal space contributions
"""

import numpy as np
import pygcmc
from pygcmc import MCState, MCAtom, MCResidue, MCForceField
from pygcmc import initializePMEParameters, computeSystemEnergyPME

try:
    from openmm import *
    from openmm.app import *
    from openmm.unit import *
    OPENMM_AVAILABLE = True
except ImportError:
    OPENMM_AVAILABLE = False


def create_medium_system():
    """Create the medium complexity system with 20 atoms"""
    # System parameters
    box_size = 3.5  # nm
    cutoff = 1.2    # nm
    
    # Create state
    state = MCState()
    state.info.box = [box_size, box_size, box_size]
    state.info.cutoff = cutoff
    
    # Force field - pure electrostatics (no LJ)
    ff = MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.ljEps = [0.0]
    ff.ljSigma = [0.3]
    state.forcefield = ff
    
    atoms = []
    residues = []
    atom_idx = 0
    
    # Add 12 ions
    ion_positions = [
        [1.0, 1.0, 1.0], [2.5, 1.0, 1.0], [1.0, 2.5, 1.0], [2.5, 2.5, 1.0],
        [1.0, 1.0, 2.5], [2.5, 1.0, 2.5], [1.0, 2.5, 2.5], [2.5, 2.5, 2.5],
        [1.75, 1.75, 1.0], [1.75, 1.75, 2.5], [1.0, 1.75, 1.75], [2.5, 1.75, 1.75]
    ]
    ion_charges = [1.0, -1.0, -1.0, 1.0, -1.0, 1.0, 1.0, -1.0, 1.0, -1.0, -1.0, 1.0]
    
    for i, (pos, charge) in enumerate(zip(ion_positions, ion_charges)):
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
    
    # Add ligand (8 atoms)
    ligand_start = atom_idx
    ligand_center = box_size / 2.0
    ligand_atoms = [
        ([0.0, 0.0, 0.0], -0.3),
        ([0.15, 0.0, 0.0], 0.1),
        ([-0.15, 0.0, 0.0], 0.1),
        ([0.0, 0.15, 0.0], 0.1),
        ([0.0, -0.15, 0.0], -0.4),
        ([0.0, 0.0, 0.15], 0.2),
        ([0.0, 0.0, -0.15], 0.1),
        ([0.2, 0.2, 0.0], 0.1),
    ]
    
    for rel_pos, charge in ligand_atoms:
        atom = MCAtom()
        atom.x = ligand_center + rel_pos[0]
        atom.y = ligand_center + rel_pos[1]
        atom.z = ligand_center + rel_pos[2]
        atom.charge = charge
        atom.type = 0
        atoms.append(atom)
    
    res = MCResidue()
    res.active = True
    res.fixed = False
    res.atomStart = ligand_start
    res.atomCount = len(ligand_atoms)
    res.type = 0
    residues.append(res)
    
    state.atoms = atoms
    state.activeAtomCount = len(atoms)
    state.residues = residues
    state.activeResidueCount = len(residues)
    
    return state


def analyze_pme_components():
    """Analyze PME energy components and compare with OpenMM"""
    
    print("PME COMPONENT ANALYSIS FOR MEDIUM COMPLEXITY SYSTEM")
    print("="*70)
    
    state = create_medium_system()
    print(f"\nSystem: {state.activeAtomCount} atoms in {state.info.box[0]} nm box")
    
    # Test parameters
    alpha_values = [2.5, 3.0, 3.5, 4.0]
    mesh_size = [32, 32, 32]
    spline_order = 4
    
    print("\nPyGCMC PME Components:")
    print(f"{'Alpha':>6} {'Real':>12} {'Recip':>12} {'Self':>12} {'Total':>12} {'Real%':>8} {'Recip%':>8}")
    print("-"*86)
    
    for alpha in alpha_values:
        # PyGCMC calculation
        pygcmc.setPMEParameters(alpha, mesh_size, spline_order)
        initializePMEParameters(state.info.cutoff, state.info.box, alpha, mesh_size, spline_order)
        computeSystemEnergyPME(state)
        
        real = state.ewald_energy.get('real_space', 0.0)
        recip = state.ewald_energy.get('reciprocal', 0.0)
        self = state.ewald_energy.get('self', 0.0)
        total = real + recip + self
        
        # Calculate percentages
        if total != 0:
            real_pct = abs(real/total) * 100
            recip_pct = abs(recip/total) * 100
        else:
            real_pct = recip_pct = 0
        
        print(f"{alpha:6.1f} {real:12.4f} {recip:12.4f} {self:12.4f} {total:12.4f} {real_pct:7.1f}% {recip_pct:7.1f}%")
    
    # OpenMM comparison for best alpha
    if OPENMM_AVAILABLE:
        print("\nOpenMM Comparison (alpha=3.5):")
        alpha = 3.5
        
        # PyGCMC
        pygcmc.setPMEParameters(alpha, mesh_size, spline_order)
        initializePMEParameters(state.info.cutoff, state.info.box, alpha, mesh_size, spline_order)
        computeSystemEnergyPME(state)
        
        pygcmc_total = (state.ewald_energy.get('real_space', 0.0) + 
                       state.ewald_energy.get('reciprocal', 0.0) + 
                       state.ewald_energy.get('self', 0.0))
        
        # OpenMM
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
        
        for atom in state.atoms:
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
        openmm_energy = energy_state.getPotentialEnergy().value_in_unit(kilojoule_per_mole)
        
        print(f"  PyGCMC: {pygcmc_total:.4f} kJ/mol")
        print(f"  OpenMM: {openmm_energy:.4f} kJ/mol")
        print(f"  Difference: {abs(pygcmc_total - openmm_energy):.4f} kJ/mol ({abs(pygcmc_total - openmm_energy)/abs(openmm_energy)*100:.2f}%)")



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


def analyze_distance_distribution():
    """Analyze charge-charge distances in the system"""
    
    print("\nDISTANCE DISTRIBUTION ANALYSIS")
    print("="*70)
    
    state = create_medium_system()
    
    # Calculate all pairwise distances
    distances = []
    energies = []
    
    for i in range(state.activeAtomCount):
        for j in range(i+1, state.activeAtomCount):
            dx = state.atoms[i].x - state.atoms[j].x
            dy = state.atoms[i].y - state.atoms[j].y
            dz = state.atoms[i].z - state.atoms[j].z
            
            # Apply minimum image convention
            box = state.info.box
            dx -= box[0] * round(dx / box[0])
            dy -= box[1] * round(dy / box[1])
            dz -= box[2] * round(dz / box[2])
            
            dist = np.sqrt(dx*dx + dy*dy + dz*dz)
            distances.append(dist)
            
            # Calculate Coulomb energy for this pair
            if dist > 0:
                energy = 138.935 * state.atoms[i].charge * state.atoms[j].charge / dist
                energies.append(energy)
    
    # Analyze by distance ranges
    ranges = [(0, 0.5), (0.5, 1.0), (1.0, 1.2), (1.2, 2.0), (2.0, 10.0)]
    
    print("\nEnergy contribution by distance range:")
    print(f"{'Range (nm)':>15} {'Count':>8} {'Energy (kJ/mol)':>15} {'% of Total':>12}")
    print("-"*50)
    
    total_energy = sum(energies)
    
    for r_min, r_max in ranges:
        count = 0
        energy_sum = 0
        for d, e in zip(distances, energies):
            if r_min <= d < r_max:
                count += 1
                energy_sum += e
        
        pct = (energy_sum / total_energy * 100) if total_energy != 0 else 0
        print(f"{r_min:6.1f} - {r_max:6.1f} {count:8d} {energy_sum:15.2f} {pct:11.1f}%")
    
    print(f"\nTotal pairs: {len(distances)}")
    print(f"Total Coulomb energy (no cutoff): {total_energy:.2f} kJ/mol")
    print(f"Cutoff: {state.info.cutoff} nm")
    
    # Count pairs within cutoff
    within_cutoff = 0
    for d in distances:
        if d <= state.info.cutoff:
            within_cutoff += 1
    print(f"Pairs within cutoff: {within_cutoff} ({within_cutoff/len(distances)*100:.1f}%)")


def test_mesh_convergence():
    """Test how PME converges with mesh size"""
    
    print("\nMESH SIZE CONVERGENCE TEST")
    print("="*70)
    
    state = create_medium_system()
    
    alpha = 3.5
    spline_order = 4
    mesh_sizes = [16, 24, 32, 48, 64]
    
    print(f"\nAlpha = {alpha}, Spline order = {spline_order}")
    print(f"{'Mesh':>6} {'Real':>12} {'Recip':>12} {'Total':>12} {'ΔTotal':>12}")
    print("-"*54)
    
    prev_total = None
    
    for mesh in mesh_sizes:
        mesh_size = [mesh, mesh, mesh]
        
        pygcmc.setPMEParameters(alpha, mesh_size, spline_order)
        initializePMEParameters(state.info.cutoff, state.info.box, alpha, mesh_size, spline_order)
        computeSystemEnergyPME(state)
        
        real = state.ewald_energy.get('real_space', 0.0)
        recip = state.ewald_energy.get('reciprocal', 0.0)
        self = state.ewald_energy.get('self', 0.0)
        total = real + recip + self
        
        if prev_total is not None:
            delta = total - prev_total
        else:
            delta = 0
        
        print(f"{mesh:6d} {real:12.4f} {recip:12.4f} {total:12.4f} {delta:12.6f}")
        
        prev_total = total


if __name__ == "__main__":
    analyze_pme_components()
    analyze_distance_distribution()
    test_mesh_convergence()
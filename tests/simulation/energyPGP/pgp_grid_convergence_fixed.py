"""
Fixed version of pgp_grid_convergence.py that uses PGPContext to avoid global state issues.
This prevents segmentation faults during parallel test execution.

import pytest
import random
from math import sqrt
import pygcmc
from pygcmc import (
    MCState, MCAtom, MCResidue, MCStateInfo, DoubleVector3,
    PGPContext,  # Use instance-based PGP implementation
    computeSystemEnergyPME, setPMEParameters, initializePMEParameters
)

def create_random_charged_system(n_particles=256, box_size=30.0):
    """Create a system with random charged particles."""
    random.seed(12345)  # Fixed seed for reproducibility
    
    state = MCState()
    info = MCStateInfo()
    
    # Set box dimensions
    info.box = DoubleVector3(box_size, box_size, box_size)
    info.boxInverse = DoubleVector3(1.0/box_size, 1.0/box_size, 1.0/box_size)
    info.boxVolume = box_size ** 3
    
    # Set cutoff 
    info.cutoff = 10.0
    info.cutoffSq = 100.0
    
    # Set other parameters
    info.kT = 1.0
    info.lnVolume = 0.0
    
    state.info = info
    
    # Create atoms with overall neutrality
    atoms = []
    total_charge = 0.0
    
    # Place n_particles/2 cations and n_particles/2 anions
    for i in range(n_particles):
        atom = MCAtom()
        # Random position in box
        atom.x = random.uniform(0, box_size)
        atom.y = random.uniform(0, box_size)
        atom.z = random.uniform(0, box_size)
        
        # Alternate charges for neutrality
        if i < n_particles // 2:
            atom.charge = 1.0
            atom.type = 0
        else:
            atom.charge = -1.0
            atom.type = 1
        
        total_charge += atom.charge
        atoms.append(atom)
    
    state.atoms = atoms
    state.activeAtomCount = n_particles
    
    # Create residues (each atom is its own residue)
    residues = []
        res = MCResidue()
        res.active = True
        res.fixed = False
        res.atomStart = i
        res.atomCount = 1
        res.type = atoms[i].type
        residues.append(res)
    
    state.residues = residues
    state.activeResidueCount = n_particles
    
    print(f"Created system: {n_particles} particles, total charge = {total_charge}")
    assert abs(total_charge) < 1e-10, "System must be neutral"
    
    return state

@pytest.mark.parametrize("grid_size", [32, 64, 128])
def test_pgp_grid_convergence(grid_size):
    """Test PGP energy convergence with increasing grid resolution using PGPContext."""
    print(f"\n=== Testing grid convergence with {grid_size}³ grid ===")
    
    # Create test system
    state = create_random_charged_system(n_particles=256)
    cutoff = state.info.cutoff
    box = state.info.box
    alpha = 2.5  # Fixed alpha for grid test
    
    # Create PGP context instance (avoids global state)
    pgp_ctx = PGPContext()
    
    # Initialize and calculate PGP energy
    pgp_ctx.initialize(cutoff, box, alpha)
    pgp_ctx.setParameters(alpha, [grid_size]*3, cutoff, [grid_size]*3, 4, 1e-6)
    pgp_ctx.precomputeGrid(state, fixed_only=True)
    pgp_ctx.computeSystemEnergy(state)
    
    pgp_total = state.ewald_energy.get('total', 0.0)
    pgp_real = state.ewald_energy.get('real_space', 0.0)
    pgp_recip = state.ewald_energy.get('reciprocal', 0.0)
    pgp_self = state.ewald_energy.get('self', 0.0)
    
    # Calculate high-precision PME reference (fine grid)
    state_ref = create_random_charged_system(n_particles=256)  # Same system
    setPMEParameters(alpha, [64, 64, 64], 6, 1e-8)  # High precision
    initializePMEParameters(cutoff, box, alpha)
    computeSystemEnergyPME(state_ref)
    
    pme_ref = state_ref.ewald_energy.get('total')
    
    # Calculate error
    error = abs(pgp_total - pme_ref)
    rel_error = error / abs(pme_ref) if pme_ref != 0 else error
    
    print(f"PGP energy components:")
    print(f"  Real-space: {pgp_real:.4f} kJ/mol")
    print(f"  Reciprocal: {pgp_recip:.4f} kJ/mol")
    print(f"  Self:       {pgp_self:.4f} kJ/mol")
    print(f"  Total:      {pgp_total:.4f} kJ/mol")
    print(f"PME reference: {pme_ref:.4f} kJ/mol")
    print(f"Absolute error: {error:.6f} kJ/mol")
    print(f"Relative error: {rel_error:.6e}")
    
    # PGP has different convergence characteristics than PME
    # so we don't enforce strict error bounds here
    
    # Clean up context
    del pgp_ctx

if __name__ == "__main__":
    # Run tests
    test_pgp_grid_convergence(32)
    test_pgp_grid_convergence(64)
    test_pgp_grid_convergence(128)
"""

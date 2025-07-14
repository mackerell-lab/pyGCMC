# tests/simulation/energyPGP/vspme_energy_diagnostics.py
"""Energy calculation diagnostic helper functions for combined tests."""

import math
import pygcmc
from . import pgp_wrapper
from .pgp_wrapper import computeSystemVdwEnergyCutoff, computeSystemEnergyPME, computeSystemEnergyPGP, computeMovementEnergyPME
from .pgp_wrapper import computeMovementEnergyPGP, setPMEParameters, setPGPParameters, initializePMEParameters
from .pgp_wrapper import precomputeGridPotential
from pygcmc import MCState, MCAtom, MCResidue
from pygcmc import MCForceField, MCMovementResidueInfo

def create_diagnostic_test_system():
    """Create a simple two-atom test system for energy calculation diagnostics."""
    # Create a simple test system
    box_size = 4.0  # nm
    cutoff = 1.2    # nm
    
    # Create MCState
    state = MCState()
    state.info.box = [box_size, box_size, box_size]
    state.info.setTemperature(300.0)
    state.info.cutoff = cutoff
    
    # Set force field parameters
    ff = MCForceField()
    ff.numTotalTypes = 2
    sigma = 0.4  # nm
    eps = 0.02   # kJ/mol
    ff.ljSigma = [sigma, sigma, sigma, sigma]
    ff.ljEps = [eps, eps, eps, eps]
    state.forcefield = ff
    
    # Create atoms
    atoms = []
    
    # Fixed atom
    fixed_atom = MCAtom()
    fixed_atom.x = 2.0
    fixed_atom.y = 2.0
    fixed_atom.z = 2.0
    fixed_atom.charge = 1.0
    fixed_atom.type = 0
    atoms.append(fixed_atom)
    
    # Moving atom - initial position
    moving_atom = MCAtom()
    moving_atom.x = 3.0  # Initial distance 1.0nm
    moving_atom.y = 2.0
    moving_atom.z = 2.0
    moving_atom.charge = -1.0
    moving_atom.type = 1
    atoms.append(moving_atom)
    
    # Create residues
    residues = []
    
    # Fixed residue
    fixed_res = MCResidue()
    fixed_res.atomStart = 0
    fixed_res.atomCount = 1
    fixed_res.active = True
    fixed_res.fixed = True
    fixed_res.energy_vdw = 0.0
    fixed_res.energy_elec = 0.0
    residues.append(fixed_res)
    
    # Moving residue
    moving_res = MCResidue()
    moving_res.atomStart = 1
    moving_res.atomCount = 1
    moving_res.active = True
    moving_res.fixed = False
    moving_res.energy_vdw = 0.0
    moving_res.energy_elec = 0.0
    residues.append(moving_res)
    
    # Set system
    state.atoms = atoms
    state.residues = residues
    state.activeAtomCount = len(atoms)
    state.activeResidueCount = len(residues)
    
    return state, box_size, cutoff, sigma, eps

def setup_energy_parameters(state, box_size, cutoff):
    """Initialize PME and PGP parameters for energy calculations."""
    mesh_size = [16, 16, 16]
    box = [box_size, box_size, box_size]
    alpha = 0.25
    
    setPMEParameters(alpha=alpha, meshSize=mesh_size, splineOrder=4, tolerance=1e-5)
    initializePMEParameters(cutoff, box, alpha)
    setPGPParameters(alpha, mesh_size, state.info.cutoff, mesh_size, 4, 1e-6)
    precomputeGridPotential(state)
    
    # Set moving residue information
    state.movementResidues.clear()
    movement_info = MCMovementResidueInfo()
    movement_info.startIndex = 1
    movement_info.activeCount = 1
    state.movementResidues.append(movement_info)
    
    return alpha, mesh_size

def analyze_initial_state_energies(state, sigma, eps):
    """Analyze initial state energy components using different calculation methods."""
    print("\n--- Phase 1: Initial State Detailed Check ---")
    print("\n1.1 Calculate LJ energy using computeSystemVdwEnergyCutoff")
    
    # Reset energies
    for res in state.residues:
        res.energy_vdw = 0.0
        res.energy_elec = 0.0
    
    # Direct LJ energy calculation
    computeSystemVdwEnergyCutoff(state)
    
    # Print LJ energy for each residue
    vdw_total = 0.0
    for i, res in enumerate(state.residues):
        vdw_total += res.energy_vdw
        print(f"  Residue {i} LJ energy: {res.energy_vdw:.6f} kJ/mol")
    print(f"  Total LJ energy: {vdw_total:.6f} kJ/mol")
    
    print("\n1.2 Calculate system total energy using computeSystemEnergyPME")
    
    # Reset energies
    for res in state.residues:
        res.energy_vdw = 0.0
        res.energy_elec = 0.0
    
    # Calculate total energy using PME
    computeSystemEnergyPME(state)
    
    # View results
    pme_total_energy = state.ewald_energy.get("total", 0.0)
    pme_reciprocal = state.ewald_energy.get("reciprocal", 0.0)
    pme_real_space = state.ewald_energy.get("real_space", 0.0)
    pme_self = state.ewald_energy.get("self", 0.0)
    
    # Check residue LJ energies again
    pme_vdw_total = 0.0
    for i, res in enumerate(state.residues):
        pme_vdw_total += res.energy_vdw
        print(f"  Residue {i} LJ energy after PME: {res.energy_vdw:.6f} kJ/mol")
    
    # Print all energy components
    print(f"  PME total energy: {pme_total_energy:.6f} kJ/mol")
    print(f"  PME reciprocal energy: {pme_reciprocal:.6f} kJ/mol")
    print(f"  PME real space energy: {pme_real_space:.6f} kJ/mol")
    print(f"  PME self correction: {pme_self:.6f} kJ/mol")
    print(f"  PME calculated LJ energy: {pme_vdw_total:.6f} kJ/mol")
    print(f"  PME energy component sum: {pme_reciprocal + pme_real_space + pme_self + pme_vdw_total:.6f} kJ/mol")
    
    # Check if ewald_energy contains LJ energy
    print(f"  ewald_energy dictionary keys: {list(state.ewald_energy.keys())}")
    
    print("\n1.3 Calculate system total energy using computeSystemEnergyPGP")
    
    # Reset energies
    for res in state.residues:
        res.energy_vdw = 0.0
        res.energy_elec = 0.0
    
    # Calculate total energy using PGP
    computeSystemEnergyPGP(state)
    
    # View results
    pgp_total_energy = state.ewald_energy.get("total", 0.0)
    pgp_reciprocal = state.ewald_energy.get("reciprocal", 0.0)
    pgp_real_space = state.ewald_energy.get("real_space", 0.0)
    pgp_self = state.ewald_energy.get("self", 0.0)
    
    # Check residue LJ energies again
    pgp_vdw_total = 0.0
    for i, res in enumerate(state.residues):
        pgp_vdw_total += res.energy_vdw
        print(f"  Residue {i} LJ energy after PGP: {res.energy_vdw:.6f} kJ/mol")
    
    # Print all energy components
    print(f"  PGP total energy: {pgp_total_energy:.6f} kJ/mol")
    print(f"  PGP reciprocal energy: {pgp_reciprocal:.6f} kJ/mol")
    print(f"  PGP real space energy: {pgp_real_space:.6f} kJ/mol")
    print(f"  PGP self correction: {pgp_self:.6f} kJ/mol")
    print(f"  PGP calculated LJ energy: {pgp_vdw_total:.6f} kJ/mol")
    print(f"  PGP energy component sum: {pgp_reciprocal + pgp_real_space + pgp_self + pgp_vdw_total:.6f} kJ/mol")
    
    return vdw_total


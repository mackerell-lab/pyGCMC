"""
Detailed analysis of PME error sources in medium complexity systems

This script analyzes why medium complexity systems (20 atoms) have 4-6% error
by separating real space (short-range) and reciprocal space (long-range) contributions.
"""

import pytest
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pygcmc
from pygcmc import MCState, MCAtom, MCResidue, MCForceField
from pygcmc import initializePMEParameters, computeSystemEnergyPME, computeSystemEnergyPBC

try:
    from openmm import *
    from openmm.app import *
    from openmm.unit import *
    OPENMM_AVAILABLE = True
except ImportError:
    OPENMM_AVAILABLE = False

from pme_medium_complexity import create_medium_complexity_system, calculate_openmm_energy_medium


def analyze_energy_components(state):
    """Analyze different energy components for the medium complexity system"""
    
    print("\n" + "="*80)
    print("DETAILED ENERGY COMPONENT ANALYSIS")
    print("="*80)
    
    # System info
    print(f"\nSystem: {state.activeAtomCount} atoms")
    print(f"Box: [{state.info.box[0]}, {state.info.box[1]}, {state.info.box[2]}] nm")
    print(f"Cutoff: {state.info.cutoff} nm")
    
    # Count atom types and distances
    analyze_system_structure(state)
    
    # 1. Direct calculation (all interactions within cutoff)
    print("\n1. CUTOFF CALCULATION (pairwise interactions within cutoff):")
    # Calculate energy directly using pairwise summation
    cutoff_energy = 0.0
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
            
            dist = math.sqrt(dx*dx + dy*dy + dz*dz)
            if dist <= state.info.cutoff and dist > 0:
                # Coulomb energy
                cutoff_energy += 138.935 * state.atoms[i].charge * state.atoms[j].charge / dist
    
    print(f"   Total cutoff energy: {cutoff_energy:.4f} kJ/mol")
    
    # 2. PME calculation with different parameters
    print("\n2. PME CALCULATIONS:")
    
    test_configs = [
        (2.5, [32, 32, 32], 4),
        (3.0, [32, 32, 32], 4),
        (3.5, [32, 32, 32], 4),  # Best from previous test
        (4.0, [32, 32, 32], 4),
    ]
    
    for alpha, mesh_size, spline_order in test_configs:
        print(f"\n   Alpha={alpha}, Mesh={mesh_size[0]}x{mesh_size[1]}x{mesh_size[2]}:")
        
        # PyGCMC PME
        pygcmc.setPMEParameters(alpha, mesh_size, spline_order)
        initializePMEParameters(state.info.cutoff, state.info.box, alpha, mesh_size, spline_order)
        computeSystemEnergyPME(state)
        
        real_space = state.ewald_energy.get('real_space', 0.0)
        reciprocal = state.ewald_energy.get('reciprocal', 0.0)
        self_energy = state.ewald_energy.get('self', 0.0)
        pme_total = real_space + reciprocal + self_energy
        
        print(f"      Real space:  {real_space:12.4f} kJ/mol")
        print(f"      Reciprocal:  {reciprocal:12.4f} kJ/mol")
        print(f"      Self:        {self_energy:12.4f} kJ/mol")
        print(f"      PME Total:   {pme_total:12.4f} kJ/mol")
        
        # Compare with OpenMM if available
        if OPENMM_AVAILABLE:
            openmm_energy = calculate_openmm_energy_medium(state, alpha)
            if openmm_energy is not None:
                diff = abs(pme_total - openmm_energy)
                rel_diff = diff / abs(openmm_energy) * 100
                print(f"      OpenMM:      {openmm_energy:12.4f} kJ/mol")
                print(f"      Difference:  {diff:12.4f} kJ/mol ({rel_diff:.2f}%)")
        
        # Analyze error distribution
        print(f"      Real/Total ratio: {abs(real_space/pme_total)*100:.1f}%")
        print(f"      Recip/Total ratio: {abs(reciprocal/pme_total)*100:.1f}%")
        print(f"      Self/Total ratio: {abs(self_energy/pme_total)*100:.1f}%")


def analyze_system_structure(state):
    """Analyze the structure of the system"""
    
    # Count charges
    positive_charges = 0
    negative_charges = 0
    total_charge = 0.0
    
    for i in range(state.activeAtomCount):
        charge = state.atoms[i].charge
        total_charge += charge
        if charge > 0:
            positive_charges += 1
        elif charge < 0:
            negative_charges += 1
    
    print(f"\nCharge distribution:")
    print(f"  Positive charges: {positive_charges}")
    print(f"  Negative charges: {negative_charges}")
    print(f"  Total charge: {total_charge:.6f}")
    
    # Calculate distances between all pairs
    distances = []
    close_pairs = 0  # < 0.5 nm
    medium_pairs = 0  # 0.5-1.0 nm
    far_pairs = 0     # > 1.0 nm
    
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
            
            dist = math.sqrt(dx*dx + dy*dy + dz*dz)
            distances.append(dist)
            
            if dist < 0.5:
                close_pairs += 1
            elif dist < 1.0:
                medium_pairs += 1
            else:
                far_pairs += 1
    
    print(f"\nDistance distribution ({len(distances)} pairs):")
    print(f"  Close (<0.5 nm): {close_pairs} pairs")
    print(f"  Medium (0.5-1.0 nm): {medium_pairs} pairs")
    print(f"  Far (>1.0 nm): {far_pairs} pairs")
    print(f"  Min distance: {min(distances):.3f} nm")
    print(f"  Max distance: {max(distances):.3f} nm")
    print(f"  Cutoff: {state.info.cutoff} nm")



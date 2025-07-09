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


def test_distance_dependent_errors():
    """Test how error varies with distance cutoffs"""
    
    state = create_medium_complexity_system()
    
    print("\n" + "="*80)
    print("DISTANCE-DEPENDENT ERROR ANALYSIS")
    print("="*80)
    
    # Test with different cutoffs
    cutoffs = [0.8, 1.0, 1.2, 1.4, 1.6]
    
    for cutoff in cutoffs:
        print(f"\nCutoff = {cutoff} nm:")
        state.info.cutoff = cutoff
        
        # Cutoff calculation
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
                if dist <= cutoff and dist > 0:
                    cutoff_energy += 138.935 * state.atoms[i].charge * state.atoms[j].charge / dist
        
        print(f"  Cutoff energy: {cutoff_energy:.4f} kJ/mol")
        
        # PME with optimal parameters
        alpha = 3.5
        mesh_size = [32, 32, 32]
        spline_order = 4
        
        pygcmc.setPMEParameters(alpha, mesh_size, spline_order)
        initializePMEParameters(cutoff, state.info.box, alpha, mesh_size, spline_order)
        computeSystemEnergyPME(state)
        
        pme_total = (state.ewald_energy.get('real_space', 0.0) + 
                     state.ewald_energy.get('reciprocal', 0.0) + 
                     state.ewald_energy.get('self', 0.0))
        
        print(f"  PME energy: {pme_total:.4f} kJ/mol")
        print(f"  Difference: {abs(cutoff_energy - pme_total):.4f} kJ/mol")


def test_real_space_accuracy():
    """Test real space calculation accuracy by comparing with direct calculation at short range"""
    
    state = create_medium_complexity_system()
    
    print("\n" + "="*80)
    print("REAL SPACE ACCURACY TEST")
    print("="*80)
    
    # Use small cutoff to isolate real space contribution
    state.info.cutoff = 0.6  # Very small cutoff
    
    # Cutoff calculation with small cutoff
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
                cutoff_energy += 138.935 * state.atoms[i].charge * state.atoms[j].charge / dist
    
    print(f"\nCutoff energy (cutoff={state.info.cutoff} nm): {cutoff_energy:.4f} kJ/mol")
    
    # PME real space should match direct for small cutoff
    alpha = 5.0  # Large alpha for sharp cutoff
    mesh_size = [32, 32, 32]
    spline_order = 4
    
    pygcmc.setPMEParameters(alpha, mesh_size, spline_order)
    initializePMEParameters(state.info.cutoff, state.info.box, alpha, mesh_size, spline_order)
    computeSystemEnergyPME(state)
    
    real_space = state.ewald_energy.get('real_space', 0.0)
    print(f"PME real space: {real_space:.4f} kJ/mol")
    print(f"Difference: {abs(cutoff_energy - real_space):.4f} kJ/mol")
    
    if abs(cutoff_energy) > 1e-6:
        rel_diff = abs(cutoff_energy - real_space) / abs(cutoff_energy) * 100
        print(f"Relative difference: {rel_diff:.2f}%")


def test_reciprocal_space_convergence():
    """Test reciprocal space convergence with mesh size"""
    
    state = create_medium_complexity_system()
    
    print("\n" + "="*80)
    print("RECIPROCAL SPACE CONVERGENCE TEST")
    print("="*80)
    
    alpha = 3.5
    spline_order = 4
    mesh_sizes = [[16, 16, 16], [24, 24, 24], [32, 32, 32], [48, 48, 48], [64, 64, 64]]
    
    print(f"\nAlpha = {alpha}, testing different mesh sizes:")
    print(f"{'Mesh':>12} {'Real':>12} {'Reciprocal':>12} {'Self':>12} {'Total':>12}")
    print("-" * 60)
    
    reciprocal_values = []
    
    for mesh_size in mesh_sizes:
        pygcmc.setPMEParameters(alpha, mesh_size, spline_order)
        initializePMEParameters(state.info.cutoff, state.info.box, alpha, mesh_size, spline_order)
        computeSystemEnergyPME(state)
        
        real = state.ewald_energy.get('real_space', 0.0)
        recip = state.ewald_energy.get('reciprocal', 0.0)
        self = state.ewald_energy.get('self', 0.0)
        total = real + recip + self
        
        mesh_str = f"{mesh_size[0]}x{mesh_size[1]}x{mesh_size[2]}"
        print(f"{mesh_str:>12} {real:12.4f} {recip:12.4f} {self:12.4f} {total:12.4f}")
        
        reciprocal_values.append(recip)
    
    # Check convergence
    print("\nReciprocal space convergence:")
    for i in range(1, len(reciprocal_values)):
        diff = abs(reciprocal_values[i] - reciprocal_values[i-1])
        print(f"  {mesh_sizes[i-1][0]}→{mesh_sizes[i][0]}: Δ = {diff:.6f} kJ/mol")


if __name__ == "__main__":
    # Create system once
    state = create_medium_complexity_system()
    
    # Run all analyses
    analyze_energy_components(state)
    test_distance_dependent_errors()
    test_real_space_accuracy()
    test_reciprocal_space_convergence()
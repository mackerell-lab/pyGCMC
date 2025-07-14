# tests/simulation/energyPGP/vspme_two_atom.py
"""Simple two-atom system energy comparison test."""

import pytest
import math
import pygcmc
from . import pgp_wrapper
from .pgp_wrapper import setPGPParameters, precomputeGridPotential, computeSystemEnergyPGP
from .pgp_wrapper import computeMovementEnergyPGP, computeSystemEnergyPME, computeMovementEnergyPME
from .pgp_wrapper import computeSystemVdwEnergyCutoff
import sys
from .vspme_two_atom_helpers import create_simple_two_atom_state
from pygcmc import MCState, MCInfo, MCAtom
from pygcmc import MCResidue, MCForceField, MCMovementResidueInfo

def test_simple_two_atom_system():
    """Test comparison of theoretical energy vs calculated energy in simple two-atom system"""
    print("\n--- Test: Simple Two-Atom System Energy Calculation Comparison (Independent States) ---")
    
    # Set LJ parameters
    sigma = 0.400  # nm
    epsilon = 0.020  # kJ/mol
    print(f"LJ parameters: sigma = {sigma:.3f} nm, epsilon = {epsilon:.3f} kJ/mol\n")
    
    # Test multiple distance points - use more reasonable distance points, avoid extremely short distances
    distances = [1.0, 0.8, 0.7, 0.6, 0.5, 0.45]  # nm
    
    # Table headers
    print("Compare energy calculations at different distances:")
    headers = ["Dist(nm)", "Theor LJ", "Theor Coul", "Theor Total", "Direct LJ", "PME Total", "PME LJ", "PGP Total", "PGP LJ", "PME Delta", "PGP Delta"] 
    fmt = "{:10.4f} | {:12.4f} | {:12.4f} | {:12.4f} | {:12.4f} | {:12.4f} | {:12.4f} | {:12.4f} | {:12.4f} | {:12.4f} | {:12.4f}"
    
    # Print table header
    print("    " + " | ".join(headers))
    print("    " + " | ".join(["-" * 10] * len(headers)))
    
    # Record energy at previous distance point for calculating energy changes
    prev_pme_total = None
    prev_pgp_total = None
    prev_distance = None
    first_point = True
    
    for distance in distances:
        try:
            # Create new MCState, one independent system state for each distance
            state = create_simple_two_atom_state(distance, 4.0, 1.2)
            
            # ------------ Theoretical energy calculation ------------
            # Calculate theoretical LJ energy
            r = distance
            r6 = (sigma / r) ** 6
            r12 = r6 * r6
            theoretical_lj = 4 * epsilon * (r12 - r6)
            
            # Calculate theoretical Coulomb energy
            q1 = 1.0  # Fixed ion charge
            q2 = -1.0  # Moving ion charge
            coulomb = 138.935456  # kJ·mol^-1·nm·e^-2
            theoretical_coulomb = coulomb * q1 * q2 / r
            
            # Total theoretical energy
            theoretical_total = theoretical_lj + theoretical_coulomb
            
            # ------------ Direct LJ energy calculation ------------
            # First reset all energies
            for res in state.residues:
                res.energy_vdw = 0.0
                res.energy_elec = 0.0
                
            # Calculate LJ energy
            direct_lj_result = computeSystemVdwEnergyCutoff(state)
            
            # Get sum of LJ energies for each residue
            direct_lj = 0.0
            for res in state.residues:
                direct_lj += res.energy_vdw
            
            # ------------ Calculate system energy using PME ------------
            try:
                # Calculate PME system energy
                pme_result = computeSystemEnergyPME(state)
                
                # Get total energy and energy components
                pme_total = state.ewald_energy.get("total", 0.0)
                
                # Get sum of LJ energies for each residue
                pme_lj = 0.0
                for res in state.residues:
                    pme_lj += res.energy_vdw
            except Exception as e:
                print(f"PME system energy calculation error: {e}")
                pme_total = float('nan')
                pme_lj = float('nan')
                
            # ------------ Calculate system energy using PGP ------------
            try:
                # First reset all energies
                for res in state.residues:
                    res.energy_vdw = 0.0
                    res.energy_elec = 0.0
                    
                # Reset PGP parameters and precompute potential field
                meshSize = [16, 16, 16]
                potentialGridSize = [16, 16, 16]
                setPGPParameters(0.2, meshSize, 1.2, 
                                potentialGridSize, 4, 1e-4)
                precomputeGridPotential(state, True)
                
                # Calculate PGP system energy
                pgp_result = computeSystemEnergyPGP(state)
                
                # Get total energy
                pgp_total = state.ewald_energy.get("total", 0.0)
                
                # Get sum of LJ energies for each residue
                pgp_lj = 0.0
                for res in state.residues:
                    pgp_lj += res.energy_vdw
            except Exception as e:
                print(f"PGP system energy calculation error: {e}")
                pgp_total = float('nan')
                pgp_lj = float('nan')
                
            # ------------ Calculate energy changes ------------
            # Set moving part
            state.movementResidues.clear()
            movement_info = MCMovementResidueInfo()
            movement_info.startIndex = 1  # Moving ion residue index
            movement_info.activeCount = 1
            state.movementResidues.append(movement_info)
            
            # Calculate PME movement energy change
            try:
                # Calculate PME movement energy
                pme_movement_result = computeMovementEnergyPME(state)
                pme_delta = state.ewald_energy.get("total", 0.0)
                
                # If this is the first point, movement energy should be close to 0 or close to absolute energy
                if first_point:
                    if abs(pme_delta) < 1e-6:
                        print("Note: Movement function returns energy change values, initial point is 0")
                    first_point = False
            except Exception as e:
                print(f"PME movement energy calculation error: {e}")
                pme_delta = float('nan')
                
            # Calculate PGP movement energy change
            try:
                # Calculate PGP movement energy
                pgp_movement_result = computeMovementEnergyPGP(state)
                pgp_delta = state.ewald_energy.get("total", 0.0)
            except Exception as e:
                print(f"PGP movement energy calculation error: {e}")
                pgp_delta = float('nan')
            
            # Calculate energy difference from previous point - better understand movement function behavior
            if prev_distance is not None:
                expected_pme_delta = pme_total - prev_pme_total if prev_pme_total is not None else None
                expected_pgp_delta = pgp_total - prev_pgp_total if prev_pgp_total is not None else None
                
                if expected_pme_delta is not None and expected_pgp_delta is not None:
                    delta_ratio_pme = abs(pme_delta / expected_pme_delta) if abs(expected_pme_delta) > 1e-6 else float('nan')
                    delta_ratio_pgp = abs(pgp_delta / expected_pgp_delta) if abs(expected_pgp_delta) > 1e-6 else float('nan')
                    print(f"    Distance change: {prev_distance:.4f} → {distance:.4f} nm")
                    print(f"    Expected PME energy change: {expected_pme_delta:.4f}, actual: {pme_delta:.4f}, ratio: {delta_ratio_pme:.4f}")
                    print(f"    Expected PGP energy change: {expected_pgp_delta:.4f}, actual: {pgp_delta:.4f}, ratio: {delta_ratio_pgp:.4f}")
            
            # Record current energy as reference for next point
            prev_pme_total = pme_total
            prev_pgp_total = pgp_total
            prev_distance = distance
            
            # Print result row
            print("    " + fmt.format(
                distance, theoretical_lj, theoretical_coulomb, theoretical_total,
                direct_lj, pme_total, pme_lj, pgp_total, pgp_lj, pme_delta, pgp_delta
            ))
            
        except Exception as e:
            print(f"Error occurred while processing distance {distance} nm: {e}")
            import traceback
            traceback.print_exc()
            continue
    
    print("\nEnergy Analysis:")
    print("1. Direct LJ calculation (direct_lj) should closely match theoretical LJ")
    print("2. PME/PGP system total energy should include Coulomb and LJ energies")
    print("3. LJ components in PME/PGP should be consistent with direct LJ calculation")
    print("4. PME/PGP Movement functions calculate energy changes, not absolute energy values")
    print("5. Energy differences returned by PME/PGP Movement should be close to total energy differences caused by distance changes")
    
    print("\n--- Test Complete ---")

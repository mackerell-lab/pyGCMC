#!/usr/bin/env python
"""Test Thole screening with water molecules"""

import sys
import numpy as np
sys.path.append('/home/zhaomt/gcmc/test107/pygcmc_dev/build/modules/bindings')

import pygcmc

def test_thole_water_screening():
    """Test Thole screening between water molecules"""
    print("=== Testing Thole Screening with Water ===\n")
    
    # Initialize
    pygcmc.initializeDrudeForce()
    
    # Create state
    state = pygcmc.MCState()
    state.info.cutoff = 2.0
    state.info.box = [10.0, 10.0, 10.0]
    
    # Force field - just O and D types
    state.forcefield.numTotalTypes = 2
    state.forcefield.numMovementTypes = 2
    
    # No LJ for this test
    state.forcefield.ljSigma = [0.0] * 4
    state.forcefield.ljEps = [0.0] * 4
    
    # SWM4-NDP parameters
    qO_core = 1.66260
    qD = -1.76260
    alpha_nm3 = 0.0013
    thole = 1.3  # From CHARMM parameters
    
    print("Water parameters (SWM4-NDP):")
    print(f"  O core charge: {qO_core}")
    print(f"  Drude charge: {qD}")
    print(f"  Net O charge: {qO_core + qD:.4f}")
    print(f"  Polarizability: {alpha_nm3*1000} Å³")
    print(f"  Thole parameter: {thole}")
    print()
    
    # Test at different distances
    distances = [0.25, 0.3, 0.35, 0.4, 0.5, 0.6]  # nm
    
    for dist in distances:
        # Clear and reinitialize
        pygcmc.clearDrudeForce()
        pygcmc.initializeDrudeForce()
        
        # Create two oxygen atoms with Drude particles
        atoms = []
        
        # Water 1 - just O and D
        o1 = pygcmc.MCAtom()
        o1.x = 0.0
        o1.y = 0.0
        o1.z = 0.0
        o1.charge = qO_core
        o1.type = 0
        atoms.append(o1)
        
        d1 = pygcmc.MCAtom()
        d1.x = 0.0
        d1.y = 0.0
        d1.z = 0.0
        d1.charge = qD
        d1.type = 1
        atoms.append(d1)
        
        # Water 2 - just O and D
        o2 = pygcmc.MCAtom()
        o2.x = dist
        o2.y = 0.0
        o2.z = 0.0
        o2.charge = qO_core
        o2.type = 0
        atoms.append(o2)
        
        d2 = pygcmc.MCAtom()
        d2.x = dist
        d2.y = 0.0
        d2.z = 0.0
        d2.charge = qD
        d2.type = 1
        atoms.append(d2)
        
        state.atoms = atoms
        state.activeAtomCount = 4
        
        # Two residues
        res1 = pygcmc.MCResidue()
        res1.atomStart = 0
        res1.atomCount = 2
        res1.active = True
        res1.type = 0
        
        res2 = pygcmc.MCResidue()
        res2.atomStart = 2
        res2.atomCount = 2
        res2.active = True
        res2.type = 0
        
        state.residues = [res1, res2]
        state.activeResidueCount = 2
        
        # Add Drude particles
        pygcmc.addDrudeParticle(
            drudeIndex=1,
            parentIndex=0,
            charge=qD,
            polarizability=alpha_nm3
        )
        
        pygcmc.addDrudeParticle(
            drudeIndex=3,
            parentIndex=2,
            charge=qD,
            polarizability=alpha_nm3
        )
        
        # Test WITH and WITHOUT screening
        params = pygcmc.DrudeSCFParams()
        params.tolerance = 1e-8
        params.maxIterations = 100
        params.dampingFactor = 0.5
        params.maxDrudeDistance = 0.02
        pygcmc.setDrudeSCFParameters(params)
        
        # First WITHOUT screening
        energy_no_screen, _ = pygcmc.computeSystemEnergyDrude(state)
        
        # Get Drude positions after SCF
        d1_pos_no_screen = [state.atoms[1].x, state.atoms[1].y, state.atoms[1].z]
        d2_pos_no_screen = [state.atoms[3].x, state.atoms[3].y, state.atoms[3].z]
        
        # Reset Drude positions
        state.atoms[1].x = 0.0
        state.atoms[3].x = dist
        
        # Now WITH screening
        pygcmc.addDrudeScreenedPair(
            dipole1=0,  # First Drude particle
            dipole2=1,  # Second Drude particle
            thole=thole
        )
        
        energy_screened, _ = pygcmc.computeSystemEnergyDrude(state)
        
        # Get Drude positions after SCF with screening
        d1_pos_screen = [state.atoms[1].x, state.atoms[1].y, state.atoms[1].z]
        d2_pos_screen = [state.atoms[3].x, state.atoms[3].y, state.atoms[3].z]
        
        # Calculate u parameter
        u = dist * thole / (alpha_nm3 ** (1.0/6.0))
        screening_factor = 1.0 - (1.0 + 0.5 * u) * np.exp(-u)
        
        print(f"Distance: {dist*10:.1f} Å")
        print(f"  u parameter: {u:.3f}")
        print(f"  Screening factor: {screening_factor:.3f}")
        print(f"  Energy (no screening): {energy_no_screen:.3f} kJ/mol")
        print(f"  Energy (with screening): {energy_screened:.3f} kJ/mol")
        print(f"  Energy difference: {energy_screened - energy_no_screen:.3f} kJ/mol")
        
        # Compare Drude displacements
        d1_disp_no = np.sqrt(sum(d**2 for d in d1_pos_no_screen))
        d1_disp_yes = np.sqrt(sum(d**2 for d in d1_pos_screen))
        d2_disp_no = np.sqrt((d2_pos_no_screen[0]-dist)**2 + d2_pos_no_screen[1]**2 + d2_pos_no_screen[2]**2)
        d2_disp_yes = np.sqrt((d2_pos_screen[0]-dist)**2 + d2_pos_screen[1]**2 + d2_pos_screen[2]**2)
        
        print(f"  Drude 1 displacement: {d1_disp_no*1000:.3f} pm → {d1_disp_yes*1000:.3f} pm")
        print(f"  Drude 2 displacement: {d2_disp_no*1000:.3f} pm → {d2_disp_yes*1000:.3f} pm")
        print()

def test_screening_formula():
    """Verify the Thole screening formula"""
    print("\n=== Verifying Thole Screening Formula ===\n")
    
    # Test the screening function shape
    u_values = np.linspace(0, 10, 100)
    screening_values = []
    
    for u in u_values:
        if u == 0:
            screening = 0.0
        else:
            screening = 1.0 - (1.0 + 0.5 * u) * np.exp(-u)
        screening_values.append(screening)
    
    # Check key properties
    print("Screening function properties:")
    print(f"  S(0) = {screening_values[0]:.3f} (should be 0)")
    print(f"  S(∞) → {screening_values[-1]:.3f} (should approach 1)")
    
    # Find where screening = 0.5
    for i, s in enumerate(screening_values):
        if s > 0.5:
            u_half = u_values[i]
            break
    
    print(f"  S(u) = 0.5 at u ≈ {u_half:.2f}")
    print(f"  S(u) > 0.9 for u > {u_values[np.where(np.array(screening_values) > 0.9)[0][0]]:.1f}")
    
    # For SWM4-NDP with thole=1.3 and α=0.0013 nm³
    alpha = 0.0013
    thole = 1.3
    r_half = u_half / (thole / (alpha ** (1.0/6.0)))
    print(f"\nFor SWM4-NDP (thole={thole}, α={alpha} nm³):")
    print(f"  50% screening at r ≈ {r_half*10:.1f} Å")

if __name__ == "__main__":
    test_thole_water_screening()
    test_screening_formula()
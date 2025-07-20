#!/usr/bin/env python
"""Test Thole screening implementation"""

import sys
import numpy as np
sys.path.append('/home/zhaomt/gcmc/test107/pygcmc_dev/build/modules/bindings')

import pygcmc

def test_thole_screening():
    """Test Thole damped dipole-dipole interactions"""
    print("=== Testing Thole Screening Implementation ===\n")
    
    # Initialize
    pygcmc.initializeDrudeForce()
    
    # Create state
    state = pygcmc.MCState()
    state.info.cutoff = 2.0
    state.info.box = [10.0, 10.0, 10.0]
    
    # Force field
    state.forcefield.numTotalTypes = 2
    state.forcefield.numMovementTypes = 2
    
    # No LJ for this test
    state.forcefield.ljSigma = [0.0] * 4
    state.forcefield.ljEps = [0.0] * 4
    
    # Parameters
    charge_core = 1.0  # Core charge
    charge_drude = -1.0  # Drude charge (net neutral)
    polarizability = 0.001  # nm³
    thole = 2.5  # Typical Thole parameter
    
    print("Test parameters:")
    print(f"  Core charge: {charge_core}")
    print(f"  Drude charge: {charge_drude}")
    print(f"  Polarizability: {polarizability} nm³")
    print(f"  Thole parameter: {thole}")
    print()
    
    # Test at different distances
    distances = [0.3, 0.4, 0.5, 0.6, 0.8, 1.0]  # nm
    
    results = []
    
    for dist in distances:
        # Clear previous setup
        pygcmc.clearDrudeForce()
        pygcmc.initializeDrudeForce()
        
        # Create two dipoles
        atoms = []
        
        # Dipole 1
        core1 = pygcmc.MCAtom()
        core1.x = 0.0
        core1.y = 0.0
        core1.z = 0.0
        core1.charge = charge_core
        core1.type = 0
        atoms.append(core1)
        
        drude1 = pygcmc.MCAtom()
        drude1.x = 0.0
        drude1.y = 0.0
        drude1.z = 0.0
        drude1.charge = charge_drude
        drude1.type = 1
        atoms.append(drude1)
        
        # Dipole 2
        core2 = pygcmc.MCAtom()
        core2.x = dist
        core2.y = 0.0
        core2.z = 0.0
        core2.charge = charge_core
        core2.type = 0
        atoms.append(core2)
        
        drude2 = pygcmc.MCAtom()
        drude2.x = dist
        drude2.y = 0.0
        drude2.z = 0.0
        drude2.charge = charge_drude
        drude2.type = 1
        atoms.append(drude2)
        
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
            charge=charge_drude,
            polarizability=polarizability
        )
        
        pygcmc.addDrudeParticle(
            drudeIndex=3,
            parentIndex=2,
            charge=charge_drude,
            polarizability=polarizability
        )
        
        # Add screened pair
        pygcmc.addDrudeScreenedPair(
            dipole1=0,  # First Drude particle
            dipole2=1,  # Second Drude particle
            thole=thole
        )
        
        # Set SCF parameters
        params = pygcmc.DrudeSCFParams()
        params.tolerance = 1e-8
        params.maxIterations = 100
        params.dampingFactor = 0.5
        params.maxDrudeDistance = 0.02
        pygcmc.setDrudeSCFParameters(params)
        
        # Calculate energy with screening
        energy_screened, _ = pygcmc.computeSystemEnergyDrude(state)
        
        # Also calculate unscreened energy for comparison
        # This would be the standard Coulomb interaction
        pygcmc.computeSystemEnergyCutoff(state)
        coulomb_energy = sum(res.energy_elec for res in state.residues) / 2
        
        # Calculate screening factor
        u = dist * thole / (polarizability ** (1.0/3.0))
        screening_factor = 1.0 - (1.0 + 0.5 * u) * np.exp(-u)
        
        results.append({
            'distance': dist,
            'energy_screened': energy_screened,
            'coulomb_energy': coulomb_energy,
            'u': u,
            'screening': screening_factor
        })
        
        print(f"Distance: {dist*10:.1f} Å")
        print(f"  u parameter: {u:.3f}")
        print(f"  Screening factor: {screening_factor:.3f}")
        print(f"  Screened energy: {energy_screened:.3f} kJ/mol")
        print(f"  Coulomb energy: {coulomb_energy:.3f} kJ/mol")
        print()
    
    # Verify screening behavior
    print("\nVerification:")
    print("-" * 50)
    
    # At large distances, screening should approach 1
    if results[-1]['screening'] > 0.95:
        print("✓ Screening approaches 1 at large distances")
    else:
        print("✗ Screening should approach 1 at large distances")
    
    # At small distances, screening should be significant
    if results[0]['screening'] < 0.5:
        print("✓ Significant screening at small distances")
    else:
        print("✗ Should have significant screening at small distances")
    
    # Energy should be reduced by screening
    all_reduced = all(r['energy_screened'] < r['coulomb_energy'] * 0.9 
                     for r in results if r['distance'] < 0.5)
    if all_reduced:
        print("✓ Screened energy is reduced at short distances")
    else:
        print("✗ Screened energy should be reduced at short distances")

def test_thole_parameter_dependence():
    """Test how Thole parameter affects screening"""
    print("\n\n=== Testing Thole Parameter Dependence ===\n")
    
    # Initialize
    pygcmc.clearDrudeForce()
    pygcmc.initializeDrudeForce()
    
    # Fixed distance
    distance = 0.4  # nm
    polarizability = 0.001  # nm³
    
    # Test different Thole parameters
    thole_values = [0.0, 1.0, 2.0, 2.5, 3.0, 4.0]
    
    print(f"Fixed distance: {distance*10:.1f} Å")
    print(f"Polarizability: {polarizability} nm³")
    print("\nThole parameter effect:")
    print("-" * 40)
    
    for thole in thole_values:
        u = distance * thole / (polarizability ** (1.0/3.0))
        
        if thole == 0:
            screening = 0.0  # No screening
        else:
            screening = 1.0 - (1.0 + 0.5 * u) * np.exp(-u)
        
        print(f"Thole = {thole:.1f}: u = {u:.3f}, screening = {screening:.3f}")
    
    print("\nConclusion: Larger Thole parameter → stronger screening")

if __name__ == "__main__":
    test_thole_screening()
    test_thole_parameter_dependence()
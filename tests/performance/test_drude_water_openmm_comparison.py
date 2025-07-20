#!/usr/bin/env python
"""Compare PyGCMC Drude water model with OpenMM parameters and expected results"""

import sys
import numpy as np
import math
sys.path.append('/home/zhaomt/gcmc/test107/pygcmc_dev/build/modules/bindings')

import pygcmc

# OpenMM SWM4-NDP parameters from TestReferenceDrudeKernel.cpp
OPENMM_PARAMS = {
    'oxygen_mass': 15.6,  # 15.99943 - 0.4 (Drude mass)
    'drude_mass': 0.4,
    'drude_charge': -1.71636,
    'polarizability': 0.000978253,  # nm^3
    'thole': 1.3,
    'OH_bond': 0.09572,  # nm
    'HOH_angle': 104.52 * math.pi / 180,  # radians
    'hydrogen_charge': 0.52422,
    'oxygen_charge': 1.71636,  # Core charge before adding Drude
    'msite_weights': {
        'O': 0.786646558,
        'H1': 0.106676721,
        'H2': 0.106676721
    }
}

def create_swm4_water_openmm_style(x, y, z):
    """Create SWM4-NDP water using OpenMM parameters"""
    atoms = []
    
    # Oxygen (with reduced mass for Drude)
    o = pygcmc.MCAtom()
    o.x = x
    o.y = y
    o.z = z
    o.charge = OPENMM_PARAMS['oxygen_charge']
    o.type = 0
    atoms.append(o)
    
    # Drude on oxygen
    d = pygcmc.MCAtom()
    d.x = x
    d.y = y
    d.z = z
    d.charge = OPENMM_PARAMS['drude_charge']
    d.type = 1
    atoms.append(d)
    
    # Hydrogen 1
    h1 = pygcmc.MCAtom()
    h1.x = x + OPENMM_PARAMS['OH_bond'] * math.cos(OPENMM_PARAMS['HOH_angle']/2)
    h1.y = y + OPENMM_PARAMS['OH_bond'] * math.sin(OPENMM_PARAMS['HOH_angle']/2)
    h1.z = z
    h1.charge = OPENMM_PARAMS['hydrogen_charge']
    h1.type = 2
    atoms.append(h1)
    
    # Hydrogen 2
    h2 = pygcmc.MCAtom()
    h2.x = x + OPENMM_PARAMS['OH_bond'] * math.cos(OPENMM_PARAMS['HOH_angle']/2)
    h2.y = y - OPENMM_PARAMS['OH_bond'] * math.sin(OPENMM_PARAMS['HOH_angle']/2)
    h2.z = z
    h2.charge = OPENMM_PARAMS['hydrogen_charge']
    h2.type = 2
    atoms.append(h2)
    
    # M-site (virtual site)
    # Position based on weighted average
    weights = OPENMM_PARAMS['msite_weights']
    m = pygcmc.MCAtom()
    m.x = weights['O'] * o.x + weights['H1'] * h1.x + weights['H2'] * h2.x
    m.y = weights['O'] * o.y + weights['H1'] * h1.y + weights['H2'] * h2.y
    m.z = weights['O'] * o.z + weights['H1'] * h1.z + weights['H2'] * h2.z
    # M-site charge to make molecule neutral
    m.charge = -(o.charge + d.charge + h1.charge + h2.charge)
    m.type = 3
    atoms.append(m)
    
    # Verify charge neutrality
    total_charge = sum(atom.charge for atom in atoms)
    print(f"Water total charge: {total_charge:.6f} (should be ~0)")
    
    return atoms

def test_single_water_properties():
    """Test single water molecule properties against OpenMM"""
    print("=== Testing Single SWM4-NDP Water (OpenMM Parameters) ===\n")
    
    # Create state
    state = pygcmc.MCState()
    state.info.box = [10.0, 10.0, 10.0]
    state.info.cutoff = 5.0
    
    # Create water
    water_atoms = create_swm4_water_openmm_style(5.0, 5.0, 5.0)
    state.atoms = water_atoms
    state.activeAtomCount = 5
    
    # Create residue
    res = pygcmc.MCResidue()
    res.atomStart = 0
    res.atomCount = 5
    res.active = True
    res.type = 0
    state.residues = [res]
    state.activeResidueCount = 1
    
    # Set force field (no LJ for this test)
    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 4
    ff.numMovementTypes = 4
    ff.ljSigma = [0.0] * 16
    ff.ljEps = [0.0] * 16
    state.forcefield = ff
    
    # Create Drude force
    drude_force = pygcmc.DrudeForce()
    
    # Add Drude particle with OpenMM parameters
    drude_idx = drude_force.addParticle(
        drudeIndex=1,
        parentIndex=0,
        aniso1Index=-1,
        aniso2Index=-1,
        aniso3Index=-1,
        aniso4Index=-1,
        charge=OPENMM_PARAMS['drude_charge'],
        polarizability=OPENMM_PARAMS['polarizability'],
        aniso12=1.0,
        aniso34=1.0
    )
    
    # Set SCF parameters similar to OpenMM
    scf_params = pygcmc.DrudeSCFParams()
    scf_params.tolerance = 1e-8
    scf_params.maxIterations = 100
    scf_params.maxDrudeDistance = 0.02  # OpenMM default
    drude_force.setSCFParameters(scf_params)
    
    # Calculate energy
    energy = drude_force.calculateEnergySCF(state)
    
    print(f"Single water energy: {energy:.6f} kJ/mol")
    print(f"Expected: 0.0 kJ/mol (no intramolecular interactions)\n")
    
    # Check Drude position
    dx = state.atoms[1].x - state.atoms[0].x
    dy = state.atoms[1].y - state.atoms[0].y
    dz = state.atoms[1].z - state.atoms[0].z
    dist = math.sqrt(dx*dx + dy*dy + dz*dz)
    print(f"Drude displacement: {dist*1000:.3f} pm")
    print(f"Expected: < 0.1 pm (no external field)\n")
    
    # Calculate dipole moment
    dipole_x = dipole_y = dipole_z = 0.0
    com_x = com_y = com_z = 0.0
    total_mass = 15.99943 + 0.4 + 1.00783*2  # O + D + 2H
    
    # Calculate center of mass (excluding virtual site)
    masses = [15.6, 0.4, 1.00783, 1.00783, 0]  # O, D, H, H, M
    for i in range(4):  # Exclude M-site
        com_x += masses[i] * state.atoms[i].x
        com_y += masses[i] * state.atoms[i].y
        com_z += masses[i] * state.atoms[i].z
    com_x /= total_mass
    com_y /= total_mass
    com_z /= total_mass
    
    # Calculate dipole
    for atom in state.atoms:
        dipole_x += atom.charge * (atom.x - com_x)
        dipole_y += atom.charge * (atom.y - com_y)
        dipole_z += atom.charge * (atom.z - com_z)
    
    # Convert to Debye (e*nm to Debye: multiply by 48.0321)
    dipole_mag = math.sqrt(dipole_x**2 + dipole_y**2 + dipole_z**2)
    dipole_debye = dipole_mag * 48.0321
    
    print(f"Water dipole moment: {dipole_debye:.3f} D")
    print(f"Expected (SWM4-NDP): ~2.4 D\n")

def test_water_dimer_interaction():
    """Test water dimer interaction with OpenMM parameters"""
    print("=== Testing Water Dimer Interaction ===\n")
    
    # Create state
    state = pygcmc.MCState()
    state.info.box = [10.0, 10.0, 10.0]
    state.info.cutoff = 5.0
    
    # Create two waters separated by 0.28 nm (O-O distance)
    water1_atoms = create_swm4_water_openmm_style(5.0, 5.0, 5.0)
    water2_atoms = create_swm4_water_openmm_style(5.28, 5.0, 5.0)
    
    state.atoms = water1_atoms + water2_atoms
    state.activeAtomCount = 10
    
    # Create residues
    res1 = pygcmc.MCResidue()
    res1.atomStart = 0
    res1.atomCount = 5
    res1.active = True
    res1.type = 0
    
    res2 = pygcmc.MCResidue()
    res2.atomStart = 5
    res2.atomCount = 5
    res2.active = True
    res2.type = 0
    
    state.residues = [res1, res2]
    state.activeResidueCount = 2
    
    # Set force field with LJ parameters
    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 4
    ff.numMovementTypes = 4
    
    # OpenMM uses sigma=0.318395 nm, epsilon=0.88257 kJ/mol for O-O
    ljSigma = [0.0] * 16
    ljEps = [0.0] * 16
    
    # Only O-O interaction
    ljSigma[0] = 0.318395  # O-O sigma
    ljEps[0] = 0.88257     # O-O epsilon
    
    ff.ljSigma = ljSigma
    ff.ljEps = ljEps
    state.forcefield = ff
    
    # Create Drude force
    drude_force = pygcmc.DrudeForce()
    
    # Add Drude particles for both waters
    for i in [0, 1]:
        water_offset = i * 5
        drude_force.addParticle(
            drudeIndex=water_offset + 1,
            parentIndex=water_offset + 0,
            aniso1Index=-1,
            aniso2Index=-1,
            aniso3Index=-1,
            aniso4Index=-1,
            charge=OPENMM_PARAMS['drude_charge'],
            polarizability=OPENMM_PARAMS['polarizability'],
            aniso12=1.0,
            aniso34=1.0
        )
    
    # Add Thole screening between the two Drude particles
    drude_force.addScreenedPair(0, 1, OPENMM_PARAMS['thole'])
    
    # Calculate energy
    energy_drude = drude_force.calculateEnergySCF(state)
    
    # Calculate LJ energy
    pygcmc.computeSystemEnergyCutoff(state)
    energy_lj = sum(res.energy_vdw for res in state.residues) / 2  # Avoid double counting
    
    # Calculate Coulomb energy
    energy_coulomb = sum(res.energy_elec for res in state.residues) / 2
    
    total_energy = energy_lj + energy_coulomb + energy_drude
    
    print(f"O-O distance: 0.28 nm")
    print(f"LJ energy: {energy_lj:.3f} kJ/mol")
    print(f"Coulomb energy: {energy_coulomb:.3f} kJ/mol")
    print(f"Drude energy: {energy_drude:.3f} kJ/mol")
    print(f"Total energy: {total_energy:.3f} kJ/mol")
    print(f"Expected: ~-20 to -25 kJ/mol for hydrogen-bonded dimer\n")
    
    # Check Drude displacements
    for i in [0, 1]:
        water_offset = i * 5
        dx = state.atoms[water_offset + 1].x - state.atoms[water_offset + 0].x
        dy = state.atoms[water_offset + 1].y - state.atoms[water_offset + 0].y
        dz = state.atoms[water_offset + 1].z - state.atoms[water_offset + 0].z
        dist = math.sqrt(dx*dx + dy*dy + dz*dz)
        print(f"Water {i+1} Drude displacement: {dist*1000:.3f} pm")

def test_force_constant_calculation():
    """Verify force constant calculation matches OpenMM"""
    print("\n=== Testing Force Constant Calculation ===\n")
    
    # OpenMM parameters
    charge = OPENMM_PARAMS['drude_charge']
    pol = OPENMM_PARAMS['polarizability']
    
    # Calculate force constant
    ONE_4PI_EPS0 = 138.935456  # kJ/mol·nm·e^-2
    k_calc = ONE_4PI_EPS0 * charge * charge / pol
    
    print(f"Drude charge: {charge} e")
    print(f"Polarizability: {pol} nm³ = {pol*1000} Å³")
    print(f"Calculated k: {k_calc:.1f} kJ/mol/nm²")
    print(f"In kcal/mol/Å²: {k_calc / 4.184 / 100:.1f}")
    
    # OpenMM's expected value
    k_expected = 100000.0  # kcal/mol/Å²
    k_expected_kj = k_expected * 4.184 * 100  # kJ/mol/nm²
    
    print(f"\nExpected k: {k_expected} kcal/mol/Å²")
    print(f"Expected k: {k_expected_kj:.1f} kJ/mol/nm²")
    print(f"Ratio: {k_calc/k_expected_kj:.3f}")

def test_thole_screening_effect():
    """Test Thole screening between water molecules"""
    print("\n=== Testing Thole Screening Effect ===\n")
    
    distances = [0.25, 0.30, 0.35, 0.40, 0.50]  # nm
    
    for dist in distances:
        # Calculate Thole parameter
        pol = OPENMM_PARAMS['polarizability']
        thole = OPENMM_PARAMS['thole']
        
        # u = r * thole / (pol^(1/3))
        # Note: OpenMM uses (pol1*pol2)^(1/6) for mixed polarizabilities
        u = dist * thole / (pol ** (1.0/3.0))
        
        # Screening function: S(u) = 1 - (1 + u/2) * exp(-u)
        if u > 0:
            screening = 1.0 - (1.0 + 0.5 * u) * np.exp(-u)
        else:
            screening = 0.0
        
        print(f"Distance: {dist*10:.1f} Å")
        print(f"  u parameter: {u:.3f}")
        print(f"  Screening factor: {screening:.3f}")
        print(f"  Effective interaction: {(1-screening)*100:.1f}% of unscreened")
        print()

if __name__ == "__main__":
    test_single_water_properties()
    test_water_dimer_interaction()
    test_force_constant_calculation()
    test_thole_screening_effect()
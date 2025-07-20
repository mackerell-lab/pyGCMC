#!/usr/bin/env python
"""Debug OPT3 perturbation calculation issues"""

import sys
import numpy as np
sys.path.append('/home/zhaomt/gcmc/test107/pygcmc_dev/build/modules/bindings')

import pygcmc

def create_simple_water_system(n_waters=2):
    """Create a simple water system for debugging"""
    state = pygcmc.MCState()
    state.info.box = [3.0, 3.0, 3.0]
    state.info.cutoff = 1.5
    
    atoms = []
    
    # SWM4-NDP water parameters
    charges = [1.71636, -1.71636, 0.55733, 0.55733, -1.11466]
    types = [0, 1, 2, 2, 3]
    
    # Place waters in a line with 0.3 nm spacing
    for i in range(n_waters):
        x_offset = 1.5 + i * 0.3
        
        positions = [
            [x_offset, 1.5, 1.5],          # O
            [x_offset, 1.5, 1.5],          # D (initially at parent)
            [x_offset + 0.09572, 1.5, 1.5], # H1
            [x_offset - 0.03, 1.59, 1.5],   # H2
            [x_offset + 0.015, 1.511, 1.5]  # M-site
        ]
        
        for j in range(5):
            a = pygcmc.MCAtom()
            a.x, a.y, a.z = positions[j]
            a.charge = charges[j]
            a.type = types[j]
            atoms.append(a)
    
    state.atoms = atoms
    state.activeAtomCount = len(atoms)
    
    # Create residues
    residues = []
    for i in range(n_waters):
        res = pygcmc.MCResidue()
        res.atomStart = 5 * i
        res.atomCount = 5
        res.active = True
        res.type = 0
        residues.append(res)
    
    state.residues = residues
    state.activeResidueCount = len(residues)
    
    # Set force field
    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 4
    ff.numMovementTypes = 4
    ff.ljSigma = [0.318395] + [0.0] * 15
    ff.ljEps = [0.88257] + [0.0] * 15
    state.forcefield = ff
    
    return state

def analyze_electric_field(state, drude_force):
    """Analyze electric field at Drude positions"""
    print("=== Electric Field Analysis ===\n")
    
    # Manual calculation of electric field at first Drude
    drude_idx = 1  # First Drude particle
    drude_pos = [state.atoms[drude_idx].x, 
                 state.atoms[drude_idx].y, 
                 state.atoms[drude_idx].z]
    
    print(f"Drude position: {drude_pos}")
    
    # Calculate field from each charge
    E_total = [0.0, 0.0, 0.0]
    ONE_4PI_EPS0 = 138.935456  # kJ/mol·nm·e^-2
    
    print("\nContributions to electric field:")
    for i, atom in enumerate(state.atoms):
        if i == drude_idx:
            continue  # Skip self
            
        # Skip other Drude particles for static field
        if atom.type == 1:  # Drude type
            continue
            
        dx = atom.x - drude_pos[0]
        dy = atom.y - drude_pos[1]
        dz = atom.z - drude_pos[2]
        
        # Apply PBC
        if state.info.box[0] > 0:
            dx -= state.info.box[0] * round(dx / state.info.box[0])
            dy -= state.info.box[1] * round(dy / state.info.box[1])
            dz -= state.info.box[2] * round(dz / state.info.box[2])
        
        r2 = dx*dx + dy*dy + dz*dz
        r = np.sqrt(r2)
        
        if r < 1e-10:
            continue
            
        # E = k * q / r^2 * r_hat
        E_mag = ONE_4PI_EPS0 * atom.charge / r2
        Ex = E_mag * dx / r
        Ey = E_mag * dy / r
        Ez = E_mag * dz / r
        
        E_total[0] += Ex
        E_total[1] += Ey
        E_total[2] += Ez
        
        if abs(atom.charge) > 0.01:  # Only print significant charges
            print(f"  Atom {i} (type {atom.type}, q={atom.charge:+.3f}): "
                  f"r={r:.4f} nm, E=({Ex:+.2e}, {Ey:+.2e}, {Ez:+.2e})")
    
    E_mag = np.sqrt(E_total[0]**2 + E_total[1]**2 + E_total[2]**2)
    print(f"\nTotal electric field: ({E_total[0]:.2e}, {E_total[1]:.2e}, {E_total[2]:.2e})")
    print(f"Magnitude: {E_mag:.2e} kJ/(mol·nm·e)")
    
    # Expected Drude displacement
    q_drude = -1.71636
    k_drude = drude_force.getOPT3Coefficients().c0  # Just to access the force constant
    # Actually get it from the calculation
    polarizability = 0.000978253
    k_drude = ONE_4PI_EPS0 * q_drude * q_drude / polarizability
    
    print(f"\nDrude parameters:")
    print(f"  Charge: {q_drude}")
    print(f"  Force constant: {k_drude:.1f} kJ/mol/nm²")
    print(f"  Polarizability: {polarizability:.6f} nm³")
    
    displacement = E_mag * abs(q_drude) / k_drude
    print(f"\nExpected displacement: {displacement:.6f} nm ({displacement*10:.4f} Å)")

def test_perturbation_orders():
    """Test the perturbation order calculations"""
    print("=== Testing Perturbation Orders ===\n")
    
    # Create simple 2-water system
    state = create_simple_water_system(2)
    
    # Create DrudeForce
    drude_force = pygcmc.DrudeForce()
    
    # Add Drude particles
    for i in range(2):
        drude_force.addParticle(
            drudeIndex=5*i + 1,
            parentIndex=5*i,
            aniso1Index=-1, aniso2Index=-1,
            aniso3Index=-1, aniso4Index=-1,
            charge=-1.71636,
            polarizability=0.000978253,
            aniso12=0.0, aniso34=0.0
        )
    
    # Add screening between the two waters
    drude_force.addScreenedPair(0, 1, 1.3)
    
    # First analyze the electric field
    analyze_electric_field(state, drude_force)
    
    # Run SCF to see what the actual converged position is
    print("\n\n=== SCF Convergence ===")
    drude_force.setUseOPT3(False)
    energy_scf = drude_force.calculateEnergySCF(state)
    print(f"SCF Energy: {energy_scf:.6f} kJ/mol")
    
    # Check Drude positions after SCF
    print("\nDrude positions after SCF:")
    for i in range(2):
        drude_idx = 5*i + 1
        parent_idx = 5*i
        dx = state.atoms[drude_idx].x - state.atoms[parent_idx].x
        dy = state.atoms[drude_idx].y - state.atoms[parent_idx].y
        dz = state.atoms[drude_idx].z - state.atoms[parent_idx].z
        displacement = np.sqrt(dx*dx + dy*dy + dz*dz)
        print(f"  Drude {i}: displacement = {displacement:.6f} nm ({displacement*10:.4f} Å)")
    
    # Now collect training data
    print("\n\n=== Collecting Training Data ===")
    training_data = drude_force.collectTrainingData(state)
    
    print(f"\nNumber of Drude particles: {len(training_data.r0)}")
    
    for i in range(len(training_data.r0)):
        print(f"\nDrude particle {i}:")
        r0 = training_data.r0[i]
        r1 = training_data.r1[i]
        r2 = training_data.r2[i]
        r3 = training_data.r3[i]
        r_scf = training_data.r_scf[i]
        
        print(f"  |r0|   = {r0.norm():.6f} nm ({r0.norm()*10:.4f} Å)")
        print(f"  |r1|   = {r1.norm():.6f} nm ({r1.norm()*10:.4f} Å)")
        print(f"  |r2|   = {r2.norm():.6f} nm ({r2.norm()*10:.4f} Å)")
        print(f"  |r3|   = {r3.norm():.6f} nm ({r3.norm()*10:.4f} Å)")
        print(f"  |r_scf| = {r_scf.norm():.6f} nm ({r_scf.norm()*10:.4f} Å)")
        
        # Check if orders are decreasing
        print(f"\n  Convergence check:")
        print(f"    |r1|/|r0| = {r1.norm()/r0.norm() if r0.norm() > 0 else 'inf':.4f}")
        print(f"    |r2|/|r1| = {r2.norm()/r1.norm() if r1.norm() > 0 else 'inf':.4f}")
        print(f"    |r3|/|r2| = {r3.norm()/r2.norm() if r2.norm() > 0 else 'inf':.4f}")

def main():
    """Run debugging tests"""
    test_perturbation_orders()
    
    print("\n\n=== Analysis ===")
    print("\nThe issue is that r0 is way too large (~18 Å instead of ~0.01-0.1 Å).")
    print("This suggests the electric field calculation or force constant is wrong.")
    print("\nPossible causes:")
    print("1. Electric field units mismatch")
    print("2. Force constant calculation error")
    print("3. Missing residue exclusions in field calculation")

if __name__ == "__main__":
    main()
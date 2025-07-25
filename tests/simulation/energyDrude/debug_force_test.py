"""
Debug test to check force calculation issue
"""

import pygcmc
import math


def test_simple_drude_force():
    """Simple test to debug force issue"""
    
    # Minimal system
    state = pygcmc.MCState()
    state.info.box = [10.0, 10.0, 10.0]
    state.info.cutoff = 4.0  # Ensure cutoff is set
    
    # Force field
    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 3  # Parent, Drude, External
    ff.numMovementTypes = 3
    ff.ljEps = [0.0, 0.0, 0.0]
    ff.ljSigma = [0.1, 0.1, 0.1]
    state.forcefield = ff
    
    # Just parent and Drude
    parent = pygcmc.MCAtom()
    parent.x = parent.y = parent.z = 0.0
    parent.charge = 1.0
    parent.type = 0
    
    drude = pygcmc.MCAtom()
    drude.x = drude.y = drude.z = 0.0  # Start at parent position
    drude.charge = -1.0
    drude.type = 1
    
    state.atoms = [parent, drude]
    state.activeAtomCount = 2
    
    # Single residue
    res = pygcmc.MCResidue()
    res.atomStart = 0
    res.atomCount = 2
    res.active = True
    res.type = 0
    state.residues = [res]
    state.activeResidueCount = 1
    
    # Setup Drude
    pygcmc.DrudeComplete.clear()
    
    particle = pygcmc.DrudeParticle()
    particle.drudeIndex = 1
    particle.parentIndex = 0
    particle.charge = -1.0
    particle.polarizability = 0.001  # nm^3
    particle.computeSpringConstants()
    
    print(f"Spring constant: {particle.kSpring:.6f} kJ/mol/nm^2")
    print(f"Polarizability: {particle.polarizability:.6f} nm^3")
    print(f"Drude charge: {particle.charge:.6f}")
    print(f"Parent index: {particle.parentIndex}, Drude index: {particle.drudeIndex}")
    
    pygcmc.DrudeComplete.addParticle(particle)
    print(f"Number of Drude particles: {pygcmc.DrudeComplete.getNumParticles()}")
    
    # SCF parameters
    params = pygcmc.DrudeSCFParams()
    params.tolerance = 0.1
    params.maxIterations = 100
    params.enableHardWall = False
    params.dampingFactor = 0.5
    pygcmc.DrudeComplete.setParameters(params)
    
    # Algorithm is already SCF by default
    
    # Test 1: No external field
    print("\nTest 1: No external field")
    print(f"Initial Drude position: ({state.atoms[1].x:.6f}, {state.atoms[1].y:.6f}, {state.atoms[1].z:.6f})")
    energy = pygcmc.DrudeComplete.calculateEnergy(state)
    print(f"Energy: {energy:.10f} kJ/mol")
    print(f"Final Drude position: ({state.atoms[1].x:.6f}, {state.atoms[1].y:.6f}, {state.atoms[1].z:.6f})")
    
    # Manually displace Drude to check energy calculation
    print("\nManually displacing Drude by 0.01 nm in z:")
    state.atoms[1].z = 0.01
    energy_displaced = pygcmc.DrudeComplete.calculateEnergy(state)
    print(f"Energy with displacement: {energy_displaced:.6f} kJ/mol")
    expected_energy = 0.5 * particle.kSpring * 0.01 * 0.01
    print(f"Expected harmonic energy: {expected_energy:.6f} kJ/mol")
    
    # Reset for next test
    state.atoms[1].z = 0.0
    
    # Test 2: With external charge
    print("\nTest 2: Adding external charge")
    external = pygcmc.MCAtom()
    external.x = 2.0
    external.y = external.z = 0.0
    external.charge = 1.0
    external.type = 0
    state.atoms.append(external)
    state.activeAtomCount = 3
    
    res_ext = pygcmc.MCResidue()
    res_ext.atomStart = 2
    res_ext.atomCount = 1
    res_ext.active = True
    res_ext.type = 1
    state.residues.append(res_ext)
    state.activeResidueCount = 2
    
    energy2 = pygcmc.DrudeComplete.calculateEnergy(state)
    print(f"Energy with external: {energy2:.6f} kJ/mol")
    print(f"Drude position: ({state.atoms[1].x:.6f}, {state.atoms[1].y:.6f}, {state.atoms[1].z:.6f})")
    
    # Test 3: Manual force check
    print("\nTest 3: Checking forces")
    # The Drude should have moved towards the external positive charge
    drude_disp = state.atoms[1].x - state.atoms[0].x
    print(f"Drude displacement in x: {drude_disp:.6f} nm")
    
    # Estimate force on Drude from spring
    spring_force = -particle.kSpring * drude_disp
    print(f"Spring force on Drude: {spring_force:.6f} kJ/mol/nm")
    
    # Estimate electrostatic force from external charge
    r = 2.0 - state.atoms[1].x  # Distance from Drude to external
    if abs(r) > 1e-6:
        elec_force = 138.935456 * (-1.0) * 1.0 / (r * r)  # ONE_4PI_EPS0 * q1 * q2 / r^2
        print(f"Electrostatic force from external: {elec_force:.6f} kJ/mol/nm")
        print(f"Net force (should be ~0): {spring_force + elec_force:.6f} kJ/mol/nm")
    
    pygcmc.DrudeComplete.clear()


if __name__ == "__main__":
    test_simple_drude_force()
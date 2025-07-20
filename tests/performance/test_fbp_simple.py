#!/usr/bin/env python3
"""
Simple test of Force Balance Predictor (FBP) algorithm
"""

import numpy as np
import pygcmc
import time

def create_water_system(n_waters=10):
    """Create a simple water system"""
    atoms = []
    residues = []
    
    # Simple grid placement
    spacing = 0.4  # nm
    
    for i in range(n_waters):
        # Base position
        x = (i % 5) * spacing
        y = (i // 5) * spacing
        z = 0.0
        
        # Oxygen
        atom = pygcmc.MCAtom()
        atom.x = x
        atom.y = y
        atom.z = z
        atom.charge = 1.71636
        atom.type = 0
        atoms.append(atom)
        
        # Drude on oxygen
        drude = pygcmc.MCAtom()
        drude.x = x
        drude.y = y
        drude.z = z
        drude.charge = -1.71636
        drude.type = 1
        atoms.append(drude)
        
        # Hydrogen 1
        h1 = pygcmc.MCAtom()
        h1.x = x + 0.09572
        h1.y = y
        h1.z = z
        h1.charge = 0.55733
        h1.type = 2
        atoms.append(h1)
        
        # Hydrogen 2
        h2 = pygcmc.MCAtom()
        h2.x = x - 0.04786
        h2.y = y + 0.08288
        h2.z = z
        h2.charge = 0.55733
        h2.type = 2
        atoms.append(h2)
        
        # M-site
        m = pygcmc.MCAtom()
        m.x = x
        m.y = y - 0.024034
        m.z = z
        m.charge = -1.11466
        m.type = 3
        atoms.append(m)
        
        # Create residue
        res = pygcmc.MCResidue()
        res.atomStart = 5 * i
        res.atomCount = 5
        res.active = True
        res.type = 0
        residues.append(res)
    
    # Create state
    state = pygcmc.MCState()
    state.atoms = atoms
    state.residues = residues
    state.activeAtomCount = len(atoms)
    state.activeResidueCount = len(residues)
    state.activeTypeCount = 4
    
    # Box info
    state.info.box = np.array([10.0, 10.0, 10.0])
    state.info.cutoff = 4.5
    
    # Force field
    state.forcefield.numTotalTypes = 4
    state.forcefield.numMovementTypes = 4
    state.forcefield.ljSigma = [0.318395, 0.0, 0.0, 0.0]  # Only oxygen has LJ
    state.forcefield.ljEps = [0.88257, 0.0, 0.0, 0.0]
    
    return state, n_waters

def test_fbp_algorithm():
    """Test FBP algorithm performance"""
    print("Testing Force Balance Predictor Algorithm")
    print("="*60)
    
    # Test different system sizes
    sizes = [10, 20, 40, 80]
    
    for n_waters in sizes:
        print(f"\n{n_waters} waters:")
        
        # Create system
        state, _ = create_water_system(n_waters)
        
        # Initialize Drude force
        pygcmc.initializeDrudeForce()
        
        # SWM4-NDP parameters
        charge = -1.71636
        ONE_4PI_EPS0 = 138.935456
        k = 418400.0  # kJ/mol/nm^2
        polarizability = ONE_4PI_EPS0 * charge * charge / k
        thole = 1.3
        
        # Setup force objects for each algorithm
        algorithms = [
            (pygcmc.DrudeAlgorithm.SCF, "SCF"),
            (pygcmc.DrudeAlgorithm.OPT3, "OPT3"), 
            (pygcmc.DrudeAlgorithm.SmartOPT3, "Smart OPT3"),
            (pygcmc.DrudeAlgorithm.FBP, "FBP")
        ]
        
        ref_energy = None
        
        for algo, name in algorithms:
            # Create force object
            force = pygcmc.DrudeForce()
            
            # Add Drude particles
            for i in range(n_waters):
                force.addParticle(
                    drudeIndex=5*i + 1,
                    parentIndex=5*i,
                    aniso1Index=-1, aniso2Index=-1,
                    aniso3Index=-1, aniso4Index=-1,
                    charge=charge,
                    polarizability=polarizability,
                    aniso12=1.0, aniso34=1.0
                )
            
            # Add screened pairs (all pairs)
            for i in range(n_waters):
                for j in range(i+1, n_waters):
                    force.addScreenedPair(i, j, thole)
            
            # Set parameters
            params = pygcmc.DrudeSCFParams()
            if algo == pygcmc.DrudeAlgorithm.SCF and ref_energy is None:
                # Tight convergence for reference
                params.tolerance = 0.01
                params.maxIterations = 200
            else:
                # Standard parameters
                params.tolerance = 1.0
                params.maxIterations = 50
                
            params.dampingFactor = 0.5
            params.forceCutoff = 10.0
            params.maxDrudeDistance = 0.02
            force.setSCFParameters(params)
            force.setAlgorithm(algo)
            
            # Reset Drude positions to parent positions
            for i in range(n_waters):
                state.atoms[5*i + 1].x = state.atoms[5*i].x
                state.atoms[5*i + 1].y = state.atoms[5*i].y
                state.atoms[5*i + 1].z = state.atoms[5*i].z
            
            # Time the calculation
            start = time.time()
            energy = force.calculateEnergySCF(state)
            elapsed = time.time() - start
            
            if ref_energy is None:
                ref_energy = energy
                
            error = abs(energy - ref_energy) / abs(ref_energy) * 100.0 if ref_energy != 0 else 0.0
            
            print(f"  {name:<12}: Time={elapsed*1000:6.2f}ms, Energy={energy:10.2f}, Error={error:6.3f}%")

if __name__ == "__main__":
    test_fbp_algorithm()
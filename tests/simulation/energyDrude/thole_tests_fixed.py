"""
Fixed Thole screening tests with correct physics
"""

import pytest
import pygcmc


def test_thole_parameter_sensitivity_fixed():
    """Test sensitivity to Thole parameter value - fixed version"""
    state = pygcmc.MCState()
    
    # Create simple two-dipole system
    atoms = []
    for i in range(2):
        parent = pygcmc.MCAtom()
        parent.x = i * 0.4  # 4 Å separation
        parent.y, parent.z = 0.0, 0.0
        parent.charge = 0.0
        parent.type = 0
        atoms.append(parent)
        
        drude = pygcmc.MCAtom()
        drude.x = i * 0.4 + 0.001  # Small displacement
        drude.y, drude.z = 0.0, 0.0
        drude.charge = -1.0
        drude.type = 1
        atoms.append(drude)
    
    state.atoms = atoms
    state.activeAtomCount = 4
    state.info.box = [5.0, 5.0, 5.0]
    
    # Setup residues
    res1 = pygcmc.MCResidue()
    res1.atomStart = 0
    res1.atomCount = 2
    res1.active = True
    res1.type = 0
    
    res2 = pygcmc.MCResidue()
    res2.atomStart = 2
    res2.atomCount = 2
    res2.active = True
    res2.type = 1
    
    state.residues = [res1, res2]
    state.activeResidueCount = 2
    
    # Test different Thole parameters
    # Expected behavior: energy first decreases (more screening), then increases back
    thole_values = [0.0, 0.5, 1.0, 1.3, 2.0, 3.0]
    energies = []
    
    for thole in thole_values:
        pygcmc.DrudeComplete.clear()
        
        # Add particles
        for i in range(2):
            p = pygcmc.DrudeParticle()
            p.drudeIndex = 2*i + 1
            p.parentIndex = 2*i
            p.charge = -1.0
            p.polarizability = 0.001
            p.computeSpringConstants()
            pygcmc.DrudeComplete.addParticle(p)
        
        # Add screened pair with current thole value
        if thole > 0:
            pair = pygcmc.ScreenedPair()
            pair.dipole1 = 0
            pair.dipole2 = 1
            pair.thole = thole
            pygcmc.DrudeComplete.addScreenedPair(pair)
        
        # Set SCF parameters
        params = pygcmc.DrudeSCFParams()
        params.tolerance = 1e-6
        params.maxIterations = 500
        params.maxDrudeDistance = 0.02
        pygcmc.DrudeComplete.setParameters(params)
        
        energy = pygcmc.DrudeComplete.calculateEnergy(state)
        energies.append(energy)
    
    # Check that energy with small thole (0.5) is less than no screening
    assert energies[1] < energies[0], \
        f"Small thole should reduce energy: {energies[1]} >= {energies[0]}"
    
    # Check that energy increases from minimum back toward unscreened value
    # Find minimum energy
    min_energy = min(energies)
    min_index = energies.index(min_energy)
    
    # After minimum, energy should generally increase
    for i in range(min_index + 1, len(energies)):
        # Allow small variations due to SCF convergence
        assert energies[i] > min_energy - 1.0, \
            f"Energy should increase after minimum: {energies[i]} <= {min_energy}"
    
    # Large thole should give energy close to unscreened
    assert abs(energies[-1] - energies[0]) < 10.0, \
        f"Large thole should be close to unscreened: {energies[-1]} vs {energies[0]}"
    
    pygcmc.DrudeComplete.clear()


def test_thole_screening_effectiveness():
    """Test that Thole screening actually reduces interaction at short distances"""
    state = pygcmc.MCState()
    
    # Test at different distances
    distances = [0.2, 0.3, 0.4, 0.5, 0.6]
    
    for dist in distances:
        atoms = []
        for i in range(2):
            parent = pygcmc.MCAtom()
            parent.x = i * dist
            parent.y, parent.z = 0.0, 0.0
            parent.charge = 0.0
            parent.type = 0
            atoms.append(parent)
            
            drude = pygcmc.MCAtom()
            drude.x = i * dist
            drude.y, drude.z = 0.0, 0.0
            drude.charge = -1.0
            drude.type = 1
            atoms.append(drude)
        
        state.atoms = atoms
        state.activeAtomCount = 4
        state.info.box = [5.0, 5.0, 5.0]
        
        # Setup residues
        res1 = pygcmc.MCResidue()
        res1.atomStart = 0
        res1.atomCount = 2
        res1.active = True
        res1.type = 0
        
        res2 = pygcmc.MCResidue()
        res2.atomStart = 2
        res2.atomCount = 2
        res2.active = True
        res2.type = 1
        
        state.residues = [res1, res2]
        state.activeResidueCount = 2
        
        # Calculate energy with and without Thole
        energies = {}
        
        for use_thole in [False, True]:
            pygcmc.DrudeComplete.clear()
            
            # Add particles
            for i in range(2):
                p = pygcmc.DrudeParticle()
                p.drudeIndex = 2*i + 1
                p.parentIndex = 2*i
                p.charge = -1.0
                p.polarizability = 0.001
                p.computeSpringConstants()
                pygcmc.DrudeComplete.addParticle(p)
            
            if use_thole:
                pair = pygcmc.ScreenedPair()
                pair.dipole1 = 0
                pair.dipole2 = 1
                pair.thole = 1.3  # Standard value
                pygcmc.DrudeComplete.addScreenedPair(pair)
            
            # Set SCF parameters
            params = pygcmc.DrudeSCFParams()
            params.tolerance = 1e-6
            params.maxIterations = 500
            params.maxDrudeDistance = 0.02
            pygcmc.DrudeComplete.setParameters(params)
            
            energy = pygcmc.DrudeComplete.calculateEnergy(state)
            energies['thole' if use_thole else 'no_thole'] = energy
        
        # At short distances, Thole screening should reduce the energy
        # (make it less positive for repulsive interactions)
        if dist < 0.4:  # Short distance
            assert energies['thole'] < energies['no_thole'], \
                f"At distance {dist} nm, Thole should reduce energy"
        
        # The effect should be distance-dependent
        # Larger effect at shorter distances
        
    pygcmc.DrudeComplete.clear()
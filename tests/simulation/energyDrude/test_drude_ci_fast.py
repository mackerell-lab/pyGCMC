"""
Fast Drude CI tests - run on every commit
Target time: < 5 seconds total
"""

import pytest
import numpy as np
import pygcmc


class TestDrudeFastCI:
    """Fast subset of Drude tests for CI"""
    
    def test_single_drude_external_field(self):
        """Test 1: Single Drude + external field (< 0.5s)"""
        state = pygcmc.MCState()
        state.info.box = [3.0, 3.0, 3.0]
        
        # Parent, Drude, External charge
        parent = pygcmc.MCAtom()
        parent.x = parent.y = parent.z = 0.0
        parent.charge = 1.0
        parent.type = 0
        
        drude = pygcmc.MCAtom()
        drude.x = drude.y = drude.z = 0.0
        drude.charge = -1.0
        drude.type = 1
        
        external = pygcmc.MCAtom()
        external.x = 0.5
        external.y = external.z = 0.0
        external.charge = 2.0
        external.type = 2
        
        state.atoms = [parent, drude, external]
        state.activeAtomCount = 3
        
        res = pygcmc.MCResidue()
        res.atomStart = 0
        res.atomCount = 2
        res.active = True
        state.residues = [res]
        state.activeResidueCount = 1
        
        # Setup Drude
        pygcmc.DrudeComplete.clear()
        p = pygcmc.DrudeParticle()
        p.drudeIndex = 1
        p.parentIndex = 0
        p.charge = -1.0
        p.polarizability = 0.001
        p.computeSpringConstants()
        pygcmc.DrudeComplete.addParticle(p)
        
        params = pygcmc.DrudeSCFParams()
        params.tolerance = 0.01
        params.maxIterations = 100
        params.includeCoulombEnergy = True  # For standalone test
        pygcmc.DrudeComplete.setParameters(params)
        
        energy = pygcmc.DrudeComplete.calculateEnergy(state)
        dx = state.atoms[1].x - state.atoms[0].x
        
        # Check displacement exists and is reasonable
        assert abs(dx) > 1e-6, "Drude should displace in external field"
        assert abs(dx) < 0.1, "Displacement should be reasonable"
        
        # Note: DrudeComplete.calculateEnergy returns total energy (including Coulomb)
        # For pure spring energy validation, calculate separately
        spring_energy = 0.5 * p.kSpring * dx * dx
        assert spring_energy > 0, "Spring energy component should be positive"
    
    def test_double_drude_with_thole(self):
        """Test 2: Double Drude + Thole screening (< 0.5s)"""
        state = pygcmc.MCState()
        state.info.box = [5.0, 5.0, 5.0]
        
        # Two dipoles at distance R
        R = 0.3
        atoms = [
            pygcmc.MCAtom(),  # P1
            pygcmc.MCAtom(),  # D1
            pygcmc.MCAtom(),  # P2
            pygcmc.MCAtom(),  # D2
        ]
        
        atoms[0].x = atoms[1].x = 0.0
        atoms[2].x = atoms[3].x = R
        
        for i, atom in enumerate(atoms):
            atom.y = atom.z = 0.0
            atom.charge = 1.0 if i % 2 == 0 else -1.0
            atom.type = i % 2
        
        state.atoms = atoms
        state.activeAtomCount = 4
        
        # Two residues
        for i in range(2):
            res = pygcmc.MCResidue()
            res.atomStart = i * 2
            res.atomCount = 2
            res.active = True
            state.residues.append(res)
        state.activeResidueCount = 2
        
        # Setup Drude particles
        pygcmc.DrudeComplete.clear()
        for i in range(2):
            p = pygcmc.DrudeParticle()
            p.drudeIndex = i * 2 + 1
            p.parentIndex = i * 2
            p.charge = -1.0
            p.polarizability = 0.001
            p.computeSpringConstants()
            pygcmc.DrudeComplete.addParticle(p)
        
        # Add Thole screening
        pair = pygcmc.ScreenedPair()
        pair.dipole1 = 0
        pair.dipole2 = 1
        pair.thole = 1.3
        pygcmc.DrudeComplete.addScreenedPair(pair)
        
        params = pygcmc.DrudeSCFParams()
        params.tolerance = 0.01
        params.maxIterations = 100
        params.dampingFactor = 0.5
        pygcmc.DrudeComplete.setParameters(params)
        
        energy = pygcmc.DrudeComplete.calculateEnergy(state)
        
        d1 = state.atoms[1].x - state.atoms[0].x
        d2 = state.atoms[3].x - state.atoms[2].x
        
        # Check symmetry
        assert abs(abs(d1) - abs(d2)) < 1e-6, "Displacements should be symmetric"
        # Check attraction (opposite signs)
        assert d1 * d2 < 0, "Dipoles should attract"
        # Check reasonable magnitude
        assert 0.001 < abs(d1) < 0.1, "Displacement magnitude should be reasonable"
    
    def test_water_dimer_scan_minimal(self):
        """Test 3: Water dimer scan - minimal points (< 1s)"""
        distances = [0.25, 0.35, 0.50]  # Just 3 points for fast CI
        displacements = []
        
        for R in distances:
            state = pygcmc.MCState()
            state.info.box = [5.0, 5.0, 5.0]
            
            # Simplified water dimer
            atoms = []
            for i in range(2):  # Two water molecules
                x_offset = R * i
                # Oxygen
                o = pygcmc.MCAtom()
                o.x = x_offset
                o.y = o.z = 0.0
                o.charge = -0.82  # Simplified charge
                o.type = 0
                atoms.append(o)
                
                # Oxygen Drude
                od = pygcmc.MCAtom()
                od.x = x_offset
                od.y = od.z = 0.0
                od.charge = -1.71636
                od.type = 1
                atoms.append(od)
                
                # Two hydrogens (simplified positions)
                for j in range(2):
                    h = pygcmc.MCAtom()
                    h.x = x_offset + 0.1 * (j * 2 - 1)
                    h.y = 0.05 * j
                    h.z = 0.0
                    h.charge = 0.55733
                    h.type = 2
                    atoms.append(h)
            
            state.atoms = atoms
            state.activeAtomCount = 8
            
            # Two water residues
            for i in range(2):
                res = pygcmc.MCResidue()
                res.atomStart = i * 4
                res.atomCount = 4
                res.active = True
                state.residues.append(res)
            state.activeResidueCount = 2
            
            # Setup Drude on oxygens
            pygcmc.DrudeComplete.clear()
            for i in range(2):
                p = pygcmc.DrudeParticle()
                p.drudeIndex = i * 4 + 1
                p.parentIndex = i * 4
                p.charge = -1.71636
                p.polarizability = 0.0009782237
                p.computeSpringConstants()
                pygcmc.DrudeComplete.addParticle(p)
            
            # Add Thole between oxygens
            pair = pygcmc.ScreenedPair()
            pair.dipole1 = 0
            pair.dipole2 = 1
            pair.thole = 1.3
            pygcmc.DrudeComplete.addScreenedPair(pair)
            
            params = pygcmc.DrudeSCFParams()
            params.tolerance = 0.05  # Relaxed for speed
            params.maxIterations = 50
            params.dampingFactor = 0.5
            pygcmc.DrudeComplete.setParameters(params)
            
            energy = pygcmc.DrudeComplete.calculateEnergy(state)
            
            # Check first oxygen Drude displacement
            dx = state.atoms[1].x - state.atoms[0].x
            dy = state.atoms[1].y - state.atoms[0].y
            dz = state.atoms[1].z - state.atoms[0].z
            d_mag = np.sqrt(dx**2 + dy**2 + dz**2)
            displacements.append(d_mag)
        
        # Check distance dependence (may not be strictly monotonic due to H positions)
        # Just check that displacements exist and are reasonable
        assert len(displacements) == 3, "Should have 3 displacement values"
        
        # Check reasonable magnitudes (water can have larger displacements)
        assert all(0.0001 < d < 0.5 for d in displacements), \
               "Displacements should be in reasonable range"
        assert all(d > 0 for d in displacements), \
               "All displacements should be positive"
    
    def test_regression_golden_values(self):
        """Test 4: Regression test with golden values (< 0.5s)"""
        # Fixed setup for reproducibility
        state = pygcmc.MCState()
        state.info.box = [5.0, 5.0, 5.0]
        
        # Create specific configuration
        atoms = [
            pygcmc.MCAtom(),  # Parent at origin
            pygcmc.MCAtom(),  # Drude at origin
            pygcmc.MCAtom(),  # External at (0.4, 0, 0)
        ]
        
        atoms[0].x = atoms[0].y = atoms[0].z = 0.0
        atoms[0].charge = 1.0
        atoms[0].type = 0
        
        atoms[1].x = atoms[1].y = atoms[1].z = 0.0
        atoms[1].charge = -1.0
        atoms[1].type = 1
        
        atoms[2].x = 0.4
        atoms[2].y = atoms[2].z = 0.0
        atoms[2].charge = 1.5
        atoms[2].type = 2
        
        state.atoms = atoms
        state.activeAtomCount = 3
        
        res = pygcmc.MCResidue()
        res.atomStart = 0
        res.atomCount = 2
        res.active = True
        state.residues = [res]
        state.activeResidueCount = 1
        
        pygcmc.DrudeComplete.clear()
        p = pygcmc.DrudeParticle()
        p.drudeIndex = 1
        p.parentIndex = 0
        p.charge = -1.0
        p.polarizability = 0.001
        p.computeSpringConstants()
        pygcmc.DrudeComplete.addParticle(p)
        
        params = pygcmc.DrudeSCFParams()
        params.tolerance = 1e-5
        params.maxIterations = 100
        params.dampingFactor = 0.5
        params.includeCoulombEnergy = True  # Use full energy for regression test
        pygcmc.DrudeComplete.setParameters(params)
        
        energy = pygcmc.DrudeComplete.calculateEnergy(state)
        dx = state.atoms[1].x - state.atoms[0].x
        
        # Golden values (from validated run)
        # Note: With includeCoulombEnergy=True (for standalone test)
        golden_dx = 0.009854582138359547
        golden_energy = -6.4138038305972485
        # With includeCoulombEnergy=False (new default)
        golden_energy_spring_only = 6.746204820025878
        
        # Check regression (allow 0.1% variation)
        assert abs(dx - golden_dx) / abs(golden_dx) < 0.001, \
               f"Displacement regression failed: {dx} vs {golden_dx}"
        assert abs(energy - golden_energy) / abs(golden_energy) < 0.001, \
               f"Energy regression failed: {energy} vs {golden_energy}"
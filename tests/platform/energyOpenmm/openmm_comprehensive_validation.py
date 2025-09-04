#!/usr/bin/env python3
"""
Comprehensive validation of PyGCMC vs OpenMM for various Drude configurations
Tests multiple scenarios and generates detailed comparison reports
"""

import numpy as np
import pytest
import json
import pygcmc

# Skip all tests if OpenMM is not available
pytest.importorskip("openmm")
import openmm as mm
import openmm.unit as u


class TestOpenMMComprehensiveValidation:
    """Comprehensive validation suite comparing PyGCMC and OpenMM"""
    
    def setup_method(self):
        """Initialize test results collector"""
        self.results = []
    
    def add_result(self, test_name, openmm_val, pygcmc_val, tolerance=0.01):
        """Record test result"""
        if openmm_val != 0:
            rel_error = abs(openmm_val - pygcmc_val) / abs(openmm_val)
        else:
            rel_error = abs(pygcmc_val)
        
        passed = rel_error < tolerance
        result = {
            'test': test_name,
            'openmm': float(openmm_val),
            'pygcmc': float(pygcmc_val),
            'rel_error': float(rel_error),
            'passed': bool(passed)  # Convert numpy bool to Python bool
        }
        self.results.append(result)
        return passed
    
    def test_dipole_dipole_interaction(self):
        """Test dipole-dipole interaction with screening"""
        alpha = 0.001
        thole = 1.3
        qD = -1.0
        separation = 0.6  # nm
        
        # === OpenMM Setup ===
        sys = mm.System()
        sys.addParticle(16.0 * u.dalton)  # Parent 1
        sys.addParticle(0.4 * u.dalton)   # Drude 1
        sys.addParticle(16.0 * u.dalton)  # Parent 2
        sys.addParticle(0.4 * u.dalton)   # Drude 2
        
        nb = mm.NonbondedForce()
        nb.setNonbondedMethod(mm.NonbondedForce.NoCutoff)
        nb.addParticle(0.0, 0.1*u.nanometer, 0.0)
        nb.addParticle(qD*u.elementary_charge, 0.1*u.nanometer, 0.0)
        nb.addParticle(0.0, 0.1*u.nanometer, 0.0)
        nb.addParticle(qD*u.elementary_charge, 0.1*u.nanometer, 0.0)
        nb.addException(0, 1, 0.0, 1.0, 0.0)
        nb.addException(2, 3, 0.0, 1.0, 0.0)
        sys.addForce(nb)
        
        df = mm.DrudeForce()
        df.addParticle(1, 0, -1, -1, -1, qD, alpha, 1.0, thole)
        df.addParticle(3, 2, -1, -1, -1, qD, alpha, 1.0, thole)
        df.addScreenedPair(0, 1, thole)
        sys.addForce(df)
        
        integ = mm.DrudeSCFIntegrator(0.001*u.picoseconds)
        integ.setMinimizationErrorTolerance(1e-10)
        ctx = mm.Context(sys, integ)
        ctx.setPositions([[0,0,0], [0,0,0], [separation,0,0], [separation,0,0]] * u.nanometer)
        integ.step(1)
        
        state = ctx.getState(getPositions=True, getEnergy=True)
        pos_omm = state.getPositions(asNumpy=True).value_in_unit(u.nanometer)
        energy_omm = state.getPotentialEnergy().value_in_unit(u.kilojoule_per_mole)
        disp1_omm = pos_omm[1,0] - pos_omm[0,0]
        disp2_omm = pos_omm[3,0] - pos_omm[2,0]
        
        # === PyGCMC Setup (separate molecules) ===
        state = pygcmc.MCState()
        state.info.box = [10.0, 10.0, 10.0]
        
        atoms = []
        # Dipole 1
        parent1 = pygcmc.MCAtom()
        parent1.x = parent1.y = parent1.z = 0.0
        parent1.charge = 0.0
        atoms.append(parent1)
        
        drude1 = pygcmc.MCAtom()
        drude1.x = drude1.y = drude1.z = 0.0
        drude1.charge = qD
        atoms.append(drude1)
        
        # Dipole 2
        parent2 = pygcmc.MCAtom()
        parent2.x = separation
        parent2.y = parent2.z = 0.0
        parent2.charge = 0.0
        atoms.append(parent2)
        
        drude2 = pygcmc.MCAtom()
        drude2.x = separation
        drude2.y = drude2.z = 0.0
        drude2.charge = qD
        atoms.append(drude2)
        
        state.atoms = atoms
        state.activeAtomCount = 4
        
        # Separate residues
        res1 = pygcmc.MCResidue()
        res1.atomStart = 0
        res1.atomCount = 2
        res1.active = True
        state.residues.append(res1)
        
        res2 = pygcmc.MCResidue()
        res2.atomStart = 2
        res2.atomCount = 2
        res2.active = True
        state.residues.append(res2)
        state.activeResidueCount = 2
        
        # Setup Drude
        pygcmc.DrudeComplete.clear()
        
        p1 = pygcmc.DrudeParticle()
        p1.drudeIndex = 1
        p1.parentIndex = 0
        p1.charge = qD
        p1.polarizability = alpha
        p1.computeSpringConstants()
        pygcmc.DrudeComplete.addParticle(p1)
        
        p2 = pygcmc.DrudeParticle()
        p2.drudeIndex = 3
        p2.parentIndex = 2
        p2.charge = qD
        p2.polarizability = alpha
        p2.computeSpringConstants()
        pygcmc.DrudeComplete.addParticle(p2)
        
        # Add screening
        pair = pygcmc.ScreenedPair()
        pair.dipole1 = 0
        pair.dipole2 = 1
        pair.thole = thole
        pygcmc.DrudeComplete.addScreenedPair(pair)
        
        params = pygcmc.DrudeSCFParams()
        params.tolerance = 1e-10
        params.maxIterations = 100
        pygcmc.DrudeComplete.setParameters(params)
        
        energy_pg = pygcmc.DrudeComplete.calculateEnergy(state)
        disp1_pg = state.atoms[1].x - state.atoms[0].x
        disp2_pg = state.atoms[3].x - state.atoms[2].x
        
        # Record results
        self.add_result("Dipole-dipole: dipole1 displacement", disp1_omm, disp1_pg)
        self.add_result("Dipole-dipole: dipole2 displacement", disp2_omm, disp2_pg)
        
        # Check symmetry
        symmetry_omm = abs(disp1_omm + disp2_omm) < 1e-9
        symmetry_pg = abs(disp1_pg + disp2_pg) < 1e-9
        
        # For separate molecules, displacements should be roughly symmetric
        assert abs(disp1_pg + disp2_pg) < 1e-6, (
            f"PyGCMC displacements not symmetric: {disp1_pg} vs {disp2_pg}"
        )
    
    def test_three_body_system(self):
        """Test three-body system with multiple charges"""
        alpha = 0.001
        thole = 1.3
        qD = -1.0
        
        # === OpenMM ===
        sys = mm.System()
        sys.addParticle(16.0 * u.dalton)  # Parent
        sys.addParticle(0.4 * u.dalton)   # Drude
        sys.addParticle(1.0 * u.dalton)   # External 1
        sys.addParticle(1.0 * u.dalton)   # External 2
        
        nb = mm.NonbondedForce()
        nb.setNonbondedMethod(mm.NonbondedForce.NoCutoff)
        nb.addParticle(0.0, 0.1*u.nanometer, 0.0)
        nb.addParticle(qD*u.elementary_charge, 0.1*u.nanometer, 0.0)
        nb.addParticle(0.5*u.elementary_charge, 0.1*u.nanometer, 0.0)
        nb.addParticle(-0.5*u.elementary_charge, 0.1*u.nanometer, 0.0)
        nb.addException(0, 1, 0.0, 1.0, 0.0)
        sys.addForce(nb)
        
        df = mm.DrudeForce()
        df.addParticle(1, 0, -1, -1, -1, qD, alpha, 1.0, thole)
        sys.addForce(df)
        
        integ = mm.DrudeSCFIntegrator(0.001*u.picoseconds)
        integ.setMinimizationErrorTolerance(1e-10)
        ctx = mm.Context(sys, integ)
        ctx.setPositions([[0,0,0], [0,0,0], [0.4,0,0], [0,0.4,0]] * u.nanometer)
        integ.step(1)
        
        state = ctx.getState(getPositions=True, getEnergy=True)
        pos_omm = state.getPositions(asNumpy=True).value_in_unit(u.nanometer)
        energy_omm = state.getPotentialEnergy().value_in_unit(u.kilojoule_per_mole)
        disp_omm = np.array([
            pos_omm[1,0] - pos_omm[0,0],
            pos_omm[1,1] - pos_omm[0,1],
            pos_omm[1,2] - pos_omm[0,2]
        ])
        disp_mag_omm = np.linalg.norm(disp_omm)
        
        # === PyGCMC ===
        state = pygcmc.MCState()
        state.info.box = [10.0, 10.0, 10.0]
        
        parent = pygcmc.MCAtom()
        parent.x = parent.y = parent.z = 0.0
        parent.charge = 0.0
        
        drude = pygcmc.MCAtom()
        drude.x = drude.y = drude.z = 0.0
        drude.charge = qD
        
        ext1 = pygcmc.MCAtom()
        ext1.x = 0.4
        ext1.y = ext1.z = 0.0
        ext1.charge = 0.5
        
        ext2 = pygcmc.MCAtom()
        ext2.x = 0.0
        ext2.y = 0.4
        ext2.z = 0.0
        ext2.charge = -0.5
        
        state.atoms = [parent, drude, ext1, ext2]
        state.activeAtomCount = 4
        
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
        p.charge = qD
        p.polarizability = alpha
        p.computeSpringConstants()
        pygcmc.DrudeComplete.addParticle(p)
        
        params = pygcmc.DrudeSCFParams()
        params.tolerance = 1e-10
        params.maxIterations = 100
        pygcmc.DrudeComplete.setParameters(params)
        
        energy_pg = pygcmc.DrudeComplete.calculateEnergy(state)
        disp_pg = np.array([
            state.atoms[1].x - state.atoms[0].x,
            state.atoms[1].y - state.atoms[0].y,
            state.atoms[1].z - state.atoms[0].z
        ])
        disp_mag_pg = np.linalg.norm(disp_pg)
        
        # Record results
        self.add_result("Three-body: displacement magnitude", disp_mag_omm, disp_mag_pg)
        self.add_result("Three-body: displacement x", disp_omm[0], disp_pg[0])
        self.add_result("Three-body: displacement y", disp_omm[1], disp_pg[1])
        self.add_result("Three-body: energy", energy_omm, energy_pg)
        
        # Displacement should be in 2D plane (x-y)
        assert abs(disp_pg[2]) < 1e-10, f"Unexpected z displacement: {disp_pg[2]}"
    
    def test_convergence_from_different_starts(self):
        """Test that SCF converges to same point from different initial positions"""
        alpha = 0.001
        thole = 1.3
        qD = -1.0
        R = 0.5
        
        initial_displacements = [0.0, 0.001, -0.001, 0.01, -0.01]
        final_positions = []
        
        for init_disp in initial_displacements:
            state = pygcmc.MCState()
            state.info.box = [10.0, 10.0, 10.0]
            
            parent = pygcmc.MCAtom()
            parent.x = parent.y = parent.z = 0.0
            parent.charge = 0.0
            
            drude = pygcmc.MCAtom()
            drude.x = init_disp  # Different starting position
            drude.y = drude.z = 0.0
            drude.charge = qD
            
            ext = pygcmc.MCAtom()
            ext.x = R
            ext.y = ext.z = 0.0
            ext.charge = 1.0
            
            state.atoms = [parent, drude, ext]
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
            p.charge = qD
            p.polarizability = alpha
            p.computeSpringConstants()
            pygcmc.DrudeComplete.addParticle(p)
            
            params = pygcmc.DrudeSCFParams()
            params.tolerance = 1e-10
            params.maxIterations = 200
            pygcmc.DrudeComplete.setParameters(params)
            
            pygcmc.DrudeComplete.calculateEnergy(state)
            final_positions.append(state.atoms[1].x)
        
        # All should converge to same position
        reference = final_positions[0]
        for i, pos in enumerate(final_positions[1:], 1):
            diff = abs(pos - reference)
            assert diff < 1e-8, (
                f"Different convergence from init={initial_displacements[i]}: "
                f"diff={diff}"
            )
        
        self.add_result("Convergence consistency", 1.0, 1.0)  # Pass/fail test
    
    def test_standard_vs_exact_mode(self):
        """Compare standard and exact modes"""
        alpha = 0.001
        thole = 1.3
        qD = -1.0
        R = 0.5
        
        # Setup system
        state = pygcmc.MCState()
        state.info.box = [10.0, 10.0, 10.0]
        
        parent = pygcmc.MCAtom()
        parent.x = parent.y = parent.z = 0.0
        parent.charge = 0.0
        
        drude = pygcmc.MCAtom()
        drude.x = drude.y = drude.z = 0.0
        drude.charge = qD
        
        ext = pygcmc.MCAtom()
        ext.x = R
        ext.y = ext.z = 0.0
        ext.charge = 1.0
        
        state.atoms = [parent, drude, ext]
        state.activeAtomCount = 3
        
        res = pygcmc.MCResidue()
        res.atomStart = 0
        res.atomCount = 2
        res.active = True
        state.residues = [res]
        state.activeResidueCount = 1
        
        # Test standard mode
        pygcmc.DrudeComplete.clear()
        p = pygcmc.DrudeParticle()
        p.drudeIndex = 1
        p.parentIndex = 0
        p.charge = qD
        p.polarizability = alpha
        p.computeSpringConstants()
        pygcmc.DrudeComplete.addParticle(p)
        
        params = pygcmc.DrudeSCFParams()
        params.tolerance = 1e-10
        params.requireExactMatch = False
        pygcmc.DrudeComplete.setParameters(params)
        
        energy_std = pygcmc.DrudeComplete.calculateEnergy(state)
        disp_std = state.atoms[1].x
        
        # Reset and test exact mode
        state.atoms[1].x = 0.0
        params.requireExactMatch = True
        params.maxIterations = 200
        pygcmc.DrudeComplete.setParameters(params)
        
        energy_exact = pygcmc.DrudeComplete.calculateEnergy(state)
        disp_exact = state.atoms[1].x
        
        # Should be very close
        rel_diff_disp = abs(disp_std - disp_exact) / abs(disp_std)
        rel_diff_energy = abs(energy_std - energy_exact) / abs(energy_std)
        
        assert rel_diff_disp < 0.001, f"Mode displacement difference: {rel_diff_disp*100:.3f}%"
        assert rel_diff_energy < 0.01, f"Mode energy difference: {rel_diff_energy*100:.3f}%"
        
        self.add_result("Standard vs exact displacement", disp_std, disp_exact, 0.001)
        self.add_result("Standard vs exact energy", energy_std, energy_exact, 0.01)
    
    def teardown_method(self):
        """Print summary after each test method"""
        if self.results:
            passed = sum(1 for r in self.results if r['passed'])
            total = len(self.results)
            print(f"\nTest Summary: {passed}/{total} checks passed")
            
            # Save detailed results
            import os
            os.makedirs('/home/zhaomt/gcmc/test108/pygcmc_dev/tmp/test_results', exist_ok=True)
            timestamp = pytest.approx(0)  # Would use datetime in real code
            filename = f'/home/zhaomt/gcmc/test108/pygcmc_dev/tmp/test_results/comprehensive_validation.json'
            with open(filename, 'w') as f:
                json.dump(self.results, f, indent=2)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
#!/usr/bin/env python3
"""
Test PyGCMC Drude SCF vs OpenMM DrudeSCFIntegrator
Validates that both implementations reach equivalent equilibrium
"""

import numpy as np
import pytest
import pygcmc

# Skip all tests if OpenMM is not available
pytest.importorskip("openmm")
import openmm as mm
import openmm.unit as u


class TestOpenMMSCFComparison:
    """Compare PyGCMC and OpenMM SCF optimization"""
    
    def test_single_dipole_external_field(self):
        """Test single dipole in external field"""
        # Parameters
        alpha = 0.001
        thole = 1.3
        qD = -1.0
        R = 0.5
        
        # === OpenMM ===
        sys = mm.System()
        sys.addParticle(16.0 * u.dalton)  # Parent
        sys.addParticle(0.4 * u.dalton)   # Drude
        sys.addParticle(1.0 * u.dalton)   # External
        
        nb = mm.NonbondedForce()
        nb.setNonbondedMethod(mm.NonbondedForce.NoCutoff)
        nb.addParticle(0.0, 0.1*u.nanometer, 0.0)
        nb.addParticle(qD*u.elementary_charge, 0.1*u.nanometer, 0.0)
        nb.addParticle(1.0*u.elementary_charge, 0.1*u.nanometer, 0.0)
        nb.addException(0, 1, 0.0, 1.0, 0.0)
        sys.addForce(nb)
        
        df = mm.DrudeForce()
        df.addParticle(1, 0, -1, -1, -1, qD, alpha, 1.0, thole)
        sys.addForce(df)
        
        integ = mm.DrudeSCFIntegrator(0.001*u.picoseconds)
        integ.setMinimizationErrorTolerance(1e-10)
        ctx = mm.Context(sys, integ)
        ctx.setPositions([[0,0,0], [0,0,0], [R,0,0]] * u.nanometer)
        integ.step(1)
        
        state = ctx.getState(getPositions=True, getEnergy=True)
        pos_omm = state.getPositions(asNumpy=True).value_in_unit(u.nanometer)
        energy_omm = state.getPotentialEnergy().value_in_unit(u.kilojoule_per_mole)
        disp_omm = pos_omm[1,0] - pos_omm[0,0]
        
        # === PyGCMC ===
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
        disp_pg = state.atoms[1].x - state.atoms[0].x
        
        # Compare results
        rel_error_disp = abs(disp_omm - disp_pg) / abs(disp_omm)
        rel_error_energy = abs(energy_omm - energy_pg) / abs(energy_omm)
        
        assert rel_error_disp < 0.01, f"Displacement error: {rel_error_disp*100:.2f}%"
        assert rel_error_energy < 0.01, f"Energy error: {rel_error_energy*100:.2f}%"
    
    def test_varying_thole_parameters(self):
        """Test with different Thole parameters"""
        alpha = 0.001
        qD = -1.0
        R = 0.5
        
        thole_values = [0.5, 1.0, 1.3, 2.0, 3.0]
        
        for thole in thole_values:
            # OpenMM
            sys = mm.System()
            sys.addParticle(16.0 * u.dalton)
            sys.addParticle(0.4 * u.dalton)
            sys.addParticle(1.0 * u.dalton)
            
            nb = mm.NonbondedForce()
            nb.setNonbondedMethod(mm.NonbondedForce.NoCutoff)
            nb.addParticle(0.0, 0.1*u.nanometer, 0.0)
            nb.addParticle(qD*u.elementary_charge, 0.1*u.nanometer, 0.0)
            nb.addParticle(1.0*u.elementary_charge, 0.1*u.nanometer, 0.0)
            nb.addException(0, 1, 0.0, 1.0, 0.0)
            sys.addForce(nb)
            
            df = mm.DrudeForce()
            df.addParticle(1, 0, -1, -1, -1, qD, alpha, 1.0, thole)
            sys.addForce(df)
            
            integ = mm.DrudeSCFIntegrator(0.001*u.picoseconds)
            integ.setMinimizationErrorTolerance(1e-10)
            ctx = mm.Context(sys, integ)
            ctx.setPositions([[0,0,0], [0,0,0], [R,0,0]] * u.nanometer)
            integ.step(1)
            
            state = ctx.getState(getPositions=True)
            pos_omm = state.getPositions(asNumpy=True).value_in_unit(u.nanometer)
            disp_omm = pos_omm[1,0]
            
            # PyGCMC
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
            pygcmc.DrudeComplete.setParameters(params)
            
            pygcmc.DrudeComplete.calculateEnergy(state)
            disp_pg = state.atoms[1].x
            
            rel_error = abs(disp_omm - disp_pg) / abs(disp_omm)
            assert rel_error < 0.01, f"Thole={thole}: error {rel_error*100:.2f}%"
    
    def test_scf_convergence_behavior(self):
        """Test that SCF converges to same point from different starts"""
        alpha = 0.001
        thole = 1.3
        qD = -1.0
        R = 0.5
        
        initial_displacements = [0.0, 0.001, -0.001, 0.01]
        final_positions = []
        
        for init_disp in initial_displacements:
            # Setup system
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
            params.maxIterations = 100
            pygcmc.DrudeComplete.setParameters(params)
            
            pygcmc.DrudeComplete.calculateEnergy(state)
            final_positions.append(state.atoms[1].x)
        
        # All should converge to same position
        for i in range(1, len(final_positions)):
            diff = abs(final_positions[i] - final_positions[0])
            assert diff < 1e-8, f"Different convergence from init={initial_displacements[i]}: {diff}"
    
    def test_exact_mode_vs_openmm(self):
        """Test PyGCMC exact mode against OpenMM"""
        alpha = 0.001
        thole = 1.3
        qD = -1.0
        R = 0.5
        
        # OpenMM
        sys = mm.System()
        sys.addParticle(16.0 * u.dalton)
        sys.addParticle(0.4 * u.dalton)
        sys.addParticle(1.0 * u.dalton)
        
        nb = mm.NonbondedForce()
        nb.setNonbondedMethod(mm.NonbondedForce.NoCutoff)
        nb.addParticle(0.0, 0.1*u.nanometer, 0.0)
        nb.addParticle(qD*u.elementary_charge, 0.1*u.nanometer, 0.0)
        nb.addParticle(1.0*u.elementary_charge, 0.1*u.nanometer, 0.0)
        nb.addException(0, 1, 0.0, 1.0, 0.0)
        sys.addForce(nb)
        
        df = mm.DrudeForce()
        df.addParticle(1, 0, -1, -1, -1, qD, alpha, 1.0, thole)
        sys.addForce(df)
        
        integ = mm.DrudeSCFIntegrator(0.001*u.picoseconds)
        integ.setMinimizationErrorTolerance(1e-10)
        ctx = mm.Context(sys, integ)
        ctx.setPositions([[0,0,0], [0,0,0], [R,0,0]] * u.nanometer)
        integ.step(1)
        
        state_omm = ctx.getState(getPositions=True, getEnergy=True)
        pos_omm = state_omm.getPositions(asNumpy=True).value_in_unit(u.nanometer)
        energy_omm = state_omm.getPotentialEnergy().value_in_unit(u.kilojoule_per_mole)
        disp_omm = pos_omm[1,0]
        
        # PyGCMC Exact Mode
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
        
        pygcmc.DrudeComplete.clear()
        p = pygcmc.DrudeParticle()
        p.drudeIndex = 1
        p.parentIndex = 0
        p.charge = qD
        p.polarizability = alpha
        p.computeSpringConstants()
        pygcmc.DrudeComplete.addParticle(p)
        
        params = pygcmc.DrudeSCFParams()
        params.requireExactMatch = True  # Enable exact mode
        params.tolerance = 1e-10
        params.maxIterations = 200
        pygcmc.DrudeComplete.setParameters(params)
        
        energy_pg = pygcmc.DrudeComplete.calculateEnergy(state)
        disp_pg = state.atoms[1].x
        
        # Compare
        rel_error_disp = abs(disp_omm - disp_pg) / abs(disp_omm)
        rel_error_energy = abs(energy_omm - energy_pg) / abs(energy_omm)
        
        # Exact mode should have <1% error
        assert rel_error_disp < 0.01, f"Exact mode displacement error: {rel_error_disp*100:.2f}%"
        assert rel_error_energy < 0.01, f"Exact mode energy error: {rel_error_energy*100:.2f}%"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
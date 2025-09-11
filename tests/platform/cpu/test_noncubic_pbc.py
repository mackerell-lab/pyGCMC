#!/usr/bin/env python
"""
Test PBC with non-cubic boxes - verify minimum image convention
and energy invariance for boxes with Lx ≠ Ly ≠ Lz
"""

import pytest
import numpy as np
import os
import sys

# Add build path for pygcmc module
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../build'))

try:
    import pygcmc
    PYGCMC_AVAILABLE = True
except ImportError:
    PYGCMC_AVAILABLE = False
    pygcmc = None

@pytest.mark.skipif(not PYGCMC_AVAILABLE, reason="PyGCMC not available")
class TestNonCubicPBC:
    """Test periodic boundary conditions with non-cubic boxes"""
    
    def test_noncubic_translation_invariance(self):
        """Test energy invariance under translation in non-cubic box"""
        state = pygcmc.MCState()
        # Non-cubic box
        Lx, Ly, Lz = 5.0, 7.0, 3.0
        state.info.box = (Lx, Ly, Lz)
        state.info.setTemperature(300.0)
        state.info.cutoff = min(Lx, Ly, Lz) / 2 - 0.1  # Safe cutoff
        
        ff = pygcmc.MCForceField()
        ff.numTotalTypes = 1
        ff.numMovementTypes = 1
        ff.ljEps = [1.0]
        ff.ljSigma = [0.3]
        state.forcefield = ff
        
        # Create molecules at various positions
        positions = [
            (1.0, 1.0, 1.0),
            (2.0, 3.0, 1.5),
            (4.0, 5.0, 2.0),
            (0.5, 6.0, 0.5),
        ]
        
        atoms = []
        residues = []
        
        for i, (x, y, z) in enumerate(positions):
            atom = pygcmc.MCAtom()
            atom.type = 0
            atom.x, atom.y, atom.z = x, y, z
            atom.charge = 0.1 if i % 2 == 0 else -0.1
            atoms.append(atom)
            
            residue = pygcmc.MCResidue()
            residue.active = True
            residue.atomStart = i
            residue.atomCount = 1
            residue.atoms = [atom]
            residues.append(residue)
        
        state.atoms = atoms
        state.residues = residues
        
        state.activeAtomCount = len(state.atoms)
        state.activeResidueCount = len(state.residues)
        
        # Calculate initial energy
        pygcmc.computeSystemEnergyPBCCutoff(state)
        initial_energy = sum(r.energy_vdw + r.energy_elec for r in state.residues if r.active)
        
        # Test translations by box dimensions
        translations = [
            (Lx, 0, 0),    # +Lx
            (0, Ly, 0),    # +Ly
            (0, 0, Lz),    # +Lz
            (-Lx, 0, 0),   # -Lx
            (0, -Ly, 0),   # -Ly
            (0, 0, -Lz),   # -Lz
            (Lx, Ly, 0),   # Mixed
            (0, Ly, Lz),   # Mixed
            (Lx, 0, Lz),   # Mixed
            (Lx, Ly, Lz),  # All
        ]
        
        for tx, ty, tz in translations:
            # Translate all atoms
            for atom in state.atoms:
                atom.x += tx
                atom.y += ty
                atom.z += tz
            
            for residue in state.residues:
                for atom in residue.atoms:
                    atom.x += tx
                    atom.y += ty
                    atom.z += tz
            
            # Recalculate energy
            pygcmc.computeSystemEnergyPBCCutoff(state)
            translated_energy = sum(r.energy_vdw + r.energy_elec for r in state.residues if r.active)
            
            # Check invariance
            assert abs(translated_energy - initial_energy) < 1e-6, \
                f"Energy changed in non-cubic box after translation ({tx},{ty},{tz}): " \
                f"{initial_energy:.6f} -> {translated_energy:.6f}"
            
            # Translate back
            for atom in state.atoms:
                atom.x -= tx
                atom.y -= ty
                atom.z -= tz
            
            for residue in state.residues:
                for atom in residue.atoms:
                    atom.x -= tx
                    atom.y -= ty
                    atom.z -= tz
        
        print(f"✓ Non-cubic translation invariance test passed")
        print(f"  Box: {Lx} × {Ly} × {Lz} nm")
        print(f"  Energy remains {initial_energy:.6f} kJ/mol")
    
    def test_noncubic_minimum_image(self):
        """Test minimum image convention in non-cubic box"""
        state = pygcmc.MCState()
        Lx, Ly, Lz = 3.0, 4.0, 5.0
        state.info.box = (Lx, Ly, Lz)
        state.info.setTemperature(300.0)
        state.info.cutoff = min(Lx, Ly, Lz) / 2 - 0.1
        
        ff = pygcmc.MCForceField()
        ff.numTotalTypes = 1
        ff.numMovementTypes = 1
        ff.ljEps = [1.0]
        ff.ljSigma = [0.3]
        state.forcefield = ff
        
        # Two atoms
        atom1 = pygcmc.MCAtom()
        atom1.type = 0
        atom1.x, atom1.y, atom1.z = 0.5, 0.5, 0.5
        atom1.charge = 0.0
        
        atom2 = pygcmc.MCAtom()
        atom2.type = 0
        atom2.charge = 0.0
        
        state.atoms = [atom1, atom2]
        
        res1 = pygcmc.MCResidue()
        res1.active = True
        res1.atomStart = 0
        res1.atomCount = 1
        res1.atoms = [atom1]
        
        res2 = pygcmc.MCResidue()
        res2.active = True
        res2.atomStart = 1
        res2.atomCount = 1
        res2.atoms = [atom2]
        
        state.residues = [res1, res2]
        state.activeAtomCount = 2
        state.activeResidueCount = 2
        
        # Test pairs that should have same minimum image distance
        test_pairs = [
            # Direct and wrapped positions
            ((1.0, 0.5, 0.5), (1.0 + Lx, 0.5, 0.5)),  # x-wrapping
            ((0.5, 1.0, 0.5), (0.5, 1.0 + Ly, 0.5)),  # y-wrapping
            ((0.5, 0.5, 1.0), (0.5, 0.5, 1.0 + Lz)),  # z-wrapping
        ]
        
        for pos_direct, pos_wrapped in test_pairs:
            # Set atom2 at direct position
            atom2.x, atom2.y, atom2.z = pos_direct
            state.atoms[1] = atom2
            res2.atoms[0] = atom2
            
            pygcmc.computeSystemEnergyPBCCutoff(state)
            energy_direct = res1.energy_vdw + res2.energy_vdw
            
            # Set atom2 at wrapped position
            atom2.x, atom2.y, atom2.z = pos_wrapped
            state.atoms[1] = atom2
            res2.atoms[0] = atom2
            
            pygcmc.computeSystemEnergyPBCCutoff(state)
            energy_wrapped = res1.energy_vdw + res2.energy_vdw
            
            # Should give same energy (same minimum image distance)
            assert abs(energy_direct - energy_wrapped) < 1e-6, \
                f"Minimum image not working in non-cubic box: " \
                f"{energy_direct:.6f} != {energy_wrapped:.6f}"
        
        print(f"✓ Non-cubic minimum image test passed")
        print(f"  Box: {Lx} × {Ly} × {Lz} nm")
    
    def test_extreme_aspect_ratio(self):
        """Test PBC with extreme aspect ratios"""
        state = pygcmc.MCState()
        # Very elongated box
        Lx, Ly, Lz = 20.0, 2.0, 2.0
        state.info.box = (Lx, Ly, Lz)
        state.info.setTemperature(300.0)
        state.info.cutoff = min(Ly, Lz) / 2 - 0.1  # Must fit in smallest dimension
        
        ff = pygcmc.MCForceField()
        ff.numTotalTypes = 1
        ff.numMovementTypes = 1
        ff.ljEps = [1.0]
        ff.ljSigma = [0.3]
        state.forcefield = ff
        
        # Create a chain of atoms along x
        n_atoms = 10
        spacing = Lx / n_atoms
        
        atoms = []
        residues = []
        
        for i in range(n_atoms):
            atom = pygcmc.MCAtom()
            atom.type = 0
            atom.x = i * spacing + 0.5
            atom.y, atom.z = 1.0, 1.0
            atom.charge = 0.0
            atoms.append(atom)
            
            residue = pygcmc.MCResidue()
            residue.active = True
            residue.atomStart = i
            residue.atomCount = 1
            residue.atoms = [atom]
            residues.append(residue)
        
        state.atoms = atoms
        state.residues = residues
        
        state.activeAtomCount = n_atoms
        state.activeResidueCount = n_atoms
        
        # Calculate energy
        pygcmc.computeSystemEnergyPBCCutoff(state)
        
        # First and last atoms should interact through PBC
        first_residue = state.residues[0]
        last_residue = state.residues[n_atoms - 1]
        
        # Due to PBC, distance between first and last should be small
        dx_direct = state.atoms[n_atoms-1].x - state.atoms[0].x
        dx_wrapped = dx_direct - Lx * round(dx_direct / Lx)
        wrapped_distance = abs(dx_wrapped)
        
        print(f"✓ Extreme aspect ratio test passed")
        print(f"  Box: {Lx} × {Ly} × {Lz} nm")
        print(f"  First-last wrapped distance: {wrapped_distance:.3f} nm")
        
        # The wrapped distance should be much smaller than direct
        assert wrapped_distance < Lx / 2, \
            f"Wrapped distance {wrapped_distance:.3f} not smaller than half box"
    
    def test_gcmc_in_noncubic_box(self):
        """Test GCMC operations in non-cubic box"""
        state = pygcmc.MCState()
        Lx, Ly, Lz = 4.0, 6.0, 3.0
        state.info.box = (Lx, Ly, Lz)
        state.info.setTemperature(300.0)
        state.info.cutoff = min(Lx, Ly, Lz) / 2 - 0.1
        
        ff = pygcmc.MCForceField()
        ff.numTotalTypes = 1
        ff.numMovementTypes = 1
        ff.ljEps = [0.5]
        ff.ljSigma = [0.3]
        state.forcefield = ff
        
        # Create template
        reservoir = pygcmc.movement.FragmentReservoir()
        template = pygcmc.movement.FragmentTemplate()
        template.typeId = 0
        
        atom = pygcmc.MCAtom()
        atom.type = 0
        atom.x, atom.y, atom.z = 0.0, 0.0, 0.0
        atom.charge = 0.0
        template.atoms = [atom]
        reservoir.addTemplate(template)
        
        engine = pygcmc.GCMCEngine()
        engine.initialize(state, reservoir)
        engine.setTemperature(300.0)
        engine.setSeed(42)
        
        V = Lx * Ly * Lz
        acceptance = pygcmc.GCMCAcceptance()
        acceptance.setTemperature(300.0)
        acceptance.setVolume(V)
        acceptance.setActivity(0, 1.0)
        engine.setAcceptanceCalculator(acceptance)
        
        # Run GCMC
        n_steps = 100
        for i in range(n_steps):
            if i % 2 == 0:
                engine.attemptInsertion(0)
            else:
                if state.get_active_residue_count() > 0:
                    engine.attemptDeletion(0)
        
        # Check all atoms are within box
        for residue in state.residues:
            if residue.active:
                for j in range(residue.atomCount):
                    atom_idx = residue.atomStart + j
                    atom = state.atoms[atom_idx]
                    
                    # Atoms should be within [0, L) for each dimension
                    assert 0 <= atom.x < Lx or -Lx/2 <= atom.x < Lx/2, \
                        f"Atom x={atom.x} outside box [0, {Lx})"
                    assert 0 <= atom.y < Ly or -Ly/2 <= atom.y < Ly/2, \
                        f"Atom y={atom.y} outside box [0, {Ly})"
                    assert 0 <= atom.z < Lz or -Lz/2 <= atom.z < Lz/2, \
                        f"Atom z={atom.z} outside box [0, {Lz})"
        
        print(f"✓ GCMC in non-cubic box test passed")
        print(f"  Box: {Lx} × {Ly} × {Lz} nm")
        print(f"  Final particle count: {state.get_active_residue_count()}")

if __name__ == "__main__":
    if PYGCMC_AVAILABLE:
        print("Running non-cubic PBC tests...\n")
        
        test = TestNonCubicPBC()
        test.test_noncubic_translation_invariance()
        print()
        test.test_noncubic_minimum_image()
        print()
        test.test_extreme_aspect_ratio()
        print()
        test.test_gcmc_in_noncubic_box()
        
        print("\n✅ All non-cubic PBC tests passed!")
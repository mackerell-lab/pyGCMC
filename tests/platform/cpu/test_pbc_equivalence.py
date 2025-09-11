#!/usr/bin/env python
"""
Test PBC equivalence - energy invariance under translation and mirroring
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
class TestPBCEquivalence:
    """Test periodic boundary conditions equivalence"""
    
    def test_translation_invariance(self):
        """Test that energy is invariant under whole-system translation"""
        state = pygcmc.MCState()
        box_size = 5.0
        state.info.box = (box_size, box_size, box_size)
        state.info.setTemperature(300.0)
        state.info.cutoff = 2.0
        
        # Set up force field
        ff = pygcmc.MCForceField()
        ff.numTotalTypes = 1
        ff.numMovementTypes = 1
        ff.ljEps = [1.0]
        ff.ljSigma = [0.3]
        state.forcefield = ff
        
        # Create some molecules manually
        state.atoms = []
        state.residues = []
        state.activeAtomCount = 0
        state.activeResidueCount = 0
        
        # Add a few atoms at specific positions
        positions = [
            (1.0, 1.0, 1.0),
            (2.0, 1.5, 1.0),
            (1.5, 2.0, 1.5),
            (3.0, 3.0, 3.0),
            (4.0, 4.0, 4.0)
        ]
        
        for i, (x, y, z) in enumerate(positions):
            atom = pygcmc.MCAtom()
            atom.type = 0
            atom.x, atom.y, atom.z = x, y, z
            atom.charge = 0.0
            state.atoms.append(atom)
            
            # Create residue for each atom
            residue = pygcmc.MCResidue()
            residue.active = True
            residue.atomStart = i
            residue.atomCount = 1
            residue.atoms = [atom]
            state.residues.append(residue)
        
        state.activeAtomCount = len(state.atoms)
        state.activeResidueCount = len(state.residues)
        
        # Calculate initial energy with PBC
        pygcmc.computeSystemEnergyPBCCutoff(state)
        initial_energy = sum(r.energy_vdw + r.energy_elec for r in state.residues if r.active)
        
        # Test various translations
        translations = [
            (box_size, 0, 0),      # +L in x
            (0, box_size, 0),      # +L in y
            (0, 0, box_size),      # +L in z
            (-box_size, 0, 0),     # -L in x
            (0, -box_size, 0),     # -L in y
            (0, 0, -box_size),     # -L in z
            (box_size, box_size, 0),     # +L in x,y
            (box_size, box_size, box_size)  # +L in all
        ]
        
        for tx, ty, tz in translations:
            # Translate all atoms
            for atom in state.atoms:
                atom.x += tx
                atom.y += ty
                atom.z += tz
                
            # Update residue atoms too
            for residue in state.residues:
                for atom in residue.atoms:
                    atom.x += tx
                    atom.y += ty
                    atom.z += tz
                    
            # Recalculate energy
            pygcmc.computeSystemEnergyPBCCutoff(state)
            translated_energy = sum(r.energy_vdw + r.energy_elec for r in state.residues if r.active)
            
            # Check energy invariance
            assert abs(translated_energy - initial_energy) < 1e-6, \
                f"Energy changed after translation ({tx},{ty},{tz}): {initial_energy:.6f} -> {translated_energy:.6f}"
            
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
                
        print(f"✓ Translation invariance test passed")
        print(f"  Energy remains {initial_energy:.6f} kJ/mol under all translations")
    
    def test_mirror_symmetry(self):
        """Test that energy is invariant under mirroring"""
        state = pygcmc.MCState()
        box_size = 5.0
        state.info.box = (box_size, box_size, box_size)
        state.info.setTemperature(300.0)
        state.info.cutoff = 2.0
        
        ff = pygcmc.MCForceField()
        ff.numTotalTypes = 1
        ff.numMovementTypes = 1
        ff.ljEps = [1.0]
        ff.ljSigma = [0.3]
        state.forcefield = ff
        
        # Create atoms
        positions = [
            (1.0, 1.0, 1.0),
            (2.0, 1.5, 1.0),
            (1.5, 2.0, 1.5)
        ]
        
        state.atoms = []
        state.residues = []
        
        for i, (x, y, z) in enumerate(positions):
            atom = pygcmc.MCAtom()
            atom.type = 0
            atom.x, atom.y, atom.z = x, y, z
            atom.charge = 0.0
            state.atoms.append(atom)
            
            residue = pygcmc.MCResidue()
            residue.active = True
            residue.atomStart = i
            residue.atomCount = 1
            residue.atoms = [atom]
            state.residues.append(residue)
        
        state.activeAtomCount = len(state.atoms)
        state.activeResidueCount = len(state.residues)
        
        # Calculate initial energy
        pygcmc.computeSystemEnergyPBCCutoff(state)
        initial_energy = sum(r.energy_vdw + r.energy_elec for r in state.residues if r.active)
        
        # Mirror in x (x -> -x)
        for atom in state.atoms:
            atom.x = -atom.x
        for residue in state.residues:
            for atom in residue.atoms:
                atom.x = -atom.x
            
        pygcmc.computeSystemEnergyPBCCutoff(state)
        mirrored_energy = sum(r.energy_vdw + r.energy_elec for r in state.residues if r.active)
        
        # For LJ interactions without charges, mirroring should preserve energy
        assert abs(mirrored_energy - initial_energy) < 1e-6, \
            f"Energy changed after mirroring: {initial_energy:.6f} -> {mirrored_energy:.6f}"
        
        print(f"✓ Mirror symmetry test passed")
        print(f"  Energy remains {initial_energy:.6f} kJ/mol under mirroring")
    
    def test_cutoff_continuity(self):
        """Test that energy is continuous near cutoff"""
        state = pygcmc.MCState()
        state.info.box = (10.0, 10.0, 10.0)
        state.info.setTemperature(300.0)
        cutoff = 1.2  # nm
        state.info.cutoff = cutoff
        
        ff = pygcmc.MCForceField()
        ff.numTotalTypes = 1
        ff.numMovementTypes = 1
        ff.ljEps = [1.0]
        ff.ljSigma = [0.3]
        state.forcefield = ff
        
        # Two atoms at variable distance
        atom1 = pygcmc.MCAtom()
        atom1.type = 0
        atom1.x, atom1.y, atom1.z = 0.0, 0.0, 0.0
        atom1.charge = 0.0
        
        atom2 = pygcmc.MCAtom()
        atom2.type = 0
        atom2.charge = 0.0
        
        state.atoms = [atom1, atom2]
        
        # Create two residues
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
        
        # Sample energy near cutoff
        distances = np.linspace(cutoff - 0.1, cutoff + 0.1, 21)
        energies = []
        
        for r in distances:
            atom2.x = r
            atom2.y = 0.0
            atom2.z = 0.0
            state.atoms[1] = atom2
            res2.atoms[0] = atom2
            
            pygcmc.computeSystemEnergyCutoff(state)
            total_energy = res1.energy_vdw + res2.energy_vdw
            energies.append(total_energy)
        
        energies = np.array(energies)
        
        # Check that energy goes to zero at cutoff
        idx_cutoff = np.argmin(np.abs(distances - cutoff))
        energy_at_cutoff = energies[idx_cutoff]
        energy_beyond_cutoff = energies[idx_cutoff + 1:]
        
        # Energy should be zero beyond cutoff
        assert all(abs(e) < 1e-10 for e in energy_beyond_cutoff), \
            f"Energy not zero beyond cutoff: {energy_beyond_cutoff}"
        
        # Energy should be continuous (no sudden jump at cutoff)
        # This might not be perfectly smooth without switching function
        energy_before_cutoff = energies[idx_cutoff - 1]
        if abs(energy_before_cutoff) > 1e-10:
            # Check for reasonable continuity (not a huge jump)
            jump = abs(energy_at_cutoff - energy_before_cutoff)
            assert jump < abs(energy_before_cutoff), \
                f"Discontinuous jump at cutoff: {energy_before_cutoff:.6f} -> {energy_at_cutoff:.6f}"
        
        print(f"✓ Cutoff continuity test passed")
        print(f"  Energy at r=cutoff-0.05: {energies[idx_cutoff-1]:.6f} kJ/mol")
        print(f"  Energy at r=cutoff: {energy_at_cutoff:.6f} kJ/mol")
        print(f"  Energy at r=cutoff+0.05: {energies[idx_cutoff+1]:.6f} kJ/mol")
    
    def test_minimum_image_convention(self):
        """Test that minimum image convention is correctly applied"""
        state = pygcmc.MCState()
        box_size = 3.0
        state.info.box = (box_size, box_size, box_size)
        state.info.setTemperature(300.0)
        state.info.cutoff = box_size / 2 - 0.1  # Less than half box
        
        ff = pygcmc.MCForceField()
        ff.numTotalTypes = 1
        ff.numMovementTypes = 1
        ff.ljEps = [1.0]
        ff.ljSigma = [0.3]
        state.forcefield = ff
        
        # Two atoms: one at origin, one that we'll move
        atom1 = pygcmc.MCAtom()
        atom1.type = 0
        atom1.x, atom1.y, atom1.z = 0.0, 0.0, 0.0
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
        
        # Place atom2 at x = 0.5 (direct distance)
        atom2.x, atom2.y, atom2.z = 0.5, 0.0, 0.0
        state.atoms[1] = atom2
        res2.atoms[0] = atom2
        
        pygcmc.computeSystemEnergyPBCCutoff(state)
        energy_direct = res1.energy_vdw + res2.energy_vdw
        
        # Place atom2 at x = 2.5 (should use image at -0.5 via PBC)
        atom2.x = 2.5  # Direct distance is 2.5, but minimum image is 0.5
        state.atoms[1] = atom2
        res2.atoms[0] = atom2
        
        pygcmc.computeSystemEnergyPBCCutoff(state)
        energy_wrapped = res1.energy_vdw + res2.energy_vdw
        
        # Energies should be the same (same minimum image distance)
        assert abs(energy_direct - energy_wrapped) < 1e-6, \
            f"Minimum image not working: {energy_direct:.6f} != {energy_wrapped:.6f}"
        
        print(f"✓ Minimum image convention test passed")
        print(f"  Energy at x=0.5: {energy_direct:.6f} kJ/mol")
        print(f"  Energy at x=2.5 (wrapped to -0.5): {energy_wrapped:.6f} kJ/mol")

if __name__ == "__main__":
    if PYGCMC_AVAILABLE:
        print("Running PBC equivalence tests...\n")
        
        test = TestPBCEquivalence()
        test.test_translation_invariance()
        print()
        test.test_mirror_symmetry()
        print()
        test.test_cutoff_continuity()
        print()
        test.test_minimum_image_convention()
        
        print("\n✅ All PBC equivalence tests passed!")
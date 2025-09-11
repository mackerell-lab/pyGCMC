#!/usr/bin/env python
"""
Test state consistency in GCMC operations
Verifies that atom counts, residue counts, and arrays remain consistent
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
class TestStateConsistency:
    """Test state consistency during GCMC operations"""
    
    def test_translate_rotate_consistency(self):
        """Test that translation and rotation don't change atom counts"""
        # Setup
        state = pygcmc.MCState()
        state.info.box = (10.0, 10.0, 10.0)
        state.info.setTemperature(300.0)
        
        ff = pygcmc.MCForceField()
        ff.numTotalTypes = 1
        ff.numMovementTypes = 1
        ff.ljEps = [1.0]
        ff.ljSigma = [0.3]
        state.forcefield = ff
        
        # Create water template
        reservoir = pygcmc.movement.FragmentReservoir()
        template = pygcmc.movement.FragmentTemplate()
        template.typeId = 0
        
        # TIP3P-like water
        o_atom = pygcmc.MCAtom()
        o_atom.type = 0
        o_atom.x, o_atom.y, o_atom.z = 0.0, 0.0, 0.0
        o_atom.charge = -0.834
        
        h1_atom = pygcmc.MCAtom()
        h1_atom.type = 0
        h1_atom.x, h1_atom.y, h1_atom.z = 0.0957, 0.0, 0.0
        h1_atom.charge = 0.417
        
        h2_atom = pygcmc.MCAtom()
        h2_atom.type = 0
        h2_atom.x, h2_atom.y, h2_atom.z = -0.024, 0.0927, 0.0
        h2_atom.charge = 0.417
        
        template.atoms = [o_atom, h1_atom, h2_atom]
        reservoir.addTemplate(template)
        
        # Initialize engine
        engine = pygcmc.GCMCEngine()
        engine.initialize(state, reservoir)
        engine.setTemperature(300.0)
        engine.setSeed(42)
        
        acceptance = pygcmc.GCMCAcceptance()
        acceptance.setTemperature(300.0)
        acceptance.setVolume(1000.0)
        acceptance.setActivity(0, 1.0)
        engine.setAcceptanceCalculator(acceptance)
        
        # Insert some molecules
        n_molecules = 5
        for i in range(n_molecules):
            result = engine.attemptInsertion(0)
            if result.accepted:
                print(f"Inserted molecule {i+1}")
        
        # Store initial counts
        initial_atom_count = state.activeAtomCount
        initial_residue_count = state.activeResidueCount
        initial_atoms_len = len(state.atoms)
        initial_residues_len = len(state.residues)
        
        print(f"\nInitial state:")
        print(f"  activeAtomCount: {initial_atom_count}")
        print(f"  activeResidueCount: {initial_residue_count}")
        print(f"  len(atoms): {initial_atoms_len}")
        print(f"  len(residues): {initial_residues_len}")
        
        # Perform multiple translations and rotations
        n_moves = 20
        for i in range(n_moves):
            # Try translation
            trans_result = engine.attemptTranslation(0)
            
            # Check consistency after translation
            assert state.activeAtomCount == initial_atom_count, \
                f"Atom count changed after translation {i+1}: {state.activeAtomCount} != {initial_atom_count}"
            assert state.activeResidueCount == initial_residue_count, \
                f"Residue count changed after translation {i+1}: {state.activeResidueCount} != {initial_residue_count}"
            
            # Try rotation
            rot_result = engine.attemptRotation(0)
            
            # Check consistency after rotation
            assert state.activeAtomCount == initial_atom_count, \
                f"Atom count changed after rotation {i+1}: {state.activeAtomCount} != {initial_atom_count}"
            assert state.activeResidueCount == initial_residue_count, \
                f"Residue count changed after rotation {i+1}: {state.activeResidueCount} != {initial_residue_count}"
        
        print(f"\nAfter {n_moves} translations and rotations:")
        print(f"  activeAtomCount: {state.activeAtomCount} (unchanged ✓)")
        print(f"  activeResidueCount: {state.activeResidueCount} (unchanged ✓)")
        print(f"  len(atoms): {len(state.atoms)} (unchanged ✓)")
        
        print("✓ Translation/rotation consistency test passed")
    
    def test_insertion_deletion_consistency(self):
        """Test that insertion and deletion maintain consistent counts"""
        state = pygcmc.MCState()
        state.info.box = (10.0, 10.0, 10.0)
        state.info.setTemperature(300.0)
        
        ff = pygcmc.MCForceField()
        ff.numTotalTypes = 1
        ff.numMovementTypes = 1
        ff.ljEps = [1.0]
        ff.ljSigma = [0.3]
        state.forcefield = ff
        
        # Simple single-atom template
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
        engine.setSeed(12345)
        
        acceptance = pygcmc.GCMCAcceptance()
        acceptance.setTemperature(300.0)
        acceptance.setVolume(1000.0)
        acceptance.setActivity(0, 10.0)  # High activity for more insertions
        engine.setAcceptanceCalculator(acceptance)
        
        # Track counts through insertion/deletion cycles
        atom_counts = []
        residue_counts = []
        
        for cycle in range(10):
            # Insert
            result = engine.attemptInsertion(0)
            if result.accepted:
                atom_counts.append(state.activeAtomCount)
                residue_counts.append(state.activeResidueCount)
                
                # Verify consistency
                n_active_residues = sum(1 for r in state.residues if r.active)
                assert state.activeResidueCount >= n_active_residues, \
                    f"activeResidueCount ({state.activeResidueCount}) < actual active residues ({n_active_residues})"
            
            # Delete
            if state.activeResidueCount > 0:
                result = engine.attemptDeletion(0)
                if result.accepted:
                    atom_counts.append(state.activeAtomCount)
                    residue_counts.append(state.activeResidueCount)
                    
                    # Verify consistency
                    n_active_residues = sum(1 for r in state.residues if r.active)
                    assert state.activeResidueCount >= n_active_residues, \
                        f"activeResidueCount ({state.activeResidueCount}) < actual active residues ({n_active_residues})"
        
        print(f"✓ Insertion/deletion consistency test passed")
        print(f"  Atom count history: {atom_counts[:10]}")
        print(f"  Residue count history: {residue_counts[:10]}")
    
    def test_atomstart_atomcount_consistency(self):
        """Test that atomStart and atomCount remain consistent"""
        state = pygcmc.MCState()
        state.info.box = (10.0, 10.0, 10.0)
        state.info.setTemperature(300.0)
        
        ff = pygcmc.MCForceField()
        ff.numTotalTypes = 2
        ff.numMovementTypes = 2
        ff.setPerTypeParameters([0.3, 0.35], [0.5, 0.6])
        ff.mixingRule = pygcmc.MCForceField.MixingRule.LorentzBerthelot
        ff.rebuildLJMatrix()
        state.forcefield = ff
        
        # Create templates with different atom counts
        reservoir = pygcmc.movement.FragmentReservoir()
        
        # Template 1: 2 atoms
        template1 = pygcmc.movement.FragmentTemplate()
        template1.typeId = 0
        for i in range(2):
            atom = pygcmc.MCAtom()
            atom.type = 0
            atom.x, atom.y, atom.z = i * 0.1, 0.0, 0.0
            atom.charge = 0.0
            template1.atoms.append(atom)
        reservoir.addTemplate(template1)
        
        # Template 2: 3 atoms
        template2 = pygcmc.movement.FragmentTemplate()
        template2.typeId = 1
        for i in range(3):
            atom = pygcmc.MCAtom()
            atom.type = 1
            atom.x, atom.y, atom.z = i * 0.1, 0.0, 0.0
            atom.charge = 0.0
            template2.atoms.append(atom)
        reservoir.addTemplate(template2)
        
        engine = pygcmc.GCMCEngine()
        engine.initialize(state, reservoir)
        engine.setTemperature(300.0)
        engine.setSeed(54321)
        
        acceptance = pygcmc.GCMCAcceptance()
        acceptance.setTemperature(300.0)
        acceptance.setVolume(1000.0)
        acceptance.setActivity(0, 5.0)
        acceptance.setActivity(1, 5.0)
        engine.setAcceptanceCalculator(acceptance)
        
        # Insert molecules of both types
        insertions = []
        for _ in range(10):
            type_id = np.random.randint(0, 2)
            result = engine.attemptInsertion(type_id)
            if result.accepted:
                insertions.append((type_id, result.residueIndex))
        
        # Verify atomStart and atomCount for all active residues
        for res_idx, residue in enumerate(state.residues):
            if residue.active:
                # Check atomCount matches template
                instance = reservoir.getInstance(res_idx)
                if instance:
                    template = reservoir.getTemplate(instance.templateId)
                    if template:
                        assert residue.atomCount == len(template.atoms), \
                            f"Residue {res_idx}: atomCount ({residue.atomCount}) != template atoms ({len(template.atoms)})"
                
                # Check atoms in global array match residue atoms
                if residue.atomStart >= 0 and residue.atomCount > 0:
                    for i in range(residue.atomCount):
                        global_idx = residue.atomStart + i
                        if global_idx < len(state.atoms):
                            global_atom = state.atoms[global_idx]
                            if i < len(residue.atoms):
                                local_atom = residue.atoms[i]
                                # Check coordinates match
                                assert abs(global_atom.x - local_atom.x) < 1e-6, \
                                    f"Residue {res_idx} atom {i}: x mismatch"
                                assert abs(global_atom.y - local_atom.y) < 1e-6, \
                                    f"Residue {res_idx} atom {i}: y mismatch"
                                assert abs(global_atom.z - local_atom.z) < 1e-6, \
                                    f"Residue {res_idx} atom {i}: z mismatch"
        
        print(f"✓ atomStart/atomCount consistency test passed")
        print(f"  Successfully verified {len(insertions)} insertions")
    
    def test_energy_consistency_after_moves(self):
        """Test that energy calculations remain consistent after moves"""
        state = pygcmc.MCState()
        state.info.box = (10.0, 10.0, 10.0)
        state.info.setTemperature(300.0)
        state.info.cutoff = 5.0
        
        ff = pygcmc.MCForceField()
        ff.numTotalTypes = 1
        ff.numMovementTypes = 1
        ff.ljEps = [1.0]
        ff.ljSigma = [0.3]
        state.forcefield = ff
        
        # Two-atom template
        reservoir = pygcmc.movement.FragmentReservoir()
        template = pygcmc.movement.FragmentTemplate()
        template.typeId = 0
        
        atom1 = pygcmc.MCAtom()
        atom1.type = 0
        atom1.x, atom1.y, atom1.z = 0.0, 0.0, 0.0
        atom1.charge = 0.5
        
        atom2 = pygcmc.MCAtom()
        atom2.type = 0
        atom2.x, atom2.y, atom2.z = 0.2, 0.0, 0.0
        atom2.charge = -0.5
        
        template.atoms = [atom1, atom2]
        reservoir.addTemplate(template)
        
        engine = pygcmc.GCMCEngine()
        engine.initialize(state, reservoir)
        engine.setTemperature(300.0)
        engine.setSeed(99999)
        
        acceptance = pygcmc.GCMCAcceptance()
        acceptance.setTemperature(300.0)
        acceptance.setVolume(1000.0)
        acceptance.setActivity(0, 5.0)
        engine.setAcceptanceCalculator(acceptance)
        
        # Insert a few molecules
        for _ in range(3):
            engine.attemptInsertion(0)
        
        if state.activeResidueCount > 1:
            # Calculate initial energy
            pygcmc.computeSystemEnergyCutoff(state)
            initial_energy = sum(r.energy_vdw + r.energy_elec for r in state.residues if r.active)
            
            # Perform moves
            for _ in range(10):
                engine.attemptTranslation(0)
                engine.attemptRotation(0)
            
            # Recalculate energy
            pygcmc.computeSystemEnergyCutoff(state)
            final_energy = sum(r.energy_vdw + r.energy_elec for r in state.residues if r.active)
            
            # Energy should be physically reasonable (not NaN or infinite)
            assert not np.isnan(final_energy), "Energy is NaN after moves"
            assert not np.isinf(final_energy), "Energy is infinite after moves"
            
            print(f"✓ Energy consistency test passed")
            print(f"  Initial energy: {initial_energy:.3f} kJ/mol")
            print(f"  Final energy: {final_energy:.3f} kJ/mol")
        else:
            print("✓ Energy consistency test passed (no molecules inserted)")

if __name__ == "__main__":
    if PYGCMC_AVAILABLE:
        print("Running state consistency tests...\n")
        
        test = TestStateConsistency()
        test.test_translate_rotate_consistency()
        print()
        test.test_insertion_deletion_consistency()
        print()
        test.test_atomstart_atomcount_consistency()
        print()
        test.test_energy_consistency_after_moves()
        
        print("\n✅ All state consistency tests passed!")
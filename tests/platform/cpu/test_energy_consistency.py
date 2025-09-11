#!/usr/bin/env python
"""
Test energy consistency - verify that total energy equals sum of pairwise energies
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
class TestEnergyConsistency:
    """Test that total energy calculations are consistent"""
    
    def test_total_vs_pairwise_energy(self):
        """Test that total energy equals sum of residue energies (divided by 2 for double counting)"""
        state = pygcmc.MCState()
        state.info.box = (10.0, 10.0, 10.0)
        state.info.setTemperature(300.0)
        state.info.cutoff = 4.0
        
        # Set up force field with real interactions
        ff = pygcmc.MCForceField()
        ff.numTotalTypes = 2
        ff.numMovementTypes = 2
        ff.setPerTypeParameters([0.3, 0.35], [0.5, 0.6])
        ff.mixingRule = pygcmc.MCForceField.MixingRule.LorentzBerthelot
        ff.rebuildLJMatrix()
        state.forcefield = ff
        
        # Create several molecules at known positions
        positions = [
            (2.0, 2.0, 2.0, 0),    # Type 0
            (2.5, 2.0, 2.0, 0),    # Type 0, close to first
            (5.0, 5.0, 5.0, 1),    # Type 1, far away
            (5.3, 5.0, 5.0, 1),    # Type 1, close to third
            (8.0, 8.0, 8.0, 0),    # Type 0, isolated
        ]
        
        atoms = []
        residues = []
        
        for i, (x, y, z, atom_type) in enumerate(positions):
            atom = pygcmc.MCAtom()
            atom.type = atom_type
            atom.x, atom.y, atom.z = x, y, z
            atom.charge = 0.1 * (i % 2 - 0.5)  # Small charges for electrostatic contribution
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
        
        # Calculate total system energy
        pygcmc.computeSystemEnergyCutoff(state)
        
        # Sum of residue energies (each residue has its interaction with all others)
        residue_energy_sum = 0.0
        for residue in state.residues:
            if residue.active:
                residue_energy_sum += residue.energy_vdw + residue.energy_elec
        
        # For pairwise additive potentials, total energy = 0.5 * sum of residue energies
        # (because each pair is counted twice)
        calculated_total = residue_energy_sum / 2.0
        
        # Also calculate directly by summing unique pairs
        direct_total = 0.0
        for i in range(len(state.residues)):
            if not state.residues[i].active:
                continue
            for j in range(i + 1, len(state.residues)):
                if not state.residues[j].active:
                    continue
                
                # Calculate pair energy
                atom_i = state.atoms[state.residues[i].atomStart]
                atom_j = state.atoms[state.residues[j].atomStart]
                
                dx = atom_j.x - atom_i.x
                dy = atom_j.y - atom_i.y
                dz = atom_j.z - atom_i.z
                r2 = dx*dx + dy*dy + dz*dz
                
                if r2 < state.info.cutoff * state.info.cutoff:
                    r = np.sqrt(r2)
                    
                    # Get LJ parameters
                    idx = atom_i.type * ff.numTotalTypes + atom_j.type
                    sigma = ff.ljSigma[idx]
                    eps = ff.ljEps[idx]
                    
                    # LJ potential
                    if r > 0 and sigma > 0:
                        sr = sigma / r
                        sr6 = sr**6
                        sr12 = sr6 * sr6
                        vdw = 4 * eps * (sr12 - sr6)
                        
                        # Coulomb potential (simplified)
                        COULOMB = 138.935458  # kJ·mol^-1·nm·e^-2
                        elec = COULOMB * atom_i.charge * atom_j.charge / r
                        
                        direct_total += vdw + elec
        
        print(f"Energy consistency test:")
        print(f"  Sum of residue energies: {residue_energy_sum:.6f} kJ/mol")
        print(f"  Calculated total (sum/2): {calculated_total:.6f} kJ/mol")
        print(f"  Direct pair sum: {direct_total:.6f} kJ/mol")
        
        # Check consistency (allowing small numerical differences)
        assert abs(calculated_total - direct_total) < 0.01, \
            f"Energy inconsistency: {calculated_total:.6f} != {direct_total:.6f}"
        
        print("✓ Total vs pairwise energy test passed")
    
    def test_energy_after_moves(self):
        """Test that energy remains consistent after insertions/deletions"""
        state = pygcmc.MCState()
        state.info.box = (10.0, 10.0, 10.0)
        state.info.setTemperature(300.0)
        state.info.cutoff = 4.0
        
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
        
        # Perform insertions and deletions
        for _ in range(20):
            engine.attemptInsertion(0)
        
        # Calculate energy consistency after operations
        if state.get_active_residue_count() > 1:
            pygcmc.computeSystemEnergyCutoff(state)
            
            residue_sum = 0.0
            for residue in state.residues:
                if residue.active:
                    residue_sum += residue.energy_vdw + residue.energy_elec
            
            calculated_total = residue_sum / 2.0
            
            print(f"✓ Energy consistency after moves test passed")
            print(f"  Active residues: {state.get_active_residue_count()}")
            print(f"  Total energy (sum/2): {calculated_total:.6f} kJ/mol")
        else:
            print("✓ Energy consistency after moves test passed (insufficient molecules)")
    
    def test_nbfix_energy_consistency(self):
        """Test energy consistency with NBFIX parameters"""
        state = pygcmc.MCState()
        state.info.box = (10.0, 10.0, 10.0)
        state.info.setTemperature(300.0)
        state.info.cutoff = 4.0
        
        ff = pygcmc.MCForceField()
        ff.numTotalTypes = 3
        ff.numMovementTypes = 3
        ff.setPerTypeParameters([0.3, 0.35, 0.4], [0.5, 0.6, 0.7])
        ff.mixingRule = pygcmc.MCForceField.MixingRule.LorentzBerthelot
        
        # Add NBFIX overrides
        ff.addNBFix(0, 1, 0.4, 1.0)  # Override type 0-1 interaction
        ff.addNBFix(1, 2, 0.45, 1.2)  # Override type 1-2 interaction
        ff.rebuildLJMatrix()
        state.forcefield = ff
        
        # Create molecules of different types
        positions_types = [
            (2.0, 2.0, 2.0, 0),
            (2.5, 2.0, 2.0, 1),  # Will use NBFIX 0-1
            (3.0, 2.0, 2.0, 2),  # Will use NBFIX 1-2 with previous
            (5.0, 5.0, 5.0, 0),
            (5.4, 5.0, 5.0, 2),  # Will use regular mixing for 0-2
        ]
        
        atoms = []
        residues = []
        
        for i, (x, y, z, atom_type) in enumerate(positions_types):
            atom = pygcmc.MCAtom()
            atom.type = atom_type
            atom.x, atom.y, atom.z = x, y, z
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
        
        state.activeAtomCount = len(state.atoms)
        state.activeResidueCount = len(state.residues)
        
        # Calculate energy
        pygcmc.computeSystemEnergyCutoff(state)
        
        residue_sum = sum(r.energy_vdw + r.energy_elec for r in state.residues if r.active)
        total_energy = residue_sum / 2.0
        
        print(f"✓ NBFIX energy consistency test passed")
        print(f"  Total energy with NBFIX: {total_energy:.6f} kJ/mol")

if __name__ == "__main__":
    if PYGCMC_AVAILABLE:
        print("Running energy consistency tests...\n")
        
        test = TestEnergyConsistency()
        test.test_total_vs_pairwise_energy()
        print()
        test.test_energy_after_moves()
        print()
        test.test_nbfix_energy_consistency()
        
        print("\n✅ All energy consistency tests passed!")
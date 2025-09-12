#!/usr/bin/env python
"""
Test trajectory writer functionality
"""

import pytest
import numpy as np
import os
import sys
import tempfile

# Add build path for pygcmc module
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../build'))

try:
    import pygcmc
    PYGCMC_AVAILABLE = True
except ImportError:
    PYGCMC_AVAILABLE = False
    pygcmc = None

@pytest.mark.skipif(not PYGCMC_AVAILABLE, reason="PyGCMC not available")
class TestTrajectoryWriter:
    """Test trajectory writing functionality"""
    
    def setup_test_state(self):
        """Create a test MCState with some molecules"""
        state = pygcmc.MCState()
        state.info.box = (5.0, 5.0, 5.0)
        state.info.setTemperature(300.0)
        
        # Setup force field
        ff = pygcmc.MCForceField()
        ff.numTotalTypes = 2
        ff.numMovementTypes = 2
        ff.ljEps = [1.0, 0.5]
        ff.ljSigma = [0.3, 0.35]
        state.forcefield = ff
        
        # Note: atomTypes and residueTypes are managed internally by TypeMaps
        # We'll use type indices directly
        
        # Create some molecules
        atoms = []
        residues = []
        
        # Molecule 1
        atom1 = pygcmc.MCAtom()
        atom1.type = 0
        atom1.x, atom1.y, atom1.z = 1.0, 1.0, 1.0
        atom1.charge = 0.0
        atoms.append(atom1)
        
        atom2 = pygcmc.MCAtom()
        atom2.type = 1
        atom2.x, atom2.y, atom2.z = 1.5, 1.0, 1.0
        atom2.charge = -0.5
        atoms.append(atom2)
        
        res1 = pygcmc.MCResidue()
        res1.active = True
        res1.atomStart = 0
        res1.atomCount = 2
        res1.atoms = [atom1, atom2]
        res1.energy_vdw = -1.5
        res1.energy_elec = -0.5
        residues.append(res1)
        
        # Molecule 2
        atom3 = pygcmc.MCAtom()
        atom3.type = 0
        atom3.x, atom3.y, atom3.z = 3.0, 3.0, 3.0
        atom3.charge = 0.0
        atoms.append(atom3)
        
        res2 = pygcmc.MCResidue()
        res2.active = True
        res2.atomStart = 2
        res2.atomCount = 1
        res2.atoms = [atom3]
        res2.energy_vdw = -0.8
        res2.energy_elec = 0.0
        residues.append(res2)
        
        state.atoms = atoms
        state.residues = residues
        state.activeAtomCount = 3
        state.activeResidueCount = 2
        
        return state
    
    def test_pdb_writer(self):
        """Test PDB format writing"""
        state = self.setup_test_state()
        
        # Note: TypeMaps doesn't support setTypeName in current implementation
        
        with tempfile.NamedTemporaryFile(suffix='.pdb', delete=False) as f:
            filename = f.name
        
        try:
            # Write PDB file
            writer = pygcmc.TrajectoryWriter(filename, pygcmc.TrajectoryWriter.Format.PDB)
            writer.write_frame(state, 0)
            writer.write_frame(state, 1)
            writer.close()
            
            # Check file was created and has content
            assert os.path.exists(filename)
            assert os.path.getsize(filename) > 0
            
            # Read and verify content
            with open(filename, 'r') as f:
                content = f.read()
                assert 'REMARK' in content
                assert 'MODEL' in content
                assert 'ATOM' in content
                assert 'CRYST1' in content
                assert 'ENDMDL' in content
                assert 'END' in content
                
                # Check box dimensions (converted to Angstrom)
                assert '50.000' in content  # 5.0 nm = 50.0 Angstrom
                
                # Check we have 2 frames
                assert content.count('MODEL') == 2
                assert content.count('ENDMDL') == 2
        finally:
            if os.path.exists(filename):
                os.remove(filename)
    
    def test_xyz_writer(self):
        """Test XYZ format writing"""
        state = self.setup_test_state()
        
        with tempfile.NamedTemporaryFile(suffix='.xyz', delete=False) as f:
            filename = f.name
        
        try:
            writer = pygcmc.TrajectoryWriter(filename, pygcmc.TrajectoryWriter.Format.XYZ)
            writer.write_frame(state, 0)
            writer.close()
            
            # Check file content
            with open(filename, 'r') as f:
                lines = f.readlines()
                assert lines[0].strip() == '3'  # 3 atoms
                assert 'Frame 0' in lines[1]
                assert 'Energy:' in lines[1]
                
                # Check atom lines
                atom_lines = lines[2:]
                assert len(atom_lines) == 3
                for line in atom_lines:
                    parts = line.split()
                    assert len(parts) == 4  # Element x y z
        finally:
            if os.path.exists(filename):
                os.remove(filename)
    
    def test_dat_writer(self):
        """Test DAT format writing"""
        state = self.setup_test_state()
        
        with tempfile.NamedTemporaryFile(suffix='.dat', delete=False) as f:
            filename = f.name
        
        try:
            writer = pygcmc.TrajectoryWriter(filename, pygcmc.TrajectoryWriter.Format.DAT)
            
            # Write multiple frames
            for i in range(5):
                writer.write_frame(state, i)
            
            writer.close()
            
            # Check file content
            with open(filename, 'r') as f:
                lines = f.readlines()
                assert len(lines) == 5
                
                for i, line in enumerate(lines):
                    parts = line.split()
                    assert len(parts) == 4  # frame, n_molecules, energy, volume
                    assert int(parts[0]) == i  # Frame number
                    assert int(parts[1]) == 2  # 2 molecules
                    assert float(parts[3]) == 125.0  # Volume = 5^3
        finally:
            if os.path.exists(filename):
                os.remove(filename)
    
    def test_top_writer(self):
        """Test TOP format writing"""
        state = self.setup_test_state()
        
        with tempfile.NamedTemporaryFile(suffix='.top', delete=False) as f:
            filename = f.name
        
        try:
            writer = pygcmc.TrajectoryWriter(filename, pygcmc.TrajectoryWriter.Format.TOP)
            writer.write_topology(state)
            writer.close()
            
            # Check file content
            with open(filename, 'r') as f:
                content = f.read()
                assert '[ atomtypes ]' in content
                assert '[ molecules ]' in content
                assert 'Box:' in content
                assert 'Temperature:' in content
                
                # Check atom types section
                assert 'sigma' in content
                assert 'epsilon' in content
        finally:
            if os.path.exists(filename):
                os.remove(filename)
    
    def test_data_writer(self):
        """Test DataWriter for analysis output"""
        with tempfile.NamedTemporaryFile(suffix='.dat', delete=False) as f:
            filename = f.name
        
        try:
            writer = pygcmc.DataWriter(filename)
            
            # Write header
            writer.write_header(["Frame", "N_molecules", "Energy", "Volume"])
            
            # Write some data
            for i in range(10):
                writer.write_row([i, i*2, -i*0.5, 125.0])
            
            # Write comment
            writer.write_comment("End of simulation")
            
            writer.close()
            
            # Check file content
            with open(filename, 'r') as f:
                lines = f.readlines()
                assert lines[0].startswith('#')  # Header
                assert 'Frame' in lines[0]
                assert 'N_molecules' in lines[0]
                
                # Check data lines
                for i in range(1, 11):
                    parts = lines[i].split()
                    assert len(parts) == 4
                    assert float(parts[0]) == i - 1
                
                # Check comment
                assert lines[-1].startswith('# End of simulation')
        finally:
            if os.path.exists(filename):
                os.remove(filename)
    
    def test_context_manager(self):
        """Test using writers as context managers"""
        state = self.setup_test_state()
        
        with tempfile.NamedTemporaryFile(suffix='.pdb', delete=False) as f:
            filename = f.name
        
        try:
            # Use with statement
            with pygcmc.TrajectoryWriter(filename) as writer:
                writer.write_frame(state)
            
            # File should be properly closed
            assert os.path.exists(filename)
            assert os.path.getsize(filename) > 0
            
            # Check it was properly closed
            with open(filename, 'r') as f:
                content = f.read()
                assert 'END' in content  # PDB end marker
        finally:
            if os.path.exists(filename):
                os.remove(filename)

if __name__ == "__main__":
    if PYGCMC_AVAILABLE:
        print("Running trajectory writer tests...")
        
        test = TestTrajectoryWriter()
        test.test_pdb_writer()
        print("✓ PDB writer test passed")
        
        test.test_xyz_writer()
        print("✓ XYZ writer test passed")
        
        test.test_dat_writer()
        print("✓ DAT writer test passed")
        
        test.test_top_writer()
        print("✓ TOP writer test passed")
        
        test.test_data_writer()
        print("✓ Data writer test passed")
        
        test.test_context_manager()
        print("✓ Context manager test passed")
        
        print("\n✅ All trajectory writer tests passed!")
    else:
        print("PyGCMC not available, skipping tests")
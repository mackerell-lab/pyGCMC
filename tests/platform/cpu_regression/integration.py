#!/usr/bin/env python3
"""
Integration tests for gcmc_cpu using real data files
These tests demonstrate end-to-end functionality
"""

import pytest
import subprocess
import os
import tempfile
from pathlib import Path

GCMC_CPU_PATH = Path(__file__).resolve().parents[3] / "build" / "bin" / "gcmc_cpu"
TEST_DATA_DIR = Path(__file__).resolve().parents[2] / "data"


class TestGCMCIntegration:
    """Integration tests demonstrating full GCMC workflows"""
    
    @pytest.fixture
    def temp_dir(self):
        """Create a temporary directory for test outputs"""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield tmpdir
    
    def test_water_box_fragitp_only(self, temp_dir):
        """Test water box simulation using only fragitp (no PDB/TOP needed)"""
        # Create INP file for water box with absolute paths
        fragitp_path = TEST_DATA_DIR / "charmm36.ff/mol/sol.itp"
        inp_content = f"""# Water box test - fragitp mode
inp_units:nm
box: 3.0 3.0 3.0
temperature: 298.15
fragname: SOL
fragmuex: -10.5
fragconc: 55.0
mcsteps: 50
nprint: 10
fragitp: {fragitp_path}
"""
        inp_path = Path(temp_dir) / "water_box.inp"
        inp_path.write_text(inp_content)
        
        # Run simulation with verbose for detailed output
        # Use full path for prefix to ensure output goes to temp_dir
        output_prefix = str(Path(temp_dir) / "water_box")
        result = subprocess.run(
            [str(GCMC_CPU_PATH), "--inp", str(inp_path), 
             "--prefix", output_prefix, "--seed", "42", "--verbose"],
            cwd=temp_dir,  # Run from temp dir to keep outputs there
            capture_output=True,
            text=True,
            timeout=5
        )
        
        # Combine stdout and stderr for checking
        combined_output = result.stdout + result.stderr
        
        # Check successful completion
        assert result.returncode == 0, f"Failed with: {result.stderr}"
        assert "Simulation completed" in result.stdout or "Simulation completed" in combined_output
        
        # Verify key outputs (in verbose mode these appear in stderr)
        assert "Temperature: 298.15 K" in combined_output
        assert "Box size: 3 x 3 x 3 nm" in combined_output
        assert "Fragment SOL" in combined_output
        
        # Check molecule insertion occurred
        assert "SOL:" in result.stdout  # Fragment count output
        
    def test_multi_fragment_salt_water(self, temp_dir):
        """Test multi-component simulation with water and ions"""
        # Create INP file with multiple fragments using absolute paths
        sol_itp = TEST_DATA_DIR / "charmm36.ff/mol/sol.itp"
        na_itp = TEST_DATA_DIR / "charmm36.ff/mol/na.itp"
        cl_itp = TEST_DATA_DIR / "charmm36.ff/mol/cl.itp"
        inp_content = f"""# Salt water test - multiple fragments
inp_units:nm
box: 4.0 4.0 4.0
temperature: 298.15
fragname: SOL
fragmuex: -10.5
fragconc: 55.0
fragname: NA
fragmuex: -8.0
fragconc: 0.15
fragname: CL
fragmuex: -7.5
fragconc: 0.15
mcsteps: 100
nprint: 20
fragitp: {sol_itp}
fragitp: {na_itp}
fragitp: {cl_itp}
"""
        inp_path = Path(temp_dir) / "salt_water.inp"
        inp_path.write_text(inp_content)
        
        # Run simulation with output to temp_dir
        output_prefix = str(Path(temp_dir) / "salt")
        result = subprocess.run(
            [str(GCMC_CPU_PATH), "--inp", str(inp_path),
             "--prefix", output_prefix, "--seed", "123", "--verbose"],
            cwd=temp_dir,
            capture_output=True,
            text=True,
            timeout=10
        )
        
        # Check completion
        assert result.returncode == 0
        combined_output = result.stdout + result.stderr
        assert "Simulation completed" in result.stdout or "Simulation completed" in combined_output
        
        # Verify all fragments were loaded (may be in stderr in verbose mode)  
        # Check for lowercase fragment names as they appear in the logs
        assert "fragment sol" in combined_output.lower() or "sol:" in combined_output.lower()
        assert "fragment na" in combined_output.lower() or "na:" in combined_output.lower()  
        assert "fragment cl" in combined_output.lower() or "cl:" in combined_output.lower()
        
        # Check that fragment counts are reported
        output_lines = result.stdout.split('\n')
        fragment_count_lines = [l for l in output_lines if 'Fragment counts:' in l]
        assert len(fragment_count_lines) > 0
    
    def test_with_forcefield_parameters(self, temp_dir):
        """Test with force field parameters from PAR files"""
        fragitp_path = TEST_DATA_DIR / "charmm36.ff/mol/sol.itp"
        parfile_path = TEST_DATA_DIR / "charmm36.ff/toppar_water_ions.str"
        inp_content = f"""# Test with force field parameters
inp_units:nm
box: 3.0 3.0 3.0
temperature: 300.0
fragname: SOL
fragmuex: -10.0
fragconc: 50.0
mcsteps: 30
nprint: 10
fragitp: {fragitp_path}
parfile: {parfile_path}
"""
        inp_path = Path(temp_dir) / "with_forcefield.inp"
        inp_path.write_text(inp_content)
        
        # Run simulation with output to temp_dir
        output_prefix = str(Path(temp_dir) / "ff_test")
        result = subprocess.run(
            [str(GCMC_CPU_PATH), "--inp", str(inp_path),
             "--prefix", output_prefix, "--seed", "42", "--verbose"],
            cwd=temp_dir,
            capture_output=True,
            text=True,
            timeout=5
        )
        
        # Check for force field loading messages
        if "toppar_water_ions.str" in os.listdir(TEST_DATA_DIR / "charmm36.ff"):
            assert "force field parameters" in result.stdout.lower() or \
                   "Applying force field" in result.stdout
        
        # Basic completion check
        assert result.returncode == 0
        assert "mcsteps: 30" in result.stdout.lower() or "MC steps: 30" in result.stdout
    
    @pytest.mark.parametrize("box_size,expected_volume", [
        ([2.0, 2.0, 2.0], 8.0),
        ([3.0, 3.0, 3.0], 27.0),
        ([2.0, 3.0, 4.0], 24.0),
    ])
    def test_different_box_sizes(self, temp_dir, box_size, expected_volume):
        """Test various box sizes and verify volume calculation"""
        fragitp_path = TEST_DATA_DIR / "charmm36.ff/mol/sol.itp"
        inp_content = f"""# Box size test
inp_units:nm
box: {box_size[0]} {box_size[1]} {box_size[2]}
temperature: 298.15
fragname: SOL
fragmuex: -10.5
fragconc: 10.0
mcsteps: 10
nprint: 5
fragitp: {fragitp_path}
"""
        inp_path = Path(temp_dir) / "box_test.inp"
        inp_path.write_text(inp_content)
        
        output_prefix = str(Path(temp_dir) / "box")
        result = subprocess.run(
            [str(GCMC_CPU_PATH), "--inp", str(inp_path),
             "--prefix", output_prefix, "--seed", "42", "--verbose"],
            cwd=temp_dir,
            capture_output=True,
            text=True,
            timeout=5
        )
        
        assert result.returncode == 0
        # Check box size in output (may be in stderr)
        combined_output = result.stdout + result.stderr
        # Check for box size in various formats
        box_str1 = f"Box size: {box_size[0]} x {box_size[1]} x {box_size[2]} nm"
        box_str2 = f"Box size: {int(box_size[0])} x {int(box_size[1])} x {int(box_size[2])} nm"
        assert box_str1 in combined_output or box_str2 in combined_output, f"Box size not found. Output: {combined_output[:1000]}..."
        
        # Check volume - allow integer or float format
        volume_str1 = f"Box volume: {expected_volume} nm³"
        volume_str2 = f"Box volume: {int(expected_volume)} nm³"
        assert volume_str1 in combined_output or volume_str2 in combined_output
    
    def test_fragitp_only_mode(self, temp_dir):
        """Test simulation can run with only fragitp files (no PDB/TOP)"""
        fragitp_path = TEST_DATA_DIR / "charmm36.ff/mol/sol.itp"
        inp_content = f"""# Test fragitp-only mode
inp_units:nm
box: 3.0 3.0 3.0
temperature: 298.15
fragname: SOL
fragmuex: -10.5
fragconc: 55.0
mcsteps: 20
nprint: 10
fragitp: {fragitp_path}
"""
        inp_path = Path(temp_dir) / "fragitp_only.inp"
        inp_path.write_text(inp_content)
        
        output_prefix = str(Path(temp_dir) / "fragitp_only")
        result = subprocess.run(
            [str(GCMC_CPU_PATH), "--inp", str(inp_path),
             "--prefix", output_prefix, "--seed", "42", "--verbose"],
            cwd=temp_dir,
            capture_output=True,
            text=True,
            timeout=5
        )
        
        # Should complete successfully without PDB/TOP files
        assert result.returncode == 0, f"Failed with: {result.stdout}\n{result.stderr}"
        assert "Simulation completed" in result.stdout
        
        # Check that fragment was loaded from ITP
        combined = result.stdout + result.stderr
        assert "fragment sol" in combined.lower() or "loaded fragment" in combined.lower()
    
    def test_deterministic_with_seed(self, temp_dir):
        """Test that same seed produces deterministic results"""
        fragitp_path = TEST_DATA_DIR / "charmm36.ff/mol/sol.itp"
        inp_content = f"""# Deterministic test
inp_units:nm
box: 2.5 2.5 2.5
temperature: 298.15
fragname: SOL
fragmuex: -10.5
fragconc: 30.0
mcsteps: 50
nprint: 50
fragitp: {fragitp_path}
"""
        inp_path = Path(temp_dir) / "deterministic.inp"
        inp_path.write_text(inp_content)
        
        # Run twice with same seed
        results = []
        for run in range(2):
            output_prefix = str(Path(temp_dir) / f"det_{run}")
            result = subprocess.run(
                [str(GCMC_CPU_PATH), "--inp", str(inp_path),
                 "--prefix", output_prefix, "--seed", "999"],
                cwd=temp_dir,
                capture_output=True,
                text=True,
                timeout=5
            )
            results.append(result.stdout)
        
        # Extract final molecule counts from both runs
        def extract_final_count(output):
            lines = output.split('\n')
            # Look for any line with "SOL:", even if 0 molecules
            sol_counts = []
            for line in lines:
                if "SOL:" in line:
                    import re
                    match = re.search(r'SOL:\s*(\d+)', line)
                    if match:
                        sol_counts.append(int(match.group(1)))
            # Return the last count found, or 0 if none found but simulation completed
            if sol_counts:
                return sol_counts[-1]
            elif "Simulation completed" in output:
                # If simulation completed but no SOL molecules, count is 0
                return 0
            return None
        
        count1 = extract_final_count(results[0])
        count2 = extract_final_count(results[1])
        
        # With same seed, final counts should be identical
        assert count1 is not None and count2 is not None, f"Could not extract counts. Output1: {results[0][:500]}... Output2: {results[1][:500]}..."
        assert count1 == count2, f"Same seed should produce same results: {count1} != {count2}"
    
    def test_long_simulation_stability(self, temp_dir):
        """Test that longer simulations complete without errors"""
        fragitp_path = TEST_DATA_DIR / "charmm36.ff/mol/sol.itp"
        inp_content = f"""# Stability test
inp_units:nm
box: 3.0 3.0 3.0
temperature: 298.15
fragname: SOL
fragmuex: -10.5
fragconc: 55.0
mcsteps: 1000
nprint: 200
fragitp: {fragitp_path}
"""
        inp_path = Path(temp_dir) / "stability.inp"
        inp_path.write_text(inp_content)
        
        output_prefix = str(Path(temp_dir) / "stability")
        result = subprocess.run(
            [str(GCMC_CPU_PATH), "--inp", str(inp_path),
             "--prefix", output_prefix, "--seed", "42"],
            cwd=temp_dir,
            capture_output=True,
            text=True,
            timeout=30
        )
        
        assert result.returncode == 0
        assert "Simulation completed" in result.stdout
        
        # Check that statistics were printed
        assert "Acceptance:" in result.stdout
        assert "Fragment counts:" in result.stdout


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

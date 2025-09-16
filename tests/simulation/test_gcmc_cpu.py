#!/usr/bin/env python3
"""
Comprehensive test suite for gcmc_cpu executable
Tests various aspects of the GCMC simulation CLI tool
"""

import pytest
import os
import sys
import subprocess
import tempfile
import shutil
import json
import time
from pathlib import Path

# Find the gcmc_cpu executable
GCMC_CPU_PATH = Path(__file__).parent.parent.parent / "build/bin/gcmc_cpu"

class TestGCMCCPU:
    """Test suite for gcmc_cpu executable"""
    
    @pytest.fixture
    def temp_dir(self):
        """Create a temporary directory for test files"""
        temp_path = tempfile.mkdtemp(prefix="gcmc_test_")
        yield temp_path
        shutil.rmtree(temp_path, ignore_errors=True)
    
    @pytest.fixture
    def minimal_files(self, temp_dir):
        """Create minimal test files"""
        # INP file
        inp_path = Path(temp_dir) / "test.inp"
        inp_content = """# Minimal test
top:test.top
pdb:test.pdb
op_top:output.top
op_pdb:output.pdb
box_size:10.0 10.0 10.0
gc_center:5.0 5.0 5.0
cutoff:8.0
mcsteps:10
nprint:1
fragname:water
fragconc:55.0
fragmuex:-5.0
"""
        inp_path.write_text(inp_content)
        
        # PDB file
        pdb_path = Path(temp_dir) / "test.pdb"
        pdb_content = """REMARK Test
CRYST1   10.000   10.000   10.000  90.00  90.00  90.00 P 1           1
ATOM      1  O   WAT     1       5.000   5.000   5.000  1.00  0.00           O
ATOM      2  H1  WAT     1       5.750   5.000   5.000  1.00  0.00           H
ATOM      3  H2  WAT     1       4.750   5.750   5.000  1.00  0.00           H
END
"""
        pdb_path.write_text(pdb_content)
        
        # TOP file
        top_path = Path(temp_dir) / "test.top"
        top_content = """[ system ]
Test System

[ molecules ]
WAT 1
"""
        top_path.write_text(top_content)
        
        return {"inp": inp_path, "pdb": pdb_path, "top": top_path}
    
    def test_executable_exists(self):
        """Test that gcmc_cpu executable exists"""
        assert GCMC_CPU_PATH.exists(), f"gcmc_cpu not found at {GCMC_CPU_PATH}"
        assert os.access(GCMC_CPU_PATH, os.X_OK), f"gcmc_cpu is not executable"
    
    def test_help_option(self):
        """Test --help option"""
        result = subprocess.run(
            [str(GCMC_CPU_PATH), "--help"],
            capture_output=True,
            text=True
        )
        # Help might return 0 or 1 depending on implementation
        assert result.returncode in [0, 1], "Help option should return 0 or 1"
        assert "Usage:" in result.stdout
        assert "--inp" in result.stdout
        assert "--seed" in result.stdout
        assert "--verbose" in result.stdout
    
    def test_missing_input_file(self):
        """Test error handling for missing input file"""
        result = subprocess.run(
            [str(GCMC_CPU_PATH), "--inp", "nonexistent.inp"],
            capture_output=True,
            text=True
        )
        assert result.returncode != 0, "Should fail with nonexistent input"
        assert "ERROR" in result.stdout or "Failed" in result.stdout
    
    def test_minimal_simulation(self, minimal_files, temp_dir):
        """Test minimal GCMC simulation"""
        result = subprocess.run(
            [
                str(GCMC_CPU_PATH),
                "--inp", str(minimal_files["inp"]),
                "--prefix", "test",
                "--seed", "12345",
                "--no-stats"
            ],
            cwd=temp_dir,
            capture_output=True,
            text=True,
            timeout=5
        )
        
        assert result.returncode == 0, f"Simulation failed: {result.stderr}"
        # Check for either verbose or non-verbose output
        assert ("GCMC Simulation" in result.stdout or 
                "Starting GCMC" in result.stdout or
                "Fragment" in result.stdout), "No simulation output found"
        
        # Check output files
        output_file = Path(temp_dir) / "test_final.txt"
        assert output_file.exists(), "Final results file not created"
    
    def test_verbose_output(self, minimal_files, temp_dir):
        """Test verbose output mode"""
        result = subprocess.run(
            [
                str(GCMC_CPU_PATH),
                "--inp", str(minimal_files["inp"]),
                "--verbose",
                "--seed", "42"
            ],
            cwd=temp_dir,
            capture_output=True,
            text=True,
            timeout=5
        )
        
        assert "[INFO]" in result.stdout, "Verbose mode should show INFO messages"
        assert "Temperature:" in result.stdout
        assert "Box size:" in result.stdout
    
    def test_seed_reproducibility(self, minimal_files, temp_dir):
        """Test that same seed produces same results"""
        results = []
        
        for i in range(2):
            result = subprocess.run(
                [
                    str(GCMC_CPU_PATH),
                    "--inp", str(minimal_files["inp"]),
                    "--prefix", f"run{i}",
                    "--seed", "999",
                    "--no-stats"
                ],
                cwd=temp_dir,
                capture_output=True,
                text=True,
                timeout=5
            )
            assert result.returncode == 0
            
            # Read final results
            output_file = Path(temp_dir) / f"run{i}_final.txt"
            results.append(output_file.read_text())
        
        # Results should be identical with same seed
        assert results[0] == results[1], "Same seed should produce identical results"
    
    def test_statistics_collection(self, minimal_files, temp_dir):
        """Test statistics collection"""
        # Modify INP for more steps
        inp_path = minimal_files["inp"]
        content = inp_path.read_text()
        content = content.replace("mcsteps:10", "mcsteps:100")
        content = content.replace("nprint:1", "nprint:10")
        inp_path.write_text(content)
        
        result = subprocess.run(
            [
                str(GCMC_CPU_PATH),
                "--inp", str(inp_path),
                "--stats-interval", "10",
                "--seed", "123"
            ],
            cwd=temp_dir,
            capture_output=True,
            text=True,
            timeout=5
        )
        
        assert result.returncode == 0
        assert "Acceptance:" in result.stdout
        assert "Fragment counts:" in result.stdout
    
    def test_multiple_fragments(self, temp_dir):
        """Test multiple fragment types"""
        # Create INP with multiple fragments
        inp_path = Path(temp_dir) / "multi.inp"
        inp_content = """# Multiple fragments test
top:test.top
pdb:test.pdb
op_top:output.top
op_pdb:output.pdb
box_size:20.0 20.0 20.0
gc_center:10.0 10.0 10.0
cutoff:8.0
mcsteps:50
nprint:10
fragname:water,methanol,ethanol
fragconc:55.0,10.0,5.0
fragmuex:-5.0,-4.5,-4.0
mctime:1.0,0.5,0.3
"""
        inp_path.write_text(inp_content)
        
        # Create dummy PDB and TOP
        pdb_path = Path(temp_dir) / "test.pdb"
        pdb_path.write_text("""ATOM      1  O   WAT     1       0.0   0.0   0.0  1.00  0.00
END""")
        
        top_path = Path(temp_dir) / "test.top"
        top_path.write_text("[ system ]\nTest\n[ molecules ]\nWAT 1")
        
        result = subprocess.run(
            [
                str(GCMC_CPU_PATH),
                "--inp", str(inp_path),
                "--verbose",
                "--seed", "777"
            ],
            cwd=temp_dir,
            capture_output=True,
            text=True,
            timeout=5
        )
        
        assert result.returncode == 0
        assert "water" in result.stdout
        assert "methanol" in result.stdout
        assert "ethanol" in result.stdout
    
    def test_cavity_bias(self, temp_dir):
        """Test cavity bias configuration"""
        inp_path = Path(temp_dir) / "cavity.inp"
        inp_content = """# Cavity bias test
top:test.top
pdb:test.pdb
op_top:output.top
op_pdb:output.pdb
box_size:15.0 15.0 15.0
gc_center:7.5 7.5 7.5
cutoff:8.0
mcsteps:20
nprint:5
use_cavity_bias:yes
cavity_grid_dx:0.5
probe_radius:1.4
fragname:water
fragconc:55.0
fragmuex:-5.0
"""
        inp_path.write_text(inp_content)
        
        # Create dummy files
        pdb_path = Path(temp_dir) / "test.pdb"
        pdb_path.write_text("ATOM      1  O   WAT     1       0.0   0.0   0.0  1.00  0.00\nEND")
        Path(temp_dir).joinpath("test.top").write_text("[ system ]\nTest\n[ molecules ]\nWAT 1")
        
        result = subprocess.run(
            [
                str(GCMC_CPU_PATH),
                "--inp", str(inp_path),
                "--verbose"
            ],
            cwd=temp_dir,
            capture_output=True,
            text=True,
            timeout=5
        )
        
        # Should run without error even if cavity bias not fully implemented
        assert result.returncode == 0
    
    def test_output_prefix(self, minimal_files, temp_dir):
        """Test custom output prefix"""
        custom_prefix = "my_simulation"
        
        result = subprocess.run(
            [
                str(GCMC_CPU_PATH),
                "--inp", str(minimal_files["inp"]),
                "--prefix", custom_prefix,
                "--seed", "555"
            ],
            cwd=temp_dir,
            capture_output=True,
            text=True,
            timeout=5
        )
        
        assert result.returncode == 0
        
        # Check for custom prefix in output files
        expected_file = Path(temp_dir) / f"{custom_prefix}_final.txt"
        assert expected_file.exists(), f"Output file with prefix {custom_prefix} not found"
    
    def test_adaptive_sampling(self, minimal_files, temp_dir):
        """Test adaptive sampling option"""
        result = subprocess.run(
            [
                str(GCMC_CPU_PATH),
                "--inp", str(minimal_files["inp"]),
                "--adaptive",
                "--store-probabilities",
                "--seed", "888"
            ],
            cwd=temp_dir,
            capture_output=True,
            text=True,
            timeout=5
        )
        
        assert result.returncode == 0, "Adaptive sampling should work"
    
    def test_performance_timing(self, minimal_files, temp_dir):
        """Test performance reporting"""
        # Increase steps for meaningful timing
        inp_path = minimal_files["inp"]
        content = inp_path.read_text()
        content = content.replace("mcsteps:10", "mcsteps:1000")
        inp_path.write_text(content)
        
        result = subprocess.run(
            [
                str(GCMC_CPU_PATH),
                "--inp", str(inp_path),
                "--verbose",
                "--seed", "100"
            ],
            cwd=temp_dir,
            capture_output=True,
            text=True,
            timeout=10
        )
        
        assert result.returncode == 0
        assert "Performance:" in result.stdout
        assert "steps/second" in result.stdout or "steps/s" in result.stdout
    
    @pytest.mark.parametrize("seed", [1, 42, 12345, 999999])
    def test_different_seeds(self, minimal_files, temp_dir, seed):
        """Test simulation with different random seeds"""
        result = subprocess.run(
            [
                str(GCMC_CPU_PATH),
                "--inp", str(minimal_files["inp"]),
                "--prefix", f"seed_{seed}",
                "--seed", str(seed),
                "--no-stats"
            ],
            cwd=temp_dir,
            capture_output=True,
            text=True,
            timeout=5
        )
        
        assert result.returncode == 0, f"Failed with seed {seed}"
        
        # Check output file exists
        output_file = Path(temp_dir) / f"seed_{seed}_final.txt"
        assert output_file.exists()

# Additional integration tests
class TestGCMCIntegration:
    """Integration tests for gcmc_cpu"""
    
    @pytest.fixture
    def water_box_setup(self, tmp_path):
        """Setup a water box simulation"""
        # Create a more realistic water box
        pdb_content = """REMARK Water box
CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1           1
"""
        # Add 10 water molecules
        atom_id = 1
        for i in range(10):
            x = 5.0 + i * 2.5
            y = 5.0 + (i % 3) * 2.5
            z = 5.0 + (i % 2) * 2.5
            pdb_content += f"ATOM  {atom_id:5d}  O   WAT {i+1:5d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           O\n"
            atom_id += 1
            pdb_content += f"ATOM  {atom_id:5d}  H1  WAT {i+1:5d}    {x+0.75:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           H\n"
            atom_id += 1
            pdb_content += f"ATOM  {atom_id:5d}  H2  WAT {i+1:5d}    {x:8.3f}{y+0.75:8.3f}{z:8.3f}  1.00  0.00           H\n"
            atom_id += 1
        pdb_content += "END\n"
        
        pdb_path = tmp_path / "waterbox.pdb"
        pdb_path.write_text(pdb_content)
        
        # Create topology
        top_path = tmp_path / "waterbox.top"
        top_path.write_text("""[ system ]
Water Box

[ molecules ]
WAT 10
""")
        
        # Create INP file
        inp_path = tmp_path / "waterbox.inp"
        inp_path.write_text("""# Water box GCMC
top:waterbox.top
pdb:waterbox.pdb
op_top:output.top
op_pdb:output.pdb
box_size:30.0 30.0 30.0
gc_center:15.0 15.0 15.0
cutoff:10.0
mcsteps:500
nprint:50
fragname:water
fragconc:55.0
fragmuex:-5.0
use_cavity_bias:yes
cavity_grid_dx:1.0
""")
        
        return {"inp": inp_path, "pdb": pdb_path, "top": top_path, "dir": tmp_path}
    
    def test_water_box_simulation(self, water_box_setup):
        """Test a realistic water box simulation"""
        result = subprocess.run(
            [
                str(GCMC_CPU_PATH),
                "--inp", str(water_box_setup["inp"]),
                "--prefix", "waterbox",
                "--seed", "2024",
                "--verbose",
                "--stats-interval", "50"
            ],
            cwd=water_box_setup["dir"],
            capture_output=True,
            text=True,
            timeout=30
        )
        
        assert result.returncode == 0, f"Water box simulation failed: {result.stderr}"
        
        # Check final results
        output_file = water_box_setup["dir"] / "waterbox_final.txt"
        assert output_file.exists()
        
        content = output_file.read_text()
        assert "Fragment Statistics:" in content
        assert "water:" in content
        assert "Insert attempts:" in content
        assert "Delete attempts:" in content
    
    def test_long_simulation(self, tmp_path):
        """Test a longer simulation for stability"""
        # Create simple files
        pdb_path = tmp_path / "long.pdb"
        pdb_path.write_text("ATOM      1  O   WAT     1       0.0   0.0   0.0  1.00  0.00\nEND")
        
        top_path = tmp_path / "long.top"
        top_path.write_text("[ system ]\nLong\n[ molecules ]\nWAT 1")
        
        inp_path = tmp_path / "long.inp"
        inp_path.write_text("""# Long simulation test
top:long.top
pdb:long.pdb
op_top:output.top
op_pdb:output.pdb
box_size:50.0 50.0 50.0
gc_center:25.0 25.0 25.0
cutoff:12.0
mcsteps:10000
nprint:1000
fragname:water
fragconc:55.0
fragmuex:-5.0
""")
        
        start_time = time.time()
        
        result = subprocess.run(
            [
                str(GCMC_CPU_PATH),
                "--inp", str(inp_path),
                "--prefix", "long",
                "--seed", "9999",
                "--stats-interval", "1000",
                "--print-freq", "1000"
            ],
            cwd=tmp_path,
            capture_output=True,
            text=True,
            timeout=60
        )
        
        elapsed = time.time() - start_time
        
        assert result.returncode == 0, "Long simulation should complete successfully"
        assert elapsed < 60, f"Simulation took too long: {elapsed:.2f} seconds"
        
        # Check for multiple status updates
        assert result.stdout.count("Step") >= 5, "Should have multiple progress updates"

if __name__ == "__main__":
    # Run tests with pytest
    pytest.main([__file__, "-v"])
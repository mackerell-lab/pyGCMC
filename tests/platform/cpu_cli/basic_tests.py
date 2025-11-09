"""
Basic CLI tests for gcmc_cpu
"""

import pytest
import subprocess
import os


def test_gcmc_cpu_help(gcmc_cpu, temp_dir):
    """Test that gcmc_cpu shows help"""
    result = subprocess.run(
        [gcmc_cpu, "--help"],
        capture_output=True,
        text=True,
        cwd=temp_dir
    )
    
    assert result.returncode != 0  # Help exits with non-zero
    output = result.stdout + result.stderr
    assert "GCMC CPU Driver" in output
    assert "--inp" in output
    assert "--prefix" in output
    assert "--seed" in output


def test_gcmc_cpu_missing_inp(gcmc_cpu, temp_dir):
    """Test that gcmc_cpu fails without input file"""
    result = subprocess.run(
        [gcmc_cpu],
        capture_output=True,
        text=True,
        cwd=temp_dir
    )
    
    assert result.returncode != 0
    output = result.stdout + result.stderr
    assert "Required" in output or "--inp" in output


def test_gcmc_cpu_basic_run(gcmc_cpu, test_data_dir, temp_dir):
    """Test basic gcmc_cpu run with test data"""
    # Use quick version for faster testing
    inp_file = test_data_dir / "gcmc_quick.inp"
    if not inp_file.exists():
        inp_file = test_data_dir / "gcmc.inp"
    
    if not inp_file.exists():
        pytest.skip(f"Test INP file not found: {inp_file}")
    
    # Run with minimal steps for quick test
    cmd = [
        gcmc_cpu,
        "--inp", str(inp_file),
        "--prefix", os.path.join(temp_dir, "test"),
        "--seed", "42",
        "--print-freq", "100",
        "--traj-freq", "1000",
        "--checkpoint-freq", "5000"
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=temp_dir, timeout=20)
    
    # Check that it runs (may not complete successfully due to missing files)
    assert result.returncode == 0 or "Simulation completed" in result.stdout
    
    # Check for expected output patterns
    if result.returncode == 0:
        assert "Simulation completed" in result.stdout

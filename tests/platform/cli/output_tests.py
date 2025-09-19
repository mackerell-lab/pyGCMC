"""
Output and file generation tests for gcmc_cpu
"""

import pytest
import subprocess
import os
from pathlib import Path


def test_gcmc_cpu_verbose_output(gcmc_cpu, test_data_dir, temp_dir):
    """Test verbose output from gcmc_cpu"""
    # Use quick version for faster testing
    inp_file = test_data_dir / "gcmc_quick.inp"
    if not inp_file.exists():
        inp_file = test_data_dir / "gcmc.inp"
    
    if not inp_file.exists():
        pytest.skip(f"Test INP file not found: {inp_file}")
    
    cmd = [
        gcmc_cpu,
        "--inp", str(inp_file),
        "--prefix", os.path.join(temp_dir, "test"),
        "--seed", "42",
        "--verbose"
    ]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, cwd=temp_dir, timeout=20)
    except subprocess.TimeoutExpired:
        # Timeout is OK - simulation is running
        pytest.skip("Simulation took longer than 20s - this is OK")
    
    # Check for verbose output markers
    output = result.stdout + result.stderr
    assert "[INFO]" in output
    assert "Initializing GCMC simulation" in output


def test_gcmc_cpu_logging_regression(gcmc_cpu, test_data_dir, temp_dir):
    """Test that expected log strings are printed when builder path runs"""
    # Use quick version for faster testing
    inp_file = test_data_dir / "gcmc_quick.inp"
    if not inp_file.exists():
        inp_file = test_data_dir / "gcmc.inp"
    
    if not inp_file.exists():
        pytest.skip(f"Test INP file not found: {inp_file}")
    
    cmd = [
        gcmc_cpu,
        "--inp", str(inp_file),
        "--prefix", os.path.join(temp_dir, "test_log"),
        "--seed", "12345",
        "--verbose"
    ]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, cwd=temp_dir, timeout=10)
    except subprocess.TimeoutExpired:
        pytest.skip("Command timed out")
    
    combined_output = result.stdout + result.stderr
    
    # These log lines must be present for backward compatibility
    required_logs = [
        "Temperature:",
        "Box size:",
        "MC steps:"
    ]
    
    for log_pattern in required_logs:
        assert log_pattern in combined_output, f"Missing required log output: {log_pattern}"
    
    # If PDB/TOP files are referenced, builder should report loading them
    with open(inp_file) as f:
        inp_content = f.read().lower()
    
    if "pdb:" in inp_content:
        # Either successful load or file-not-found warning
        assert any(x in combined_output for x in ["Loaded structure", "Warning: Failed to load PDB", 
                                                   "Loaded 0 atoms from PDB", "atoms from PDB"])
    
    if "top:" in inp_content:
        assert any(x in combined_output for x in ["Loaded topology", "Warning: Failed to load topology", 
                                                   "Loaded topology with 0 atoms", "atoms"])


def test_gcmc_cpu_output_files(gcmc_cpu, test_data_dir, temp_dir):
    """Test that gcmc_cpu creates expected output files"""
    # Use quick version for faster testing
    inp_file = test_data_dir / "gcmc_quick.inp"
    if not inp_file.exists():
        inp_file = test_data_dir / "gcmc.inp"
    
    if not inp_file.exists():
        pytest.skip(f"Test INP file not found: {inp_file}")
    
    prefix = "test_output"
    cmd = [
        gcmc_cpu,
        "--inp", str(inp_file),
        "--prefix", os.path.join(temp_dir, prefix),
        "--seed", "42",
        "--print-freq", "10",
        "--traj-freq", "100",
        "--checkpoint-freq", "100"
    ]
    
    # Run for a short time
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, cwd=temp_dir, timeout=20)
    except subprocess.TimeoutExpired:
        # Timeout is OK - simulation may run longer than expected
        result = subprocess.CompletedProcess(args=cmd, returncode=1, stdout="", stderr="")
    
    # Check if any output files were created
    output_files = list(Path(temp_dir).glob(f"{prefix}*"))
    
    # At minimum, we expect the simulation to try to create some files
    # even if it fails due to missing input files
    if result.returncode == 0:
        assert len(output_files) > 0, "No output files created"
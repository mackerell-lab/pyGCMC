"""
Input validation and error handling tests for gcmc_cpu
"""

import pytest
import subprocess
import os
import tempfile
from pathlib import Path


def test_gcmc_cpu_deterministic_seed(gcmc_cpu, test_data_dir, temp_dir):
    """Test that using same seed gives deterministic results"""
    import re

    # Create minimal test input for fast, deterministic runs
    inp_file = os.path.join(temp_dir, "minimal_test.inp")
    pdb_file = os.path.join(temp_dir, "minimal.pdb")
    top_file = os.path.join(temp_dir, "minimal.top")
    itp_file = os.path.join(temp_dir, "water.itp")

    # Create minimal PDB (empty box)
    with open(pdb_file, "w") as f:
        f.write("TITLE     Minimal test system\n")
        f.write("CRYST1   50.000   50.000   50.000  90.00  90.00  90.00 P 1           1\n")
        f.write("END\n")

    # Create minimal TOP
    with open(top_file, "w") as f:
        f.write("[ system ]\nMinimal test\n\n[ molecules ]\n")

    # Create minimal ITP for water
    with open(itp_file, "w") as f:
        f.write("[ moleculetype ]\n")
        f.write("; Name   nrexcl\n")
        f.write("WAT      3\n\n")
        f.write("[ atoms ]\n")
        f.write("1   O     1      WAT      O     1      -0.834   15.999\n")
        f.write("2   H     1      WAT      H1    1       0.417    1.008\n")
        f.write("3   H     1      WAT      H2    1       0.417    1.008\n")

    # Create minimal INP
    with open(inp_file, "w") as f:
        f.write(f"pdb:{pdb_file}\n")
        f.write(f"top:{top_file}\n")
        f.write(f"fragitp:{itp_file}\n")
        f.write(f"op_pdb:{temp_dir}/output.pdb\n")
        f.write(f"op_top:{temp_dir}/output.top\n")
        f.write("box_size:5.0 5.0 5.0\n")
        f.write("cutoff:2.0\n")
        f.write("mcsteps:100\n")  # Very short for speed
        f.write("nprint:50\n")
        f.write("fragname:WAT\n")
        f.write("fragmuex:1.0\n")  # Positive for reasonable acceptance

    # Run twice with same seed and parse numerical results
    parsed_results = []
    for i in range(2):
        cmd = [
            gcmc_cpu,
            "--inp", inp_file,
            "--prefix", os.path.join(temp_dir, f"test_{i}"),
            "--seed", "54321",
            "--verbose"
        ]

        result = subprocess.run(cmd, capture_output=True, text=True, cwd=temp_dir, timeout=5)
        output = result.stdout + result.stderr

        # Parse final statistics - these must be identical for same seed
        stats = {}

        # Parse fragment counts
        count_match = re.search(r"Fragment counts:\s*WAT:\s*(\d+)", output)
        if count_match:
            stats["water_count"] = int(count_match.group(1))

        # Parse acceptance rates
        insert_match = re.search(r"Insert move accept:\s*([\d.]+)%", output)
        if insert_match:
            stats["insert_accept"] = float(insert_match.group(1))

        delete_match = re.search(r"Delete move accept:\s*([\d.]+)%", output)
        if delete_match:
            stats["delete_accept"] = float(delete_match.group(1))

        # Parse energy
        energy_match = re.search(r"Average energy:\s*([\d.-]+)", output)
        if energy_match:
            stats["avg_energy"] = float(energy_match.group(1))

        parsed_results.append(stats)

    # Verify we got meaningful results
    assert len(parsed_results) == 2, "Should have results from both runs"
    assert all(parsed_results), "Both runs should produce statistics"

    # Compare results - they must be IDENTICAL for same seed
    run1, run2 = parsed_results
    for key in run1:
        if key in run2:
            assert run1[key] == run2[key], \
                f"Deterministic failure: {key} differs between runs (Run1: {run1[key]}, Run2: {run2[key]})"


def test_gcmc_cpu_invalid_inp(gcmc_cpu, temp_dir):
    """Test that gcmc_cpu handles invalid input file gracefully"""
    # Create invalid INP file
    bad_inp = os.path.join(temp_dir, "bad.inp")
    with open(bad_inp, "w") as f:
        f.write("INVALID INPUT FILE\n")
    
    cmd = [
        gcmc_cpu,
        "--inp", bad_inp,
        "--prefix", os.path.join(temp_dir, "test")
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=temp_dir)
    
    # Should fail but not crash
    assert result.returncode != 0
    # Should have some error message
    assert len(result.stderr) > 0 or "ERROR" in result.stdout


def test_gcmc_cpu_parameter_validation(gcmc_cpu, test_data_dir, temp_dir):
    """Test parameter validation"""
    import re

    # Create minimal test input with known nprint value
    inp_file = os.path.join(temp_dir, "param_test.inp")
    pdb_file = os.path.join(temp_dir, "empty.pdb")
    top_file = os.path.join(temp_dir, "empty.top")
    itp_file = os.path.join(temp_dir, "water.itp")

    # Create minimal PDB
    with open(pdb_file, "w") as f:
        f.write("TITLE     Minimal\n")
        f.write("CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1           1\n")
        f.write("END\n")

    # Create minimal TOP
    with open(top_file, "w") as f:
        f.write("[ system ]\nTest\n\n[ molecules ]\n")

    # Create minimal ITP
    with open(itp_file, "w") as f:
        f.write("[ moleculetype ]\n")
        f.write("WAT      3\n\n")
        f.write("[ atoms ]\n")
        f.write("1   O     1      WAT      O     1      -0.834   15.999\n")

    # Known nprint value
    EXPECTED_NPRINT = 25

    # Create INP with specific nprint
    with open(inp_file, "w") as f:
        f.write(f"pdb:{pdb_file}\n")
        f.write(f"top:{top_file}\n")
        f.write(f"fragitp:{itp_file}\n")
        f.write(f"op_pdb:{temp_dir}/out.pdb\n")
        f.write(f"op_top:{temp_dir}/out.top\n")
        f.write("box_size:3.0 3.0 3.0\n")
        f.write("cutoff:1.0\n")
        f.write("mcsteps:100\n")
        f.write(f"nprint:{EXPECTED_NPRINT}\n")
        f.write("fragname:WAT\n")
        f.write("fragmuex:1.0\n")

    # Test 1: Invalid seed - should fail
    cmd = [gcmc_cpu, "--inp", inp_file, "--seed", "not_a_number"]
    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        cwd=temp_dir
    )
    assert result.returncode != 0, "Invalid seed should cause failure"

    # Test 2: Negative print-freq should fallback to INP's nprint
    cmd = [gcmc_cpu, "--inp", inp_file, "--print-freq", "-100", "--seed", "12345"]
    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        timeout=5,
        cwd=temp_dir
    )

    # Verify it ran successfully
    assert result.returncode == 0, "Negative print-freq should be handled gracefully"

    # Count how many "=== Step X ===" lines appear (should be mcsteps/nprint)
    step_pattern = r"=== Step (\d+) ==="
    step_matches = re.findall(step_pattern, result.stdout)

    # With 100 steps and nprint=25, we expect steps at: 0, 25, 50, 75, 100 (5 prints)
    # Or at minimum: 0 and 100 (2 prints)
    expected_min_prints = 2  # At least start and end
    expected_max_prints = (100 // EXPECTED_NPRINT) + 2  # Plus some tolerance

    assert len(step_matches) >= expected_min_prints, \
        f"Expected at least {expected_min_prints} step outputs, got {len(step_matches)}"
    assert len(step_matches) <= expected_max_prints, \
        f"Expected at most {expected_max_prints} step outputs, got {len(step_matches)}"

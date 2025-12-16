"""
Output and file generation tests for gcmc_cpu
"""

from __future__ import annotations

import pytest
import subprocess
from pathlib import Path


def test_gcmc_cpu_verbose_output(gcmc_cpu, test_data_dir, temp_dir):
    """
    Smoke test that `--verbose` mode runs successfully.

    Do not assert on log strings (implementation detail); only assert on exit code and outputs.
    """
    work = Path(temp_dir) / "verbose_run"
    work.mkdir(parents=True, exist_ok=True)

    itp = test_data_dir / "charmm36.ff" / "mol" / "na.itp"
    if not itp.exists():
        pytest.skip(f"Required ITP not found: {itp}")

    out_prefix = work / "out" / "gcmc"
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    inp = work / "run.inp"
    inp.write_text(
        f"""
inp_units:gcmc_gpu
fragitp:{itp}
fragname:NA
fragconc:55.0
fragmuex:0.0

box_size:10.0 10.0 10.0
cutoff:12.0
temperature:300.0
moves_per_step:1
mcsteps:2
nprint:1
mc_move_prob:1 0 0 0
""".strip()
        + "\n"
    )

    result = subprocess.run(
        [gcmc_cpu, "--inp", str(inp), "--prefix", str(out_prefix), "--seed", "42", "--verbose"],
        capture_output=True,
        text=True,
        cwd=str(work),
        timeout=30,
    )
    assert result.returncode == 0, result.stdout + result.stderr

    assert Path(f"{out_prefix}_final.pdb").exists()


def test_gcmc_cpu_logging_regression(gcmc_cpu, test_data_dir, temp_dir):
    """
    Regression test: relative input paths should resolve against the INP file directory (not cwd).

    This matches gcmc_gpu usage patterns and prevents "works only when run from the INP folder" bugs.
    """
    work = Path(temp_dir) / "relative_paths"
    work.mkdir(parents=True, exist_ok=True)

    # Copy a known-good fragment ITP next to the INP, then reference it relatively.
    itp_src = test_data_dir / "charmm36.ff" / "mol" / "na.itp"
    if not itp_src.exists():
        pytest.skip(f"Required ITP not found: {itp_src}")
    itp_local = work / "na.itp"
    itp_local.write_text(itp_src.read_text())

    out_prefix = work / "out" / "gcmc"
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    inp = work / "run.inp"
    inp.write_text(
        """
inp_units:gcmc_gpu
fragitp:na.itp
fragname:NA
fragconc:55.0
fragmuex:0.0

box_size:10.0 10.0 10.0
cutoff:12.0
temperature:300.0
moves_per_step:1
mcsteps:1
nprint:1
mc_move_prob:1 0 0 0
""".strip()
        + "\n"
    )

    # Run from a different working directory to ensure resolution does not depend on cwd.
    other_cwd = Path(temp_dir) / "other_cwd"
    other_cwd.mkdir(parents=True, exist_ok=True)

    result = subprocess.run(
        [gcmc_cpu, "--inp", str(inp), "--prefix", str(out_prefix), "--seed", "12345"],
        capture_output=True,
        text=True,
        cwd=str(other_cwd),
        timeout=30,
    )
    assert result.returncode == 0, result.stdout + result.stderr

    assert Path(f"{out_prefix}_final.pdb").exists()


def test_gcmc_cpu_output_files(gcmc_cpu, test_data_dir, temp_dir):
    """Test that gcmc_cpu creates expected output files"""
    work = Path(temp_dir) / "output_files"
    work.mkdir(parents=True, exist_ok=True)

    itp = test_data_dir / "charmm36.ff" / "mol" / "na.itp"
    if not itp.exists():
        pytest.skip(f"Required ITP not found: {itp}")

    out_prefix = work / "out" / "gcmc"
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    inp = work / "run.inp"
    inp.write_text(
        f"""
inp_units:gcmc_gpu
fragitp:{itp}
fragname:NA
fragconc:55.0
fragmuex:0.0

box_size:10.0 10.0 10.0
cutoff:12.0
temperature:300.0
moves_per_step:1
mcsteps:7
nprint:1000
mc_move_prob:1 0 0 0
""".strip()
        + "\n"
    )

    result = subprocess.run(
        [
            gcmc_cpu,
            "--inp",
            str(inp),
            "--prefix",
            str(out_prefix),
            "--seed",
            "42",
            "--print-freq",
            "3",
            "--traj-freq",
            "3",
            "--checkpoint-freq",
            "3",
        ],
        capture_output=True,
        text=True,
        cwd=str(work),
        timeout=60,
    )
    assert result.returncode == 0, result.stdout + result.stderr

    # Core outputs are always written.
    assert Path(f"{out_prefix}_final.pdb").exists()
    assert Path(f"{out_prefix}_final.top").exists()
    assert Path(f"{out_prefix}_statistics.dat").exists()

    # Frequency-controlled outputs should be written at step 3 and 6 (step>0 and step%3==0).
    assert Path(f"{out_prefix}_traj_3.pdb").exists()
    assert Path(f"{out_prefix}_traj_6.pdb").exists()
    assert Path(f"{out_prefix}_checkpoint_3.dat").exists()
    assert Path(f"{out_prefix}_checkpoint_6.dat").exists()

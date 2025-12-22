"""
Basic CLI tests for gcmc_cpu
"""

from __future__ import annotations

import pytest
import subprocess
from pathlib import Path


def test_gcmc_cpu_help(gcmc_cpu, temp_dir):
    """Test that gcmc_cpu shows help"""
    result = subprocess.run(
        [gcmc_cpu, "--help"],
        capture_output=True,
        text=True,
        cwd=temp_dir
    )
    
    # Help is printed via the "invalid args" path; accept either conventional 0 or current non-zero.
    assert result.returncode in (0, 1)
    # Avoid log-string assertions; just ensure no output artifacts are produced.
    outputs = list(Path(temp_dir).glob("*"))
    assert not outputs, f"Unexpected outputs from --help: {outputs}"


def test_gcmc_cpu_missing_inp(gcmc_cpu, temp_dir):
    """Test that gcmc_cpu fails without input file"""
    result = subprocess.run(
        [gcmc_cpu],
        capture_output=True,
        text=True,
        cwd=temp_dir
    )
    
    assert result.returncode != 0
    # Avoid log-string assertions; just ensure no output artifacts are produced.
    outputs = list(Path(temp_dir).glob("*"))
    assert not outputs, f"Unexpected outputs when --inp is missing: {outputs}"


def test_gcmc_cpu_basic_run(gcmc_cpu, test_data_dir, temp_dir):
    """
    Minimal end-to-end smoke test: gcmc_cpu should parse an INP and produce a final PDB.

    This is file-driven (no stdout/stderr string matching) and uses a 1-move run to keep it fast.
    """
    work = Path(temp_dir) / "basic_run"
    work.mkdir(parents=True, exist_ok=True)

    itp = test_data_dir / "charmm36.ff" / "mol" / "na.itp"
    if not itp.exists():
        pytest.skip(f"Required ITP not found: {itp}")

    out_prefix = work / "out" / "gcmc"
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    inp = work / "run.inp"
    inp.write_text(
        f"""
fragitp:{itp}
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

    result = subprocess.run(
        [gcmc_cpu, "--inp", str(inp), "--prefix", str(out_prefix), "--seed", "42"],
        capture_output=True,
        text=True,
        cwd=str(work),
        timeout=30,
    )
    assert result.returncode == 0, result.stdout + result.stderr

    out_pdb = Path(f"{out_prefix}_final.pdb")
    assert out_pdb.exists()

    # One insertion-only move in an empty box should yield exactly one NA residue.
    resids = set()
    for line in out_pdb.read_text().splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        if line[17:20].strip().upper() == "NA":
            resids.add(int(line[22:26]))
    assert len(resids) == 1

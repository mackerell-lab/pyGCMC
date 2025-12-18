"""
Input validation and error handling tests for gcmc_cpu
"""

from __future__ import annotations

import json
import pytest
import subprocess
from pathlib import Path


def test_gcmc_cpu_deterministic_seed(gcmc_cpu, test_data_dir, temp_dir):
    """Same `--seed` should produce deterministic final structures (file-driven; no log parsing)."""
    itp = test_data_dir / "charmm36.ff" / "mol" / "na.itp"
    if not itp.exists():
        pytest.skip(f"Required ITP not found: {itp}")

    def run_once(work: Path, seed: int) -> tuple[tuple[str, str, int, float, float, float], ...]:
        work.mkdir(parents=True, exist_ok=True)
        out_prefix = work / "out" / "gcmc"
        out_prefix.parent.mkdir(parents=True, exist_ok=True)

        inp = work / "run.inp"
        inp.write_text(
            f"""
fragitp:{itp}
fragname:NA
fragconc:55.0
fragmuex:0.0

box_size:30.0 30.0 30.0
cutoff:12.0
temperature:300.0
moves_per_step:1
mcsteps:10
nprint:1000
mc_move_prob:1 0 0 0
""".strip()
            + "\n"
        )

        result = subprocess.run(
            [gcmc_cpu, "--inp", str(inp), "--prefix", str(out_prefix), "--seed", str(seed)],
            capture_output=True,
            text=True,
            cwd=str(work),
            timeout=60,
        )
        assert result.returncode == 0, result.stdout + result.stderr

        out_pdb = Path(f"{out_prefix}_final.pdb")
        assert out_pdb.exists()

        atoms: list[tuple[str, str, int, float, float, float]] = []
        for line in out_pdb.read_text().splitlines():
            if not line.startswith(("ATOM", "HETATM")):
                continue
            atom_name = line[12:16].strip().upper()
            resname = line[17:20].strip().upper()
            resid = int(line[22:26])
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            atoms.append((resname, atom_name, resid, round(x, 3), round(y, 3), round(z, 3)))

        return tuple(sorted(atoms))

    run1 = Path(temp_dir) / "seed_det" / "run1"
    run2 = Path(temp_dir) / "seed_det" / "run2"
    assert run_once(run1, 54321) == run_once(run2, 54321)

    run3 = Path(temp_dir) / "seed_det" / "run3"
    run4 = Path(temp_dir) / "seed_det" / "run4"
    assert run_once(run3, 54321) != run_once(run4, 54322)


def test_gcmc_cpu_invalid_inp(gcmc_cpu, temp_dir):
    """Test that gcmc_cpu handles invalid input file gracefully"""
    # Create invalid INP file
    bad_inp = Path(temp_dir) / "bad.inp"
    bad_inp.write_text("INVALID INPUT FILE\n")
    
    cmd = [
        gcmc_cpu,
        "--inp", str(bad_inp),
        "--prefix", str(Path(temp_dir) / "test")
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=temp_dir)
    
    # Should fail but not crash
    assert result.returncode != 0
    # Should have some error message
    assert len(result.stderr) > 0 or "ERROR" in result.stdout


def test_gcmc_cpu_parameter_validation(gcmc_cpu, test_data_dir, temp_dir):
    """Test CLI arg parsing and print-freq fallback without depending on log text."""
    itp = test_data_dir / "charmm36.ff" / "mol" / "na.itp"
    if not itp.exists():
        pytest.skip(f"Required ITP not found: {itp}")

    work = Path(temp_dir) / "param_validation"
    work.mkdir(parents=True, exist_ok=True)

    out_prefix = work / "out" / "gcmc"
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    expected_nprint = 25

    inp_file = work / "param_test.inp"
    inp_file.write_text(
        f"""
fragitp:{itp}
fragname:NA
fragconc:55.0
fragmuex:0.0

box_size:10.0 10.0 10.0
cutoff:12.0
temperature:300.0
mcsteps:100
nprint:{expected_nprint}
moves_per_step:1
mc_move_prob:1 0 0 0
""".strip()
        + "\n"
    )

    # Test 1: Invalid seed - should fail
    cmd = [gcmc_cpu, "--inp", str(inp_file), "--seed", "not_a_number"]
    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        cwd=str(work),
    )
    assert result.returncode != 0, "Invalid seed should cause failure"

    # Test 2: Negative print-freq should fallback to INP's nprint
    cmd = [
        gcmc_cpu,
        "--inp",
        str(inp_file),
        "--prefix",
        str(out_prefix),
        "--print-freq",
        "-100",
        "--seed",
        "12345",
    ]
    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        timeout=5,
        cwd=str(work),
    )

    # Verify it ran successfully
    assert result.returncode == 0, "Negative print-freq should be handled gracefully"

    # Use the stable statistics DAT file rather than stdout logs.
    stats_path = Path(f"{out_prefix}_statistics.dat")
    assert stats_path.exists()
    steps = []
    for line in stats_path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        steps.append(int(parts[0]))

    # With mcsteps=100 and nprint=25, we should see data at 25,50,75 plus the final line at 100.
    assert steps and steps[-1] == 100
    assert 25 in steps
    assert 50 in steps
    assert 75 in steps
    assert all(s % expected_nprint == 0 for s in steps)


def test_dump_params_reports_unknown_inp_keys(gcmc_cpu, test_data_dir, temp_dir):
    """Unknown/unsupported INP keys must be visible via --dump-params (no stdout parsing)."""
    itp = test_data_dir / "charmm36.ff" / "mol" / "na.itp"
    if not itp.exists():
        pytest.skip(f"Required ITP not found: {itp}")

    work = Path(temp_dir) / "unknown_inp_keys"
    work.mkdir(parents=True, exist_ok=True)

    out_prefix = work / "out" / "gcmc"
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    params_json = work / "out" / "params.json"

    unknown_key = "__definitely_unknown_key__"

    inp = work / "run.inp"
    inp.write_text(
        f"""
{unknown_key}:123
fragitp:{itp}
fragname:NA
fragconc:55.0
fragmuex:0.0

box_size:10.0 10.0 10.0
cutoff:12.0
temperature:300.0
moves_per_step:1
mcsteps:0
nprint:1
mc_move_prob:1 0 0 0
""".strip()
        + "\n"
    )

    result = subprocess.run(
        [gcmc_cpu, "--inp", str(inp), "--prefix", str(out_prefix), "--dump-params", str(params_json)],
        capture_output=True,
        text=True,
        cwd=str(work),
        timeout=30,
    )
    assert result.returncode == 0, result.stdout + result.stderr

    params = json.loads(params_json.read_text())
    unknown = params["basic"]["unknown_inp_keys"]
    assert unknown_key in set(unknown)

"""
CBMC (use_conf_bias) trial-count visibility + strictness regressions.

These tests are file-driven via --dump-params and do not depend on stdout/stderr logs.
"""

from __future__ import annotations

import json
import subprocess
from pathlib import Path

import pytest


def _read_warning_codes(params: dict) -> set[str]:
    warnings = params.get("basic", {}).get("warnings", [])
    return {w.get("code", "") for w in warnings if isinstance(w, dict)}


def _write_inp(path: Path, content: str) -> None:
    path.write_text(content.strip() + "\n")


def test_dump_params_exposes_conf_bias_trials_per_fragment(gcmc_cpu, test_data_dir, temp_dir):
    """When `fragconf` is provided, dump-params must expose per-fragment CBMC trials."""
    na_itp = test_data_dir / "charmm36.ff" / "mol" / "na.itp"
    cl_itp = test_data_dir / "charmm36.ff" / "mol" / "cl.itp"
    if not na_itp.exists() or not cl_itp.exists():
        pytest.skip("Required ITP files not found under tests/data/charmm36.ff/mol/")

    work = Path(temp_dir) / "cbmc_trials_visibility" / "fragconf_list"
    work.mkdir(parents=True, exist_ok=True)

    out_prefix = work / "out" / "gcmc"
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    params_json = work / "out" / "params.json"

    inp = work / "run.inp"
    _write_inp(
        inp,
        f"""
fragitp:{na_itp}
fragitp:{cl_itp}
fragname:NA CL
fragconc:1.0 1.0
fragmuex:0.0 0.0
mctime:1 1

use_conf_bias:yes
fragconf:5 7

box_size:30.0 30.0 30.0
cutoff:12.0
temperature:300.0
moves_per_step:1
mcsteps:0
nprint:1
mc_move_prob:1 0 0 0
""",
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
    assert params.get("bias", {}).get("use_conf_bias") is True
    assert params.get("fragment", {}).get("conf_bias_trials") == [5, 7]
    assert "CBMC_TRIAL_COUNT_DEFAULTED" not in _read_warning_codes(params)
    assert "CBMC_TRIAL_COUNT_LIST_LENGTH_MISMATCH" not in _read_warning_codes(params)


def test_dump_params_uses_num_conf_bias_trial_when_fragconf_missing(gcmc_cpu, test_data_dir, temp_dir):
    """When only `num_conf_bias_trial` is provided, all fragments must inherit it."""
    na_itp = test_data_dir / "charmm36.ff" / "mol" / "na.itp"
    cl_itp = test_data_dir / "charmm36.ff" / "mol" / "cl.itp"
    if not na_itp.exists() or not cl_itp.exists():
        pytest.skip("Required ITP files not found under tests/data/charmm36.ff/mol/")

    work = Path(temp_dir) / "cbmc_trials_visibility" / "global_trials"
    work.mkdir(parents=True, exist_ok=True)

    out_prefix = work / "out" / "gcmc"
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    params_json = work / "out" / "params.json"

    inp = work / "run.inp"
    _write_inp(
        inp,
        f"""
fragitp:{na_itp}
fragitp:{cl_itp}
fragname:NA CL
fragconc:1.0 1.0
fragmuex:0.0 0.0
mctime:1 1

use_conf_bias:yes
num_conf_bias_trial:9

box_size:30.0 30.0 30.0
cutoff:12.0
temperature:300.0
moves_per_step:1
mcsteps:0
nprint:1
mc_move_prob:1 0 0 0
""",
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
    assert params.get("bias", {}).get("use_conf_bias") is True
    assert params.get("fragment", {}).get("conf_bias_trials") == [9, 9]
    assert "CBMC_TRIAL_COUNT_DEFAULTED" not in _read_warning_codes(params)


def test_strict_inp_fails_when_use_conf_bias_has_no_trial_count_source(gcmc_cpu, test_data_dir, temp_dir):
    """
    `--strict-inp` must fail-fast when use_conf_bias is enabled but no trial-count key is provided.
    This prevents silent semantic drift (e.g., unexpected defaults).
    """
    na_itp = test_data_dir / "charmm36.ff" / "mol" / "na.itp"
    if not na_itp.exists():
        pytest.skip(f"Required ITP not found: {na_itp}")

    work = Path(temp_dir) / "cbmc_trials_visibility" / "strict_missing_trials"
    work.mkdir(parents=True, exist_ok=True)

    out_prefix = work / "out" / "gcmc"
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    params_json = work / "out" / "params.json"

    inp = work / "run.inp"
    _write_inp(
        inp,
        f"""
fragitp:{na_itp}
fragname:NA
fragconc:1.0
fragmuex:0.0

use_conf_bias:yes

box_size:30.0 30.0 30.0
cutoff:12.0
temperature:300.0
moves_per_step:1
mcsteps:0
nprint:1
mc_move_prob:1 0 0 0
""",
    )

    result = subprocess.run(
        [
            gcmc_cpu,
            "--inp",
            str(inp),
            "--prefix",
            str(out_prefix),
            "--dump-params",
            str(params_json),
            "--strict-inp",
        ],
        capture_output=True,
        text=True,
        cwd=str(work),
        timeout=30,
    )
    assert result.returncode != 0

    params = json.loads(params_json.read_text())
    assert params["basic"]["unknown_inp_keys"] == []
    assert params["basic"]["ignored_inp_keys"] == []
    assert "CBMC_TRIAL_COUNT_DEFAULTED" in _read_warning_codes(params)

    assert not Path(f"{out_prefix}_final.pdb").exists()
    assert not Path(f"{out_prefix}_final.top").exists()


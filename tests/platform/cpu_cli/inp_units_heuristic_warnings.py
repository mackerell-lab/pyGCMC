"""
Unit/geometry heuristic warnings (structured via --dump-params).

These tests avoid stdout/stderr log matching and instead assert on the JSON schema.
"""

from __future__ import annotations

import json
import subprocess
from pathlib import Path

import pytest


def _read_warning_codes(params: dict) -> set[str]:
    warnings = params.get("basic", {}).get("warnings", [])
    return {w.get("code", "") for w in warnings if isinstance(w, dict)}


def test_nm_mode_suspicious_lengths_emit_unit_warnings(gcmc_cpu, test_data_dir, temp_dir):
    """`inp_units:nm` + legacy-style Å values should emit structured unit warnings."""
    itp = test_data_dir / "charmm36.ff" / "mol" / "na.itp"
    if not itp.exists():
        pytest.skip(f"Required ITP not found: {itp}")

    work = Path(temp_dir) / "inp_units_warnings_nm"
    work.mkdir(parents=True, exist_ok=True)

    out_prefix = work / "out" / "gcmc"
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    params_json = work / "out" / "params.json"

    inp = work / "run.inp"
    inp.write_text(
        f"""
inp_units:nm
fragitp:{itp}
fragname:NA
fragconc:55.0
fragmuex:0.0

box_size:30.0 30.0 30.0
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
    codes = _read_warning_codes(params)
    assert "UNIT_SUSPECT_CUTOFF_TOO_LARGE_FOR_NM" in codes
    assert "UNIT_SUSPECT_BOX_TOO_LARGE_FOR_NM" in codes


def test_strict_inp_fails_on_unit_warnings(gcmc_cpu, test_data_dir, temp_dir):
    """`--strict-inp` must fail fast if heuristic warnings are present (prevents silent unit mistakes)."""
    itp = test_data_dir / "charmm36.ff" / "mol" / "na.itp"
    if not itp.exists():
        pytest.skip(f"Required ITP not found: {itp}")

    work = Path(temp_dir) / "strict_inp_warnings"
    work.mkdir(parents=True, exist_ok=True)

    out_prefix = work / "out" / "gcmc"
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    params_json = work / "out" / "params.json"

    inp = work / "run.inp"
    inp.write_text(
        f"""
inp_units:nm
fragitp:{itp}
fragname:NA
fragconc:55.0
fragmuex:0.0

box_size:30.0 30.0 30.0
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
    assert "UNIT_SUSPECT_CUTOFF_TOO_LARGE_FOR_NM" in _read_warning_codes(params)

    assert not Path(f"{out_prefix}_final.pdb").exists()
    assert not Path(f"{out_prefix}_final.top").exists()


def test_pdb_cryst1_vs_inp_box_size_factor10_mismatch_is_reported(gcmc_cpu, test_data_dir, temp_dir):
    """If PDB CRYST1 differs from INP box_size by ~10x, emit a structured warning."""
    itp = test_data_dir / "charmm36.ff" / "mol" / "na.itp"
    if not itp.exists():
        pytest.skip(f"Required ITP not found: {itp}")

    work = Path(temp_dir) / "cryst1_box_mismatch"
    work.mkdir(parents=True, exist_ok=True)

    pdb = work / "box.pdb"
    pdb.write_text(
        """
CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1           1
ATOM      1  NA   NA  A   1       0.000   0.000   0.000  1.00  0.00          NA
END
""".lstrip()
    )

    out_prefix = work / "out" / "gcmc"
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    params_json = work / "out" / "params.json"

    # Default inp_units:auto => gcmc_gpu (Å). box_size:3.0 => 3Å=0.3nm, while CRYST1=30Å=3nm (~10x mismatch).
    inp = work / "run.inp"
    inp.write_text(
        f"""
pdb:{pdb}
fragitp:{itp}
fragname:NA
fragconc:55.0
fragmuex:0.0

box_size:3.0 3.0 3.0
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
    codes = _read_warning_codes(params)
    assert "UNIT_MISMATCH_PDB_CRYST1_VS_INP_BOX_SIZE" in codes


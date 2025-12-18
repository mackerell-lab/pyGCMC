"""
CLI regression tests for movement step-size parameters.

These are gcmc_gpu/opencl compatibility keys:
- max_translation: length-like (Å in inp_units:auto/gcmc_gpu; nm in inp_units:nm)
- max_rotation: degrees (legacy decks); converted to radians when configuring the engine

We verify this via the stable --dump-params JSON (no stdout/stderr parsing).
"""

from __future__ import annotations

import json
import math
import subprocess
from pathlib import Path

import pytest


def test_max_translation_and_rotation_keys_are_parsed_converted_and_applied(gcmc_cpu, test_data_dir, temp_dir):
    itp = test_data_dir / "charmm36.ff" / "mol" / "na.itp"
    if not itp.exists():
        pytest.skip(f"Required ITP not found: {itp}")

    work = Path(temp_dir) / "move_step_params"
    work.mkdir(parents=True, exist_ok=True)

    out_prefix = work / "out" / "gcmc"
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    params_json = work / "out" / "params.json"

    # Default inp_units:auto => gcmc_gpu-style Å for lengths.
    max_translation_angstrom = 1.5
    max_rotation_deg = 20.0

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
mcsteps:0
nprint:1
mc_move_prob:1 0 0 0

max_translation:{max_translation_angstrom}
max_rotation:{max_rotation_deg}
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
    unknown = set(params["basic"]["unknown_inp_keys"])
    assert "max_translation" not in unknown
    assert "max_rotation" not in unknown

    expected_translation_nm = max_translation_angstrom / 10.0
    assert float(params["mc"]["max_translation_nm"]) == pytest.approx(expected_translation_nm, rel=0, abs=1e-6)
    assert float(params["mc"]["max_rotation_deg"]) == pytest.approx(max_rotation_deg, rel=0, abs=1e-6)

    assert float(params["engine"]["max_translation_nm"]) == pytest.approx(expected_translation_nm, rel=0, abs=1e-6)
    expected_rotation_rad = max_rotation_deg * (math.pi / 180.0)
    assert float(params["engine"]["max_rotation_rad"]) == pytest.approx(expected_rotation_rad, rel=0, abs=1e-12)

"""
INP compatibility tests for gcmc_gpu-style units (dump-params + region conversions).
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from .inp_units_compat_helpers import _run_gcmc_cpu, _write_inp


def test_inp_units_gcmc_gpu_converts_grid_dx_and_cutoffs_and_target_volume(
    gcmc_cpu, test_data_dir, temp_dir
):
    """
    Verify unit conversions for less-visible keys via a stable JSON dump (no log matching):
    - grid_dx (Å -> nm)
    - energy_cutoff_frag (Å -> nm) and mirrored cutoff
    - target_volume (Å^3 -> nm^3)
    """
    work = Path(temp_dir) / "units_dump_params"
    work.mkdir(parents=True, exist_ok=True)

    itp = test_data_dir / "charmm36.ff" / "mol" / "na.itp"
    assert itp.exists()

    out_prefix = work / "out" / "gcmc"
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    params_json = work / "out" / "params.json"

    inp = work / "test.inp"
    _write_inp(
        inp,
        f"""
inp_units:gcmc_gpu
fragitp:{itp}
fragname:NA
fragconc:1.0
fragmuex:0.0

box_size:10.0 10.0 10.0
cutoff:12.0
grid_dx:1.0
energy_cutoff_frag:12.0
target_volume:1000.0

temperature:300.0
moves_per_step:1
mcsteps:1
nprint:1
mc_move_prob:1 0 0 0
""",
    )

    result = _run_gcmc_cpu(
        gcmc_cpu,
        workdir=work,
        inp=inp,
        out_prefix=out_prefix,
        extra_args=["--dump-params", str(params_json)],
    )
    assert result.returncode == 0, result.stdout + result.stderr

    assert params_json.exists(), "Expected --dump-params to create params.json"
    params = json.loads(params_json.read_text())

    # Internal units are nm / kJ/mol.
    assert float(params["space"]["grid_spacing_nm"]) == pytest.approx(0.1, abs=1e-6)
    assert float(params["space"]["cutoff_nm"]) == pytest.approx(1.2, abs=1e-6)
    assert float(params["energy"]["fragment_cutoff_nm"]) == pytest.approx(1.2, abs=1e-6)
    assert float(params["space"]["target_volume_nm3"]) == pytest.approx(1.0, abs=1e-6)


@pytest.mark.parametrize(
    ("region_spec_angstrom", "expected_volume_nm3"),
    [
        # Use relatively large regions to avoid occasional rejection due to hard region constraint.
        ("sphere 15.0 15.0 15.0 15.0", (4.0 / 3.0) * 3.141592653589793 * (1.5**3)),  # r=15Å=1.5nm
        ("box 2.0 2.0 2.0 28.0 28.0 28.0", 2.6**3),  # 0.2..2.8 nm inside a 3 nm box
        ("cylinder 15.0 15.0 15.0 15.0 20.0 z", 3.141592653589793 * (1.5**2) * 2.0),
    ],
)
def test_inp_units_gcmc_gpu_gcmc_region_numeric_conversion_affects_volume(
    gcmc_cpu, test_data_dir, temp_dir, region_spec_angstrom, expected_volume_nm3
):
    """
    gcmc_region numeric values are in Å for gcmc_gpu inputs, but RegionConstraint expects nm.
    Validate end-to-end by asserting the effective acceptance volume (vEff) equals the converted region volume.
    """
    work = Path(temp_dir) / "units_region" / region_spec_angstrom.split()[0]
    work.mkdir(parents=True, exist_ok=True)

    itp = test_data_dir / "charmm36.ff" / "mol" / "sol.itp"
    assert itp.exists()

    out_prefix = work / "out" / "gcmc"
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    accept_log = work / "out" / "acceptance.jsonl"

    inp = work / "test.inp"
    _write_inp(
        inp,
        f"""
inp_units:gcmc_gpu
random_seed:123
fragitp:{itp}
fragname:SOL
fragmuex:0.0

box_size:30.0 30.0 30.0
cutoff:12.0
gcmc_region:{region_spec_angstrom}

temperature:300.0
moves_per_step:1
mcsteps:1
nprint:1
mc_move_prob:1 0 0 0
""",
    )

    result = _run_gcmc_cpu(
        gcmc_cpu,
        workdir=work,
        inp=inp,
        out_prefix=out_prefix,
        extra_args=["--dump-accept", str(accept_log)],
        timeout=30,
    )
    assert result.returncode == 0, result.stdout + result.stderr

    records = [
        json.loads(line)
        for line in accept_log.read_text().splitlines()
        if line.strip()
    ]
    assert records, "Acceptance log unexpectedly empty"

    ins = next(
        (r for r in records if r.get("move") == "insertion" and str(r.get("species", "")).upper() == "SOL"),
        None,
    )
    assert ins is not None, f"Expected a SOL insertion record, got: {records[:3]}"
    assert float(ins["cavityFraction"]) == pytest.approx(1.0, abs=1e-12)
    assert float(ins["vEff"]) == pytest.approx(expected_volume_nm3, rel=1e-6, abs=1e-6)


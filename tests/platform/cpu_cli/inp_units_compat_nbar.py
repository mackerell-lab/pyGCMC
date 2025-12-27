"""
INP compatibility tests for gcmc_gpu-style units (nbar/activity modes).
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from .inp_units_compat_helpers import _first_accept_record, _run_gcmc_cpu, _write_inp


def test_nbar_volume_based_mode_sets_activity_from_concentration(gcmc_cpu, test_data_dir, temp_dir):
    """
    Default (volume-based) mode: activity should be conc(M) * 0.602214... (molecules/nm^3).
    Use a 1 nm^3 box so activity is directly observable from acceptance log 'z'.
    """
    work = Path(temp_dir) / "nbar_volume_based"
    work.mkdir(parents=True, exist_ok=True)

    itp = test_data_dir / "charmm36.ff" / "mol" / "na.itp"
    assert itp.exists()

    out_prefix = work / "out" / "gcmc"
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    accept_log = work / "out" / "acceptance.jsonl"

    inp = work / "test.inp"
    _write_inp(
        inp,
        f"""
random_seed:7
fragitp:{itp}
fragname:NA
fragconc:1.0
fragmuex:0.0

box_size:10.0 10.0 10.0
cutoff:4.0
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
    )
    assert result.returncode == 0, result.stdout + result.stderr

    rec = _first_accept_record(accept_log, move="insertion", species="NA")
    assert float(rec["z"]) == pytest.approx(0.6022, rel=1e-4, abs=1e-4)


def test_nbar_const_water_nbar_scales_all_fragments_by_fragconc(gcmc_cpu, test_data_dir, temp_dir):
    """
    const_water_nbar mode (gcmc_gpu semantics): for any fragment i,
      nbar_i = const_water_nbar / water_density * fragconc_i
      activity_i = nbar_i / V * exp(beta*muex_i)
    This must apply to non-water fragments too (not just SOL/WAT).
    """
    work = Path(temp_dir) / "nbar_const_water_all_species"
    work.mkdir(parents=True, exist_ok=True)

    itp = test_data_dir / "charmm36.ff" / "mol" / "na.itp"
    assert itp.exists()

    out_prefix = work / "out" / "gcmc"
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    accept_log = work / "out" / "acceptance.jsonl"

    inp = work / "test.inp"
    _write_inp(
        inp,
        f"""
random_seed:11
use_const_water_nbar:55
fragitp:{itp}
fragname:NA
fragconc:1.0
fragmuex:0.0

box_size:10.0 10.0 10.0
cutoff:4.0
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
    )
    assert result.returncode == 0, result.stdout + result.stderr

    rec = _first_accept_record(accept_log, move="insertion", species="NA")
    # V = 1 nm^3, water_density defaults to 55 M, muex=0 => activity = 55/55 * 1 / 1 = 1.0
    assert float(rec["z"]) == pytest.approx(1.0, abs=1e-12)


def test_nbar_const_water_nbar_scales_multiple_fragments_in_one_run(gcmc_cpu, test_data_dir, temp_dir):
    """
    Multi-fragment end-to-end check (same run):
    const_water_nbar mode must scale activity 'z' by fragconc for *each* fragment.

    Use a 1 nm^3 empty box so z is directly observable from the acceptance JSONL.
    """
    work = Path(temp_dir) / "nbar_const_water_multi_species"
    work.mkdir(parents=True, exist_ok=True)

    na_itp = test_data_dir / "charmm36.ff" / "mol" / "na.itp"
    cl_itp = test_data_dir / "charmm36.ff" / "mol" / "cl.itp"
    assert na_itp.exists()
    assert cl_itp.exists()

    out_prefix = work / "out" / "gcmc"
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    accept_log = work / "out" / "acceptance.jsonl"

    inp = work / "test.inp"
    _write_inp(
        inp,
        f"""
random_seed:77
use_const_water_nbar:55
fragitp:{na_itp}
fragitp:{cl_itp}
fragname:NA CL
fragconc:0.10 0.20
fragmuex:0.0 0.0
mctime:1.0 1.0

box_size:10.0 10.0 10.0
cutoff:4.0
temperature:300.0
moves_per_step:1
mcsteps:40
nprint:1000
mc_move_prob:1 0 0 0
""",
    )

    result = _run_gcmc_cpu(
        gcmc_cpu,
        workdir=work,
        inp=inp,
        out_prefix=out_prefix,
        extra_args=["--dump-accept", str(accept_log)],
        timeout=60,
    )
    assert result.returncode == 0, result.stdout + result.stderr

    records = [json.loads(line) for line in accept_log.read_text().splitlines() if line.strip()]
    assert records

    def first_z(species: str) -> float:
        want = species.strip().upper()
        for r in records:
            if r.get("move") != "insertion":
                continue
            if str(r.get("species", "")).strip().upper() != want:
                continue
            return float(r["z"])
        raise AssertionError(f"No insertion for {species} in acceptance log; first records: {records[:5]}")

    z_na = first_z("NA")
    z_cl = first_z("CL")

    # With V=1 nm^3 and const_water_nbar==water_density (55 M), z == fragconc.
    # Implementation stores float-like values; keep tolerances loose but meaningful.
    assert z_na == pytest.approx(0.10, rel=1e-6, abs=1e-8)
    assert z_cl == pytest.approx(0.20, rel=1e-6, abs=1e-8)
    assert (z_cl / z_na) == pytest.approx(2.0, rel=1e-6, abs=1e-8)


def test_nbar_number_water_nbar_updates_activity_after_first_insertion(gcmc_cpu, test_data_dir, temp_dir):
    """
    number_water_nbar mode: after the first accepted water insertion, waterCount==1 so
    activity for water should update to z = 1/V (with V=1 nm^3 here).
    """
    work = Path(temp_dir) / "nbar_number_water"
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
random_seed:123
use_number_water_nbar:yes
fragitp:{itp}
fragname:SOL
fragconc:55.0
fragmuex:0.0

box_size:10.0 10.0 10.0
cutoff:4.0
temperature:300.0
moves_per_step:1
mcsteps:2
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

    records = [json.loads(line) for line in accept_log.read_text().splitlines() if line.strip()]
    sol_ins = [r for r in records if r.get("move") == "insertion" and str(r.get("species", "")).upper() == "SOL"]
    assert len(sol_ins) >= 2, f"Expected >=2 SOL insertions, got {len(sol_ins)}"

    # First move: before any waters exist, leave base activity from concentration.
    assert float(sol_ins[0]["z"]) == pytest.approx(55.0 * 0.6022, rel=1e-4, abs=1e-4)

    # Second move: after 1 accepted insertion, number-water nbar sets activity to N/V (N=1, V=1 nm^3).
    assert float(sol_ins[1]["z"]) == pytest.approx(1.0, abs=1e-12)

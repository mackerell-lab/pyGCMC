"""
INP compatibility tests for gcmc_gpu-style units and legacy keys (box + μex basics).

These are end-to-end (CLI) tests that exercise InpParserGCMC via gcmc_cpu.
"""

from __future__ import annotations

import json
import subprocess
from pathlib import Path

import pytest

from .inp_units_compat_helpers import _read_cryst1_box_angstrom, _run_gcmc_cpu, _write_inp


def test_inp_units_gcmc_gpu_box_roundtrip_cryst1(gcmc_cpu, test_data_dir, temp_dir):
    """If INP is in Å, internal nm should roundtrip back to Å in output CRYST1."""
    work = Path(temp_dir) / "units_roundtrip"
    work.mkdir(parents=True, exist_ok=True)

    itp = test_data_dir / "charmm36.ff" / "mol" / "sol.itp"
    assert itp.exists()

    out_prefix = work / "out" / "gcmc"
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    inp = work / "test.inp"
    _write_inp(
        inp,
        f"""
inp_units:gcmc_gpu
fragitp:{itp}
box_size:30.0 31.0 32.0
cutoff:12.0
temperature:300.0
moves_per_step:1
mcsteps:5
nprint:1
fragname:SOL
fragmuex:0.0
attempt_prob_ins:1.0
attempt_prob_del:0.0
attempt_prob_trn:0.0
attempt_prob_rot:0.0
""",
    )

    result = subprocess.run(
        [gcmc_cpu, "--inp", str(inp), "--prefix", str(out_prefix), "--verbose"],
        cwd=str(work),
        capture_output=True,
        text=True,
        timeout=20,
    )
    assert result.returncode == 0, result.stdout + result.stderr

    out_pdb = Path(f"{out_prefix}_final.pdb")
    assert out_pdb.exists()

    lx, ly, lz = _read_cryst1_box_angstrom(out_pdb)
    assert lx == pytest.approx(30.0, abs=1e-3)
    assert ly == pytest.approx(31.0, abs=1e-3)
    assert lz == pytest.approx(32.0, abs=1e-3)


def test_inp_units_gcmc_gpu_fragmuex_kcal_to_kj_log(gcmc_cpu, test_data_dir, temp_dir):
    """fragmuex is kcal/mol in legacy mode and should be converted to kJ/mol internally."""
    work = Path(temp_dir) / "muex_units"
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
fragitp:{itp}
box_size:30.0 30.0 30.0
cutoff:12.0
temperature:300.0
moves_per_step:1
mcsteps:1
nprint:1
fragname:SOL
fragmuex:-1.0
attempt_prob_ins:1.0
attempt_prob_del:0.0
attempt_prob_trn:0.0
attempt_prob_rot:0.0
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

    assert accept_log.exists()
    records = [
        json.loads(line)
        for line in accept_log.read_text().splitlines()
        if line.strip()
    ]
    assert records, "Acceptance log unexpectedly empty"

    mu_values = [
        float(r["mu"])
        for r in records
        if str(r.get("species", "")).strip().upper() == "SOL"
    ]
    assert mu_values, f"No SOL moves found in acceptance log: {records[:3]}"
    assert mu_values[0] == pytest.approx(-4.184, rel=1e-3, abs=1e-3)


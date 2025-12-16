"""
INP compatibility tests for gcmc_gpu-style outputs (active_*.dat + muex_*.dat).
"""

from __future__ import annotations

from pathlib import Path

import pytest

from .inp_units_compat_helpers import _run_gcmc_cpu, _write_inp


def test_gcmc_gpu_active_and_muex_files_are_written(gcmc_cpu, test_data_dir, temp_dir):
    """
    gcmc_gpu-style outputs should be produced without relying on stdout:
    - active_<frag>.dat
    - muex_<frag>.dat (kcal/mol)
    Files should be written next to the chosen --prefix path.
    """
    work = Path(temp_dir) / "outputs_active_muex"
    work.mkdir(parents=True, exist_ok=True)

    itp = test_data_dir / "charmm36.ff" / "mol" / "sol.itp"
    assert itp.exists()

    out_prefix = work / "out" / "gcmc"
    out_dir = out_prefix.parent
    out_dir.mkdir(parents=True, exist_ok=True)

    inp = work / "test.inp"
    _write_inp(
        inp,
        f"""
inp_units:gcmc_gpu
random_seed:321
fragitp:{itp}
fragname:SOL
fragconc:55.0
fragmuex:-1.0

box_size:10.0 10.0 10.0
cutoff:12.0
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
        timeout=30,
    )
    assert result.returncode == 0, result.stdout + result.stderr

    active = out_dir / "active_SOL.dat"
    muex = out_dir / "muex_SOL.dat"
    assert active.exists()
    assert muex.exists()

    active_lines = [line.strip() for line in active.read_text().splitlines() if line.strip()]
    muex_lines = [line.strip() for line in muex.read_text().splitlines() if line.strip()]
    assert active_lines, "active_<frag>.dat unexpectedly empty"
    assert muex_lines, "muex_<frag>.dat unexpectedly empty"

    active_last = int(active_lines[-1])
    assert active_last >= 0

    out_pdb = Path(f"{out_prefix}_final.pdb")
    assert out_pdb.exists()
    sol_resids = set()
    for line in out_pdb.read_text().splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        resname = line[17:20].strip().upper()
        if resname != "SOL":
            continue
        sol_resids.add(int(line[22:26]))
    assert active_last == len(sol_resids)

    muex_last = float(muex_lines[-1])
    assert muex_last == pytest.approx(-1.0, abs=0.02)


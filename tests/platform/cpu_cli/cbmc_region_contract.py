"""
CBMC + region constraint contract tests (no stdout/stderr log parsing).
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

from .inp_units_compat_helpers import _run_gcmc_cpu, _write_inp


def _pdb_residue_centers_angstrom(pdb_path: Path, resname: str) -> list[tuple[float, float, float, int]]:
    want = resname.strip().upper()
    accum: dict[int, tuple[float, float, float, int]] = {}
    for line in pdb_path.read_text().splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        if line[17:20].strip().upper() != want:
            continue
        try:
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            resid = int(line[22:26])
        except ValueError as exc:
            raise AssertionError(f"Failed to parse PDB coordinates: {line}") from exc
        if resid in accum:
            sx, sy, sz, count = accum[resid]
            accum[resid] = (sx + x, sy + y, sz + z, count + 1)
        else:
            accum[resid] = (x, y, z, 1)

    centers: list[tuple[float, float, float, int]] = []
    for resid, (sx, sy, sz, count) in accum.items():
        if count <= 0:
            continue
        centers.append((sx / count, sy / count, sz / count, resid))
    return centers


def _assert_points_inside_sphere(points: list[tuple[float, float, float, int]],
                                 center: tuple[float, float, float],
                                 radius: float) -> None:
    cx, cy, cz = center
    r2 = radius * radius
    for x, y, z, resid in points:
        dx = x - cx
        dy = y - cy
        dz = z - cz
        dist2 = dx * dx + dy * dy + dz * dz
        assert dist2 <= r2 + 1e-6, (
            f"Resid {resid} center at ({x:.3f},{y:.3f},{z:.3f}) outside region sphere "
            f"(center={center}, r={radius})"
        )


def test_cbmc_insertion_respects_gcmc_region(gcmc_cpu, test_data_dir, temp_dir):
    """
    CBMC insertion must still enforce gcmc_region at the atom level.
    """
    work = Path(temp_dir) / "cbmc_region_contract"
    work.mkdir(parents=True, exist_ok=True)

    par = Path(test_data_dir) / "charmm36.ff" / "ffnonbonded.itp"
    fragitp = Path(test_data_dir) / "charmm36.ff" / "mol" / "sol.itp"
    monomerdir = Path(test_data_dir) / "charmm36.ff" / "mol"
    pdb = Path(test_data_dir) / "gcmc_examples" / "empty.pdb"
    top = Path(test_data_dir) / "gcmc_examples" / "empty.top"
    for path in (par, fragitp, monomerdir, pdb, top):
        assert path.exists(), f"Missing test data: {path}"

    center = (10.0, 10.0, 10.0)  # Å
    radius = 3.0  # Å

    out_prefix = work / "out" / "gcmc"
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    accept_log = work / "out" / "acceptance.jsonl"

    inp = work / "run.inp"
    _write_inp(
        inp,
        f"""
version:gcmc_2.0
par:{par}
fragitp:{fragitp}
monomerdir:{monomerdir}
top:{top}
pdb:{pdb}
protitp:{top}

fragname:SOL
fragconc:55.0
fragmuex:0.0
fragconf:5
use_conf_bias:yes

box_size:20.0 20.0 20.0
gcmc_region:sphere {center[0]} {center[1]} {center[2]} {radius}
cutoff:10.0
grid_dx:1.0
temperature:300.0
mcsteps:200
nprint:200
moves_per_step:1
mc_move_prob:1.0 0.0 0.0 0.0
random_seed:4242
""",
    )

    result = _run_gcmc_cpu(
        gcmc_cpu,
        workdir=work,
        inp=inp,
        out_prefix=out_prefix,
        extra_args=["--dump-accept", str(accept_log), "--store-probabilities"],
        timeout=120,
    )
    assert result.returncode == 0, result.stdout + result.stderr

    assert accept_log.exists()
    records = [json.loads(line) for line in accept_log.read_text().splitlines() if line.strip()]
    insertions = [
        r for r in records
        if str(r.get("move", "")).lower() == "insertion"
        and str(r.get("species", "")).upper() == "SOL"
        and bool(r.get("accepted"))
        and int(r.get("cbmcTrials", 1)) > 1
    ]
    assert insertions, "No accepted CBMC insertion records found"

    final_pdb = Path(f"{out_prefix}_final.pdb")
    assert final_pdb.exists()
    centers = _pdb_residue_centers_angstrom(final_pdb, resname="SOL")
    assert centers, "No SOL residues found in final PDB"

    _assert_points_inside_sphere(centers, center=center, radius=radius)

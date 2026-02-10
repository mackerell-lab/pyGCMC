"""
Shared helpers for gcmc_gpu-style INP compatibility tests.
"""

from __future__ import annotations

import json
import os
import subprocess
from pathlib import Path


def _scaled_timeout_seconds(timeout: int) -> int:
    """Scale subprocess timeout under high xdist parallelism to reduce false timeouts."""
    raw_workers = os.getenv("PYTEST_XDIST_WORKER_COUNT", "1")
    try:
        worker_count = max(1, int(raw_workers))
    except ValueError:
        worker_count = 1

    if worker_count >= 8:
        return max(timeout, timeout * 2)
    if worker_count >= 6:
        return max(timeout, int(timeout * 3 / 2))
    return timeout


def _write_inp(path: Path, content: str) -> None:
    path.write_text(content.strip() + "\n")


def _read_cryst1_box_angstrom(pdb_path: Path) -> tuple[float, float, float]:
    for line in pdb_path.read_text().splitlines():
        if line.startswith("CRYST1"):
            parts = line.split()
            assert len(parts) >= 4, f"Unexpected CRYST1 format: {line}"
            return float(parts[1]), float(parts[2]), float(parts[3])
    raise AssertionError(f"CRYST1 not found in {pdb_path}")


def _count_residues_by_resname(pdb_path: Path, resname: str) -> int:
    want = resname.strip().upper()
    resids: set[int] = set()
    for line in pdb_path.read_text().splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        if line[17:20].strip().upper() != want:
            continue
        resids.add(int(line[22:26]))
    return len(resids)


def _run_gcmc_cpu(
    gcmc_cpu: str,
    *,
    workdir: Path,
    inp: Path,
    out_prefix: Path,
    extra_args: list[str] | None = None,
    timeout: int = 20,
) -> subprocess.CompletedProcess[str]:
    args = [gcmc_cpu, "--inp", str(inp), "--prefix", str(out_prefix)]
    if extra_args:
        args.extend(extra_args)
    return subprocess.run(
        args,
        cwd=str(workdir),
        capture_output=True,
        text=True,
        timeout=_scaled_timeout_seconds(timeout),
    )


def _first_accept_record(path: Path, *, move: str, species: str) -> dict:
    records = [json.loads(line) for line in path.read_text().splitlines() if line.strip()]
    want_move = move.strip().lower()
    want_species = species.strip().upper()
    for r in records:
        if str(r.get("move", "")).strip().lower() != want_move:
            continue
        if str(r.get("species", "")).strip().upper() != want_species:
            continue
        return r
    raise AssertionError(f"No {move}/{species} record found in {path}; first records: {records[:3]}")

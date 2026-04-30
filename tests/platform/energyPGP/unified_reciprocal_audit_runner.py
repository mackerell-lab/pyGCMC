"""Smoke test entry for the unified PGP/PME reciprocal audit runner."""

from __future__ import annotations

import json
import runpy
import sys
from pathlib import Path

import pytest


def _find_workspace_root() -> Path:
    current = Path(__file__).resolve()
    for parent in current.parents:
        if (parent / "tmp" / "unified_pgp_pme_reciprocal_audit.py").is_file():
            return parent
    raise FileNotFoundError("Cannot locate workspace root containing tmp/unified_pgp_pme_reciprocal_audit.py")


@pytest.mark.parametrize(
    ("mesh", "expected_mismatch_range"),
    [
        (32, (5.0, 15.0)),   # coarse-mesh baseline monitor
        (64, (0.1, 2.0)),
        (128, (0.0, 0.8)),
    ],
)
def test_unified_pgp_pme_reciprocal_audit_runner(tmp_path, mesh, expected_mismatch_range):
    workspace_root = _find_workspace_root()
    script = workspace_root / "tmp" / "unified_pgp_pme_reciprocal_audit.py"
    output_dir = tmp_path / f"unified_audit_m{mesh}"
    trend_path = output_dir / "trend.json"

    original_argv = sys.argv[:]
    try:
        sys.argv = [
            str(script),
            "--meshes",
            str(mesh),
            "--primary-displacement",
            "0.2,-0.1,0.15",
            "--displacement-sweep",
            "0.2,-0.1,0.15;0.1,0.0,0.0",
            "--output-dir",
            str(output_dir),
            "--trend-file",
            str(trend_path),
        ]
        runpy.run_path(str(script), run_name="__main__")
    finally:
        sys.argv = original_argv

    summary_path = output_dir / "summary.json"
    report_path = output_dir / "report.md"
    assert summary_path.exists(), "Unified audit runner must emit summary.json"
    assert report_path.exists(), "Unified audit runner must emit report.md"
    assert trend_path.exists(), "Unified audit runner must emit trend file"

    summary = json.loads(summary_path.read_text())
    assert summary["checks"]["cross_identity_pass"] is True
    assert summary["checks"]["decomposition_pass"] is True
    assert summary["settings"]["meshes"] == [mesh]
    assert len(summary["cross_identity"]) == 1
    assert len(summary["delta_primary"]) == 1
    assert len(summary["delta_sweep"]) == 2

    primary = summary["delta_primary"][0]
    mismatch = float(primary["mismatch_percent_vs_cross"])
    lower, upper = expected_mismatch_range
    assert lower <= mismatch <= upper, (
        f"mesh {mesh} baseline drift: mismatch={mismatch:.6f}% not in [{lower}, {upper}]"
    )

    trend = json.loads(trend_path.read_text())
    assert "history" in trend and isinstance(trend["history"], list) and len(trend["history"]) >= 1
    last_entry = trend["history"][-1]
    assert str(mesh) in last_entry["mismatch_primary_by_mesh"]

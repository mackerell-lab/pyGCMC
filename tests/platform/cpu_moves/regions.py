"""
Test GCMC region constraint functionality
"""

import pytest
from pathlib import Path


TEST_DATA_DIR = Path(__file__).resolve().parents[2] / "data"
GCMC_CPU_PATH = Path(__file__).resolve().parents[3] / "build" / "bin" / "gcmc_cpu"


def _pdb_atom_positions_angstrom(pdb_path: Path, *, resname: str) -> list[tuple[float, float, float]]:
    positions: list[tuple[float, float, float]] = []
    for line in pdb_path.read_text().splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        if line[17:20].strip().upper() != resname.upper():
            continue
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        positions.append((x, y, z))
    return positions


def _assert_all_in_sphere(
    positions: list[tuple[float, float, float]],
    *,
    center: tuple[float, float, float],
    radius: float,
) -> None:
    cx, cy, cz = center
    r2 = radius * radius
    for x, y, z in positions:
        dx = x - cx
        dy = y - cy
        dz = z - cz
        assert (dx * dx + dy * dy + dz * dz) <= r2 + 1e-6


def _assert_all_in_box(
    positions: list[tuple[float, float, float]],
    *,
    min_corner: tuple[float, float, float],
    max_corner: tuple[float, float, float],
) -> None:
    x1, y1, z1 = min_corner
    x2, y2, z2 = max_corner
    for x, y, z in positions:
        assert x1 - 1e-6 <= x <= x2 + 1e-6
        assert y1 - 1e-6 <= y <= y2 + 1e-6
        assert z1 - 1e-6 <= z <= z2 + 1e-6


def _assert_all_in_cylinder_z(
    positions: list[tuple[float, float, float]],
    *,
    center: tuple[float, float, float],
    radius: float,
    height: float,
) -> None:
    cx, cy, cz = center
    r2 = radius * radius
    z_min = cz - 0.5 * height
    z_max = cz + 0.5 * height
    for x, y, z in positions:
        assert z_min - 1e-6 <= z <= z_max + 1e-6
        dx = x - cx
        dy = y - cy
        assert (dx * dx + dy * dy) <= r2 + 1e-6


def test_gcmc_region_sphere(tmp_path):
    """Sphere region must constrain insertion positions (file-driven, no stdout parsing)."""
    fragitp_path = TEST_DATA_DIR / "charmm36.ff/mol/na.itp"

    # Create a simple inp file with gcmc_region
    inp_content = f"""
version: gcmc_2.0
num_threads: 1
temperature: 300.0
mcsteps: 20
nprint: 20
fragname: NA
fragmuex: 100.0
fragconc: 1.0
gcmc_region: sphere 25.0 25.0 25.0 10.0
box: 50.0 50.0 50.0
cutoff: 12.0
mc_move_prob: 1 0 0 0
fragitp: {fragitp_path}
"""

    inp_file = tmp_path / "test_sphere.inp"
    inp_file.write_text(inp_content)

    # Import and run gcmc_cpu
    import subprocess

    # Run gcmc_cpu with the inp file
    gcmc_executable = GCMC_CPU_PATH
    if gcmc_executable.exists():
        out_prefix = tmp_path / "out" / "gcmc"
        out_prefix.parent.mkdir(parents=True, exist_ok=True)
        result = subprocess.run(
            [str(gcmc_executable), "--inp", str(inp_file), "--prefix", str(out_prefix), "--seed", "42"],
            cwd=tmp_path,
            capture_output=True,
            text=True
        )

        # Check that the simulation ran successfully
        assert result.returncode == 0, f"gcmc_cpu failed with error: {result.stderr}"

        final_pdb = Path(f"{out_prefix}_final.pdb")
        assert final_pdb.exists()

        positions = _pdb_atom_positions_angstrom(final_pdb, resname="NA")
        # Sphere region is small; for charged fragments (e.g., NA) acceptance may drop as N grows.
        # Still require enough samples so "region ignored" would be astronomically unlikely.
        assert len(positions) >= 5, f"Expected >=5 inserted atoms, got {len(positions)}"
        _assert_all_in_sphere(positions, center=(25.0, 25.0, 25.0), radius=10.0)
    else:
        pytest.skip("gcmc_cpu executable not found")


def test_gcmc_region_box(tmp_path):
    """Box region must constrain insertion positions (file-driven, no stdout parsing)."""
    fragitp_path = TEST_DATA_DIR / "charmm36.ff/mol/na.itp"

    # Create a simple inp file with box region
    inp_content = f"""
version: gcmc_2.0
num_threads: 1
temperature: 300.0
mcsteps: 20
nprint: 20
fragname: NA
fragmuex: 100.0
fragconc: 1.0
gcmc_region: box 10.0 10.0 10.0 40.0 40.0 40.0
box: 50.0 50.0 50.0
cutoff: 12.0
mc_move_prob: 1 0 0 0
fragitp: {fragitp_path}
"""

    inp_file = tmp_path / "test_box.inp"
    inp_file.write_text(inp_content)

    # Import and run gcmc_cpu
    import subprocess

    # Run gcmc_cpu with the inp file
    gcmc_executable = GCMC_CPU_PATH
    if gcmc_executable.exists():
        out_prefix = tmp_path / "out" / "gcmc"
        out_prefix.parent.mkdir(parents=True, exist_ok=True)
        result = subprocess.run(
            [str(gcmc_executable), "--inp", str(inp_file), "--prefix", str(out_prefix), "--seed", "42"],
            cwd=tmp_path,
            capture_output=True,
            text=True
        )

        # Check that the simulation ran successfully
        assert result.returncode == 0, f"gcmc_cpu failed with error: {result.stderr}"

        final_pdb = Path(f"{out_prefix}_final.pdb")
        assert final_pdb.exists()

        positions = _pdb_atom_positions_angstrom(final_pdb, resname="NA")
        assert len(positions) >= 10, f"Expected >=10 inserted atoms, got {len(positions)}"
        _assert_all_in_box(positions, min_corner=(10.0, 10.0, 10.0), max_corner=(40.0, 40.0, 40.0))
    else:
        pytest.skip("gcmc_cpu executable not found")


def test_gcmc_region_cylinder(tmp_path):
    """Cylinder region must constrain insertion positions (file-driven, no stdout parsing)."""
    fragitp_path = TEST_DATA_DIR / "charmm36.ff/mol/na.itp"

    # Create a simple inp file with cylinder region
    inp_content = f"""
version: gcmc_2.0
num_threads: 1
temperature: 300.0
mcsteps: 20
nprint: 20
fragname: NA
fragmuex: 100.0
fragconc: 1.0
gcmc_region: cylinder 25.0 25.0 25.0 10.0 30.0 z
box: 50.0 50.0 50.0
cutoff: 12.0
mc_move_prob: 1 0 0 0
fragitp: {fragitp_path}
"""

    inp_file = tmp_path / "test_cylinder.inp"
    inp_file.write_text(inp_content)

    # Import and run gcmc_cpu
    import subprocess

    # Run gcmc_cpu with the inp file
    gcmc_executable = GCMC_CPU_PATH
    if gcmc_executable.exists():
        out_prefix = tmp_path / "out" / "gcmc"
        out_prefix.parent.mkdir(parents=True, exist_ok=True)
        result = subprocess.run(
            [str(gcmc_executable), "--inp", str(inp_file), "--prefix", str(out_prefix), "--seed", "42"],
            cwd=tmp_path,
            capture_output=True,
            text=True
        )

        # Check that the simulation ran successfully
        assert result.returncode == 0, f"gcmc_cpu failed with error: {result.stderr}"

        final_pdb = Path(f"{out_prefix}_final.pdb")
        assert final_pdb.exists()

        positions = _pdb_atom_positions_angstrom(final_pdb, resname="NA")
        assert len(positions) >= 10, f"Expected >=10 inserted atoms, got {len(positions)}"
        _assert_all_in_cylinder_z(positions, center=(25.0, 25.0, 25.0), radius=10.0, height=30.0)
    else:
        pytest.skip("gcmc_cpu executable not found")


if __name__ == "__main__":
    raise SystemExit("Run via pytest")

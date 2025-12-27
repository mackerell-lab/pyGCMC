"""
Feature tests for gcmc_cpu.

These tests are file-driven and validate behavior via output files / structured dumps,
not stdout/stderr log text.
"""

from __future__ import annotations

import json
import math
import subprocess
import tempfile
from pathlib import Path

import pytest


GCMC_CPU_PATH = Path(__file__).resolve().parents[3] / "build" / "bin" / "gcmc_cpu"
TEST_DATA_DIR = Path(__file__).resolve().parents[2] / "data"


def _require_binary() -> None:
    if not GCMC_CPU_PATH.exists():
        pytest.skip(f"gcmc_cpu not built: {GCMC_CPU_PATH}")


def _run_gcmc(inp: Path, *, cwd: Path, prefix: Path, extra_args: list[str] | None = None, timeout: int = 120):
    _require_binary()
    args = [str(GCMC_CPU_PATH), "--inp", str(inp), "--prefix", str(prefix)]
    if extra_args:
        args.extend(extra_args)
    result = subprocess.run(
        args,
        cwd=str(cwd),
        capture_output=True,
        text=True,
        timeout=timeout,
    )
    assert result.returncode == 0, result.stdout + result.stderr


def _pdb_atom_lines(pdb_path: Path) -> list[str]:
    return [line for line in pdb_path.read_text().splitlines() if line.startswith(("ATOM", "HETATM"))]


def _group_atoms_by_residue(pdb_path: Path) -> dict[tuple[str, str], list[dict[str, float | str]]]:
    residues: dict[tuple[str, str], list[dict[str, float | str]]] = {}
    for line in _pdb_atom_lines(pdb_path):
        if len(line) < 54:
            continue
        resname = line[17:20].strip().upper()
        resid = line[22:26].strip()
        atom_name = line[12:16].strip()
        residues.setdefault((resname, resid), []).append(
            {
                "name": atom_name,
                "x": float(line[30:38]),
                "y": float(line[38:46]),
                "z": float(line[46:54]),
            }
        )
    return residues


def _dist_angstrom(a: dict[str, float | str], b: dict[str, float | str]) -> float:
    dx = float(a["x"]) - float(b["x"])
    dy = float(a["y"]) - float(b["y"])
    dz = float(a["z"]) - float(b["z"])
    return math.sqrt(dx * dx + dy * dy + dz * dz)


class TestWaterInsertion:
    """Water insertion output sanity checks."""

    @pytest.fixture
    def temp_dir(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)

    def test_inserted_water_has_reasonable_geometry(self, temp_dir: Path):
        sol_itp = TEST_DATA_DIR / "charmm36.ff" / "mol" / "sol.itp"
        assert sol_itp.exists()

        inp = temp_dir / "water_struct.inp"
        inp.write_text(
            f"""
version:gcmc_2.0
box_size:30.0 30.0 30.0
temperature:300.0
fragitp:{sol_itp}
fragname:SOL
fragconc:55.0
fragmuex:1000.0

moves_per_step:1
mcsteps:1
nprint:1

attempt_prob_ins:1.0
attempt_prob_del:0.0
attempt_prob_trn:0.0
attempt_prob_rot:0.0
""".strip()
            + "\n"
        )

        prefix = temp_dir / "out" / "water"
        prefix.parent.mkdir(parents=True, exist_ok=True)
        _run_gcmc(inp, cwd=temp_dir, prefix=prefix, extra_args=["--seed", "42"], timeout=60)

        pdb = Path(f"{prefix}_final.pdb")
        assert pdb.exists()
        residues = _group_atoms_by_residue(pdb)

        water_res = [
            atoms
            for (resname, _), atoms in residues.items()
            if resname in {"SOL", "WAT"} and len(atoms) == 3
        ]
        assert water_res, "No 3-atom water residues found in final PDB"

        for atoms in water_res:
            oxygen = next((a for a in atoms if "O" in str(a["name"]).upper()), None)
            hydrogens = [a for a in atoms if "H" in str(a["name"]).upper()]
            assert oxygen is not None
            assert len(hydrogens) == 2
            for h in hydrogens:
                dist = _dist_angstrom(oxygen, h)
                assert 0.8 < dist < 1.2, f"Unexpected O-H distance: {dist} Å"

    def test_wdens_writes_density_dat(self, temp_dir: Path):
        sol_itp = TEST_DATA_DIR / "charmm36.ff" / "mol" / "sol.itp"
        assert sol_itp.exists()

        inp = temp_dir / "density.inp"
        inp.write_text(
            f"""
version:gcmc_2.0
box_size:30.0 30.0 30.0
temperature:300.0
fragitp:{sol_itp}
fragname:SOL
fragconc:55.0
fragmuex:-10.0
wdens:1.0
mcsteps:50
nprint:10
""".strip()
            + "\n"
        )

        prefix = temp_dir / "out" / "density"
        prefix.parent.mkdir(parents=True, exist_ok=True)
        _run_gcmc(inp, cwd=temp_dir, prefix=prefix, extra_args=["--seed", "7"], timeout=120)

        density_dat = Path(f"{prefix}_density.dat")
        assert density_dat.exists()

        data_rows = [
            line.strip()
            for line in density_dat.read_text().splitlines()
            if line.strip() and not line.startswith("#")
        ]
        assert data_rows, "density.dat contains no data rows"
        cols = data_rows[-1].split()
        assert len(cols) >= 2, f"Unexpected density.dat row: {data_rows[-1]}"
        assert float(cols[1]) >= 0.0


class TestCavityBias:
    """Cavity bias parameter parsing via --dump-params."""

    @pytest.fixture
    def temp_dir(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)

    def test_cavity_bias_and_exclusion_flags_in_dump_params(self, temp_dir: Path):
        sol_itp = TEST_DATA_DIR / "charmm36.ff" / "mol" / "sol.itp"
        assert sol_itp.exists()

        inp = temp_dir / "cavity.inp"
        inp.write_text(
            f"""
version:gcmc_2.0
box_size:20.0 20.0 20.0
cutoff:10.0
temperature:300.0
fragitp:{sol_itp}
fragname:SOL
fragconc:55.0
fragmuex:-10.0

use_cavity_bias:yes
cavity_grid_dx:2.0
probe_radius:1.4
exclude_protein_volume:yes
exclude_hydrogens_from_grid:yes
use_vdw_radius_for_grid:yes

mcsteps:0
nprint:1
""".strip()
            + "\n"
        )

        prefix = temp_dir / "out" / "cavity"
        prefix.parent.mkdir(parents=True, exist_ok=True)
        params_json = temp_dir / "out" / "params.json"

        _run_gcmc(inp, cwd=temp_dir, prefix=prefix, extra_args=["--dump-params", str(params_json)], timeout=120)

        params = json.loads(params_json.read_text())
        assert params["bias"]["use_cavity_bias"] is True
        assert float(params["space"]["grid_spacing_nm"]) == pytest.approx(0.2, abs=1e-6)
        assert float(params["bias"]["probe_radius_nm"]) == pytest.approx(0.14, abs=1e-6)
        assert params["space"]["exclude_protein_volume"] is True
        assert params["space"]["exclude_hydrogens_from_grid"] is True
        assert params["space"]["use_vdw_radius_for_grid"] is True


class TestRegionConstraints:
    """GCMC region conversion checks via --dump-params."""

    @pytest.fixture
    def temp_dir(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)

    @pytest.mark.parametrize(
        ("region_spec", "expected"),
        [
            ("sphere 25.0 25.0 25.0 10.0", ("sphere", [2.5, 2.5, 2.5, 1.0])),
            ("box 10.0 10.0 10.0 40.0 40.0 40.0", ("box", [1.0, 1.0, 1.0, 4.0, 4.0, 4.0])),
            ("cylinder 25.0 25.0 25.0 10.0 30.0 z", ("cylinder", [2.5, 2.5, 2.5, 1.0, 3.0, "z"])),
        ],
    )
    def test_region_is_converted_to_nm_in_dump_params(self, temp_dir: Path, region_spec: str, expected):
        sol_itp = TEST_DATA_DIR / "charmm36.ff" / "mol" / "sol.itp"
        assert sol_itp.exists()

        inp = temp_dir / "region.inp"
        inp.write_text(
            f"""
version:gcmc_2.0
box_size:50.0 50.0 50.0
temperature:300.0
fragitp:{sol_itp}
fragname:SOL
fragconc:55.0
fragmuex:-10.0
gcmc_region:{region_spec}
mcsteps:0
nprint:1
""".strip()
            + "\n"
        )

        prefix = temp_dir / "out" / "region"
        prefix.parent.mkdir(parents=True, exist_ok=True)
        params_json = temp_dir / "out" / "params.json"

        _run_gcmc(inp, cwd=temp_dir, prefix=prefix, extra_args=["--dump-params", str(params_json)], timeout=60)

        params = json.loads(params_json.read_text())
        toks = params["space"]["gcmc_region"].split()
        assert toks[0] == expected[0]
        if expected[0] == "cylinder":
            assert [float(x) for x in toks[1:6]] == pytest.approx([float(x) for x in expected[1][:5]], abs=1e-6)
            assert toks[6] == expected[1][5]
        else:
            assert [float(x) for x in toks[1:]] == pytest.approx([float(x) for x in expected[1]], abs=1e-6)


class TestTargetControl:
    """Target water count parsing via --dump-params."""

    @pytest.fixture
    def temp_dir(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)

    def test_target_numwaters_is_exposed_in_dump_params(self, temp_dir: Path):
        sol_itp = TEST_DATA_DIR / "charmm36.ff" / "mol" / "sol.itp"
        assert sol_itp.exists()

        inp = temp_dir / "target.inp"
        inp.write_text(
            f"""
version:gcmc_2.0
box_size:30.0 30.0 30.0
temperature:300.0
fragitp:{sol_itp}
fragname:SOL
fragconc:55.0
fragmuex:-10.0
target_numwaters:50
mcsteps:0
nprint:1
""".strip()
            + "\n"
        )

        prefix = temp_dir / "out" / "target"
        prefix.parent.mkdir(parents=True, exist_ok=True)
        params_json = temp_dir / "out" / "params.json"
        _run_gcmc(inp, cwd=temp_dir, prefix=prefix, extra_args=["--dump-params", str(params_json)], timeout=60)

        params = json.loads(params_json.read_text())
        assert int(params["fragment"]["target_num_waters"]) == 50


class TestParameterParsing:
    """Parameter parsing checks via --dump-params."""

    @pytest.fixture
    def temp_dir(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)

    def test_switching_and_pairlist_frequency_in_dump_params(self, temp_dir: Path):
        sol_itp = TEST_DATA_DIR / "charmm36.ff" / "mol" / "sol.itp"
        assert sol_itp.exists()

        inp = temp_dir / "params.inp"
        inp.write_text(
            f"""
version:gcmc_2.0
box_size:30.0 30.0 30.0
temperature:298.15
fragitp:{sol_itp}
fragname:SOL
fragconc:55.0
fragmuex:-10.0

use_switching:yes
switch_r_on:8.0
switch_r_off:10.0

pairlist_freq:500
use_group_cutoff:no

mcsteps:0
nprint:1
""".strip()
            + "\n"
        )

        prefix = temp_dir / "out" / "params"
        prefix.parent.mkdir(parents=True, exist_ok=True)
        params_json = temp_dir / "out" / "params.json"
        _run_gcmc(inp, cwd=temp_dir, prefix=prefix, extra_args=["--dump-params", str(params_json)], timeout=60)

        params = json.loads(params_json.read_text())
        assert params["mc"]["use_switching"] is True
        assert float(params["mc"]["switch_r_on_nm"]) == pytest.approx(0.8, abs=1e-6)
        assert float(params["mc"]["switch_r_off_nm"]) == pytest.approx(1.0, abs=1e-6)
        assert int(params["energy"]["pairlist_freq"]) == 500
        assert params["energy"]["use_group_cutoff"] is False

    def test_invalid_missing_fragment_definitions_fail(self, temp_dir: Path):
        _require_binary()
        inp = temp_dir / "invalid.inp"
        inp.write_text(
            """
version:gcmc_2.0
box_size:30.0 30.0 30.0
temperature:300.0
mcsteps:1
nprint:1
""".strip()
            + "\n"
        )

        prefix = temp_dir / "out" / "invalid"
        prefix.parent.mkdir(parents=True, exist_ok=True)
        result = subprocess.run(
            [str(GCMC_CPU_PATH), "--inp", str(inp), "--prefix", str(prefix)],
            cwd=str(temp_dir),
            capture_output=True,
            text=True,
            timeout=60,
        )
        assert result.returncode != 0


class TestOutputValidation:
    """Output file invariants that should not depend on logging."""

    @pytest.fixture
    def temp_dir(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)

    def test_pdb_atom_numbering_is_sequential_when_atoms_present(self, temp_dir: Path):
        sol_itp = TEST_DATA_DIR / "charmm36.ff" / "mol" / "sol.itp"
        assert sol_itp.exists()

        inp = temp_dir / "numbering.inp"
        inp.write_text(
            f"""
version:gcmc_2.0
box_size:20.0 20.0 20.0
cutoff:10.0
temperature:300.0
fragitp:{sol_itp}
fragname:SOL
fragconc:55.0
fragmuex:1000.0
moves_per_step:1
mcsteps:1
nprint:1
attempt_prob_ins:1.0
attempt_prob_del:0.0
attempt_prob_trn:0.0
attempt_prob_rot:0.0
""".strip()
            + "\n"
        )

        prefix = temp_dir / "out" / "num"
        prefix.parent.mkdir(parents=True, exist_ok=True)
        _run_gcmc(inp, cwd=temp_dir, prefix=prefix, extra_args=["--seed", "111"], timeout=60)

        pdb = Path(f"{prefix}_final.pdb")
        assert pdb.exists()
        atom_lines = _pdb_atom_lines(pdb)
        assert atom_lines, "No ATOM/HETATM lines captured"
        for i, line in enumerate(atom_lines, 1):
            assert int(line[6:11].strip()) == i

    def test_top_file_statistics_section_contains_fragment_row(self, temp_dir: Path):
        sol_itp = TEST_DATA_DIR / "charmm36.ff" / "mol" / "sol.itp"
        assert sol_itp.exists()

        inp = temp_dir / "stats.inp"
        inp.write_text(
            f"""
version:gcmc_2.0
box_size:20.0 20.0 20.0
cutoff:10.0
temperature:300.0
fragitp:{sol_itp}
fragname:SOL
fragconc:55.0
fragmuex:-10.0
mcsteps:0
nprint:1
""".strip()
            + "\n"
        )

        prefix = temp_dir / "out" / "stats"
        prefix.parent.mkdir(parents=True, exist_ok=True)
        _run_gcmc(inp, cwd=temp_dir, prefix=prefix, extra_args=["--seed", "222"], timeout=60)

        top_file = Path(f"{prefix}_final.top")
        assert top_file.exists()
        top_content = top_file.read_text()

        assert "[ system ]" in top_content
        assert "[ molecules ]" in top_content
        assert "[ fragments ]" in top_content
        assert "[ statistics ]" in top_content

        frag_row = None
        in_frag = False
        for line in top_content.splitlines():
            if line.strip() == "[ fragments ]":
                in_frag = True
                continue
            if in_frag:
                if line.startswith("["):
                    break
                if not line.strip() or line.lstrip().startswith(";"):
                    continue
                cols = line.split()
                if cols and cols[0].upper() in {"WAT", "SOL"} and len(cols) >= 5:
                    frag_row = cols
                    break
        assert frag_row is not None, "Missing fragment row in [ fragments ] section"
        chem_pot_kj = float(frag_row[-1])
        assert chem_pot_kj == pytest.approx(-10.0 * 4.184, abs=0.02)

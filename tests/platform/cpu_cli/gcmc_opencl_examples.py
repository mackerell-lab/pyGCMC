"""
End-to-end compatibility checks using gcmc_opencl (gcmc_gpu-style) example decks.

These tests are file-driven and validate behavior via output files / structured dumps,
not stdout/stderr log text.
"""

from __future__ import annotations

import json
import shutil
import subprocess
from pathlib import Path

import pytest


def _read_cryst1_box_angstrom(pdb_path: Path) -> tuple[float, float, float]:
    for line in pdb_path.read_text().splitlines():
        if line.startswith("CRYST1"):
            parts = line.split()
            assert len(parts) >= 4, f"Unexpected CRYST1 format: {line}"
            return float(parts[1]), float(parts[2]), float(parts[3])
    raise AssertionError(f"CRYST1 not found in {pdb_path}")


def _count_atom_records(pdb_path: Path) -> int:
    return sum(
        1
        for line in pdb_path.read_text().splitlines()
        if line.startswith(("ATOM", "HETATM"))
    )


def _read_inp_box_size_angstrom(inp_path: Path) -> tuple[float, float, float]:
    for raw in inp_path.read_text().splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if not (line.startswith("box_size:") or line.startswith("box:")):
            continue
        parts = line.split(":", 1)[1].split()
        assert len(parts) >= 3, f"Unexpected box_size format: {raw}"
        return float(parts[0]), float(parts[1]), float(parts[2])
    raise AssertionError(f"box_size not found in {inp_path}")


def _write_smoke_inp_from_opencl_example(
    src_inp: Path,
    dst_inp: Path,
    *,
    mcsteps: int,
    overrides: dict[str, str] | None = None,
) -> None:
    overrides = overrides or {}
    lines: list[str] = []
    for raw in src_inp.read_text().splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            lines.append(raw)
            continue
        key = line.split(":", 1)[0].strip()
        if key in overrides:
            lines.append(f"{key}:{overrides[key]}")
            continue
        if key == "mcsteps":
            lines.append(f"mcsteps:{mcsteps}")
            continue
        if key == "nprint":
            lines.append(f"nprint:{max(1, mcsteps)}")
            continue
        lines.append(raw)
    dst_inp.write_text("\n".join(lines).rstrip() + "\n")


def _symlink_forcefield_dir(work: Path, test_data_dir: Path) -> None:
    link = work / "charmm36.ff"
    if link.exists() or link.is_symlink():
        link.unlink()
    link.symlink_to(test_data_dir / "charmm36.ff", target_is_directory=True)


def test_opencl_twowater_example_runs_and_preserves_pdb_cryst1(gcmc_cpu, test_data_dir, temp_dir):
    src = test_data_dir / "gcmc_opencl_examples" / "twowater"
    assert (src / "gcmc.inp").exists()
    assert (src / "twowater.top").exists()
    assert (src / "twowater.pdb").exists()

    work = Path(temp_dir) / "opencl_twowater"
    work.mkdir(parents=True, exist_ok=True)

    shutil.copy(src / "twowater.top", work / "twowater.top")
    shutil.copy(src / "twowater.pdb", work / "twowater.pdb")
    _write_smoke_inp_from_opencl_example(src / "gcmc.inp", work / "run.inp", mcsteps=25)
    _symlink_forcefield_dir(work, Path(test_data_dir))

    out_prefix = work / "out" / "gcmc"
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    params_json = work / "out" / "params.json"

    cryst_in = _read_cryst1_box_angstrom(work / "twowater.pdb")

    result = subprocess.run(
        [
            gcmc_cpu,
            "--inp",
            str(work / "run.inp"),
            "--prefix",
            str(out_prefix),
            "--seed",
            "123",
            "--dump-params",
            str(params_json),
        ],
        cwd=str(work),
        capture_output=True,
        text=True,
        timeout=60,
    )
    assert result.returncode == 0, result.stdout + result.stderr

    final_pdb = Path(f"{out_prefix}_final.pdb")
    assert final_pdb.exists()
    assert _read_cryst1_box_angstrom(final_pdb) == pytest.approx(cryst_in, abs=1e-3)

    # gcmc_gpu-style op_* outputs should exist as additional final snapshots.
    op_pdb = work / "out.pdb"
    op_top = work / "out.top"
    assert op_pdb.exists()
    assert op_top.exists()
    assert op_top.read_text().strip(), "op_top unexpectedly empty"
    assert _count_atom_records(op_pdb) == _count_atom_records(final_pdb)

    params = json.loads(params_json.read_text())
    assert params["basic"]["inp_units"] == "gcmc_gpu"


def test_opencl_waterbox_hollow_example_smoke_runs(gcmc_cpu, test_data_dir, temp_dir):
    src = test_data_dir / "gcmc_opencl_examples" / "waterbox_hollow"
    assert (src / "gcmc.inp").exists()
    assert (src / "waterbox.top").exists()
    assert (src / "waterbox.pdb").exists()

    work = Path(temp_dir) / "opencl_waterbox_hollow"
    work.mkdir(parents=True, exist_ok=True)

    shutil.copy(src / "waterbox.top", work / "waterbox.top")
    shutil.copy(src / "waterbox.pdb", work / "waterbox.pdb")
    _write_smoke_inp_from_opencl_example(src / "gcmc.inp", work / "run.inp", mcsteps=10)
    _symlink_forcefield_dir(work, Path(test_data_dir))

    out_prefix = work / "out" / "gcmc"
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    cryst_in = _read_cryst1_box_angstrom(work / "waterbox.pdb")

    result = subprocess.run(
        [
            gcmc_cpu,
            "--inp",
            str(work / "run.inp"),
            "--prefix",
            str(out_prefix),
            "--seed",
            "321",
        ],
        cwd=str(work),
        capture_output=True,
        text=True,
        timeout=120,
    )
    assert result.returncode == 0, result.stdout + result.stderr

    final_pdb = Path(f"{out_prefix}_final.pdb")
    assert final_pdb.exists()
    assert _read_cryst1_box_angstrom(final_pdb) == pytest.approx(cryst_in, abs=1e-3)

def test_opencl_waterbox_example_smoke_runs_and_enables_cavity_bias(gcmc_cpu, test_data_dir, temp_dir):
    src = test_data_dir / "gcmc_opencl_examples" / "waterbox"
    assert (src / "gcmc.inp").exists()
    assert (src / "waterbox.top").exists()
    assert (src / "waterbox.pdb").exists()

    work = Path(temp_dir) / "opencl_waterbox"
    work.mkdir(parents=True, exist_ok=True)

    shutil.copy(src / "waterbox.top", work / "waterbox.top")
    shutil.copy(src / "waterbox.pdb", work / "waterbox.pdb")
    _write_smoke_inp_from_opencl_example(src / "gcmc.inp", work / "run.inp", mcsteps=5)
    _symlink_forcefield_dir(work, Path(test_data_dir))

    out_prefix = work / "out" / "gcmc"
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    params_json = work / "out" / "params.json"

    expected_box_angstrom = _read_inp_box_size_angstrom(work / "run.inp")

    result = subprocess.run(
        [
            gcmc_cpu,
            "--inp",
            str(work / "run.inp"),
            "--prefix",
            str(out_prefix),
            "--seed",
            "2024",
            "--dump-params",
            str(params_json),
        ],
        cwd=str(work),
        capture_output=True,
        text=True,
        timeout=120,
    )
    assert result.returncode == 0, result.stdout + result.stderr

    final_pdb = Path(f"{out_prefix}_final.pdb")
    assert final_pdb.exists()
    assert _read_cryst1_box_angstrom(final_pdb) == pytest.approx(expected_box_angstrom, abs=1e-3)

    op_pdb = work / "out.pdb"
    op_top = work / "out.top"
    assert op_pdb.exists()
    assert op_top.exists()
    assert op_top.read_text().strip(), "op_top unexpectedly empty"

    params = json.loads(params_json.read_text())
    assert params["basic"]["inp_units"] == "gcmc_gpu"
    assert float(params["space"]["grid_spacing_nm"]) == pytest.approx(0.1, abs=1e-6)
    assert float(params["energy"]["fragment_cutoff_nm"]) == pytest.approx(0.8, abs=1e-6)
    assert params["bias"]["use_cavity_bias"] is True
    assert params["space"]["use_vdw_radius_for_grid"] is True
    assert params["space"]["exclude_hydrogens_from_grid"] is False


def test_opencl_benz_example_smoke_runs_and_parses_multi_fragment_lists(gcmc_cpu, test_data_dir, temp_dir):
    src = test_data_dir / "gcmc_opencl_examples" / "benz"
    assert (src / "gcmc.inp").exists()
    assert (src / "waterbox.top").exists()
    assert (src / "waterbox.pdb").exists()

    work = Path(temp_dir) / "opencl_benz"
    work.mkdir(parents=True, exist_ok=True)

    shutil.copy(src / "waterbox.top", work / "waterbox.top")
    shutil.copy(src / "waterbox.pdb", work / "waterbox.pdb")
    _write_smoke_inp_from_opencl_example(src / "gcmc.inp", work / "run.inp", mcsteps=0)
    _symlink_forcefield_dir(work, Path(test_data_dir))

    out_prefix = work / "out" / "gcmc"
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    params_json = work / "out" / "params.json"

    cryst_in = _read_cryst1_box_angstrom(work / "waterbox.pdb")

    result = subprocess.run(
        [
            gcmc_cpu,
            "--inp",
            str(work / "run.inp"),
            "--prefix",
            str(out_prefix),
            "--seed",
            "17",
            "--dump-params",
            str(params_json),
        ],
        cwd=str(work),
        capture_output=True,
        text=True,
        timeout=180,
    )
    assert result.returncode == 0, result.stdout + result.stderr

    final_pdb = Path(f"{out_prefix}_final.pdb")
    assert final_pdb.exists()
    assert _read_cryst1_box_angstrom(final_pdb) == pytest.approx(cryst_in, abs=1e-3)

    op_pdb = work / "out.pdb"
    op_top = work / "out.top"
    assert op_pdb.exists()
    assert op_top.exists()
    assert op_top.read_text().strip(), "op_top unexpectedly empty"
    assert _count_atom_records(op_pdb) == _count_atom_records(final_pdb)

    params = json.loads(params_json.read_text())
    assert params["basic"]["inp_units"] == "gcmc_gpu"
    assert float(params["space"]["grid_spacing_nm"]) == pytest.approx(0.1, abs=1e-6)
    assert float(params["energy"]["fragment_cutoff_nm"]) == pytest.approx(0.8, abs=1e-6)
    assert params["bias"]["use_cavity_bias"] is False
    assert params["bias"]["use_conf_bias"] is False
    assert [float(x) for x in params["fragment"]["conc_list_M"]] == pytest.approx([1.0, 50.0], abs=1e-6)
    expected_muex_kj = [-0.79 * 4.184, -5.6 * 4.184]
    assert [float(x) for x in params["fragment"]["muex_list_kj_mol"]] == pytest.approx(
        expected_muex_kj, abs=1e-4
    )


def test_opencl_protein_example_dump_params_parses_complex_deck(gcmc_cpu, test_data_dir, temp_dir):
    src = test_data_dir / "gcmc_opencl_examples" / "protein"
    assert (src / "gcmc.inp").exists()
    assert (src / "181L_apo_silcs.1.top").exists()
    assert (src / "181L_apo_silcs.1.pdb").exists()
    assert (src / "posre.itp").exists()
    assert (src / "posre_protein_ca.itp").exists()

    work = Path(temp_dir) / "opencl_protein"
    work.mkdir(parents=True, exist_ok=True)

    shutil.copy(src / "181L_apo_silcs.1.top", work / "181L_apo_silcs.1.top")
    shutil.copy(src / "181L_apo_silcs.1.pdb", work / "181L_apo_silcs.1.pdb")
    shutil.copy(src / "posre.itp", work / "posre.itp")
    shutil.copy(src / "posre_protein_ca.itp", work / "posre_protein_ca.itp")

    _write_smoke_inp_from_opencl_example(
        src / "gcmc.inp",
        work / "run.inp",
        mcsteps=0,
        overrides={"op_top": "out.top", "op_pdb": "out.pdb"},
    )
    _symlink_forcefield_dir(work, Path(test_data_dir))

    out_prefix = work / "out" / "gcmc"
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    params_json = work / "out" / "params.json"

    cryst_in = _read_cryst1_box_angstrom(work / "181L_apo_silcs.1.pdb")

    result = subprocess.run(
        [
            gcmc_cpu,
            "--inp",
            str(work / "run.inp"),
            "--prefix",
            str(out_prefix),
            "--seed",
            "7",
            "--dump-params",
            str(params_json),
        ],
        cwd=str(work),
        capture_output=True,
        text=True,
        timeout=240,
    )
    assert result.returncode == 0, result.stdout + result.stderr

    final_pdb = Path(f"{out_prefix}_final.pdb")
    assert final_pdb.exists()
    assert _read_cryst1_box_angstrom(final_pdb) == pytest.approx(cryst_in, abs=1e-3)

    params = json.loads(params_json.read_text())
    assert params["basic"]["inp_units"] == "gcmc_gpu"
    assert float(params["space"]["grid_spacing_nm"]) == pytest.approx(0.1, abs=1e-6)
    assert float(params["energy"]["fragment_cutoff_nm"]) == pytest.approx(0.8, abs=1e-6)
    assert float(params["energy"]["protein_cutoff_nm"]) == pytest.approx(0.8, abs=1e-6)
    assert params["bias"]["use_cavity_bias"] is True
    assert params["bias"]["use_conf_bias"] is True
    assert int(params["bias"]["num_conf_bias_trials"]) == 10
    assert params["fragment"]["use_number_water_nbar"] is True


def test_opencl_test_mg_example_smoke_runs_and_converts_cutoff_units(gcmc_cpu, test_data_dir, temp_dir):
    src = test_data_dir / "gcmc_opencl_examples" / "test"
    assert (src / "gcmc.inp").exists()
    assert (src / "solution.1.top").exists()
    assert (src / "solution.1.pdb").exists()

    work = Path(temp_dir) / "opencl_test_mg"
    work.mkdir(parents=True, exist_ok=True)

    shutil.copy(src / "solution.1.top", work / "solution.1.top")
    shutil.copy(src / "solution.1.pdb", work / "solution.1.pdb")
    _write_smoke_inp_from_opencl_example(src / "gcmc.inp", work / "run.inp", mcsteps=20)
    _symlink_forcefield_dir(work, Path(test_data_dir))

    out_prefix = work / "out" / "gcmc"
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    params_json = work / "out" / "params.json"

    cryst_in = _read_cryst1_box_angstrom(work / "solution.1.pdb")

    result = subprocess.run(
        [
            gcmc_cpu,
            "--inp",
            str(work / "run.inp"),
            "--prefix",
            str(out_prefix),
            "--seed",
            "999",
            "--dump-params",
            str(params_json),
        ],
        cwd=str(work),
        capture_output=True,
        text=True,
        timeout=60,
    )
    assert result.returncode == 0, result.stdout + result.stderr

    final_pdb = Path(f"{out_prefix}_final.pdb")
    assert final_pdb.exists()
    assert _read_cryst1_box_angstrom(final_pdb) == pytest.approx(cryst_in, abs=1e-3)

    op_pdb = work / "out.pdb"
    op_top = work / "out.top"
    assert op_pdb.exists()
    assert op_top.exists()
    assert op_top.read_text().strip(), "op_top unexpectedly empty"

    assert params_json.exists()
    params = json.loads(params_json.read_text())

    assert params["basic"]["inp_units"] == "gcmc_gpu"
    assert float(params["space"]["grid_spacing_nm"]) == pytest.approx(0.1, abs=1e-6)
    assert float(params["space"]["cutoff_nm"]) == pytest.approx(1.2, abs=1e-6)
    assert float(params["energy"]["fragment_cutoff_nm"]) == pytest.approx(1.2, abs=1e-6)
    assert params["space"]["use_vdw_radius_for_grid"] is True
    assert params["space"]["exclude_hydrogens_from_grid"] is False

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


def _extract_inp_keys(inp_path: Path) -> list[str]:
    keys: list[str] = []
    for raw in inp_path.read_text().splitlines():
        line = raw.strip()
        if not line or line.startswith("#") or ":" not in line:
            continue
        key = line.split(":", 1)[0].strip()
        if key and key not in keys:
            keys.append(key)
    return sorted(keys)


def _load_inp_key_inventory(test_data_dir: Path) -> dict:
    inventory_path = test_data_dir / "gcmc_opencl_examples" / "inp_key_inventory.json"
    return json.loads(inventory_path.read_text())


def _assert_inp_key_inventory(inp_path: Path, expected_keys: list[str]) -> None:
    actual = _extract_inp_keys(inp_path)
    assert actual == sorted(expected_keys)


def _assert_inp_key_classification(params: dict, *, example: str, inventory: dict) -> None:
    expected = inventory["examples"][example]
    expected_unknown = set(expected.get("unknown", []))
    expected_keys = set(expected["keys"])
    ignored_keys = set(inventory["ignored_keys"])
    expected_ignored = expected_keys & ignored_keys
    assert set(params["basic"]["unknown_inp_keys"]) == expected_unknown
    assert set(params["basic"]["ignored_inp_keys"]) == expected_ignored


def _sanitize_fragment_name(name: str) -> str:
    out = []
    for ch in name:
        if ch.isalnum() or ch in "_-":
            out.append(ch)
        else:
            out.append("_")
    return "".join(out)


def _extract_frag_names_and_muex(inp_path: Path) -> tuple[list[str], list[float]]:
    frag_names: list[str] = []
    frag_muex: list[float] = []
    for raw in inp_path.read_text().splitlines():
        line = raw.strip()
        if not line or line.startswith("#") or ":" not in line:
            continue
        key, value = line.split(":", 1)
        key = key.strip()
        parts = value.split()
        if key == "fragname":
            frag_names.extend(parts)
        elif key == "fragmuex":
            frag_muex.extend(float(x) for x in parts)
    assert frag_names, f"No fragname entries found in {inp_path}"
    assert frag_muex, f"No fragmuex entries found in {inp_path}"
    assert len(frag_names) == len(frag_muex)
    return frag_names, frag_muex


def _read_statistics_n_total(stats_path: Path) -> list[int]:
    n_total: list[int] = []
    for line in stats_path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 3:
            continue
        n_total.append(int(parts[2]))
    return n_total


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


def _symlink_forcefield_dir_under(work: Path, test_data_dir: Path, rel_dir: str) -> None:
    base = work / rel_dir
    base.mkdir(parents=True, exist_ok=True)
    link = base / "charmm36.ff"
    if link.exists() or link.is_symlink():
        link.unlink()
    link.symlink_to(test_data_dir / "charmm36.ff", target_is_directory=True)


def test_opencl_inp_key_inventory_matches_examples(test_data_dir):
    inventory = _load_inp_key_inventory(Path(test_data_dir))
    examples = inventory.get("examples", {})
    for name, meta in examples.items():
        inp_path = Path(test_data_dir) / "gcmc_opencl_examples" / name / meta["inp"]
        assert inp_path.exists(), f"Missing INP for {name}: {inp_path}"
        _assert_inp_key_inventory(inp_path, meta["keys"])


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
    inventory = _load_inp_key_inventory(Path(test_data_dir))
    _assert_inp_key_classification(params, example="twowater", inventory=inventory)


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
            "321",
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
    assert _read_cryst1_box_angstrom(final_pdb) == pytest.approx(cryst_in, abs=1e-3)
    params = json.loads(params_json.read_text())
    inventory = _load_inp_key_inventory(Path(test_data_dir))
    _assert_inp_key_classification(params, example="waterbox_hollow", inventory=inventory)

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
    inventory = _load_inp_key_inventory(Path(test_data_dir))
    _assert_inp_key_classification(params, example="waterbox", inventory=inventory)


def test_opencl_benz_example_smoke_runs_and_parses_multi_fragment_lists(gcmc_cpu, test_data_dir, temp_dir):
    src = test_data_dir / "gcmc_opencl_examples" / "benz"
    assert (src / "gcmc.inp").exists()
    assert (src / "waterbox.top").exists()
    assert (src / "waterbox.pdb").exists()

    work = Path(temp_dir) / "opencl_benz"
    work.mkdir(parents=True, exist_ok=True)

    shutil.copy(src / "waterbox.top", work / "waterbox.top")
    shutil.copy(src / "waterbox.pdb", work / "waterbox.pdb")
    _write_smoke_inp_from_opencl_example(
        src / "gcmc.inp",
        work / "run.inp",
        mcsteps=3,
        overrides={"nprint": "1"},
    )
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
    inventory = _load_inp_key_inventory(Path(test_data_dir))
    _assert_inp_key_classification(params, example="benz", inventory=inventory)

    frag_names, frag_muex = _extract_frag_names_and_muex(src / "gcmc.inp")
    out_dir = out_prefix.parent
    active_line_counts: list[int] = []
    active_last_counts: list[int] = []
    for frag_name, expected_muex in zip(frag_names, frag_muex):
        sanitized = _sanitize_fragment_name(frag_name)
        active_path = out_dir / f"active_{sanitized}.dat"
        muex_path = out_dir / f"muex_{sanitized}.dat"
        assert active_path.exists()
        assert muex_path.exists()

        active_lines = [line.strip() for line in active_path.read_text().splitlines() if line.strip()]
        muex_lines = [line.strip() for line in muex_path.read_text().splitlines() if line.strip()]
        assert active_lines, f"{active_path} unexpectedly empty"
        assert muex_lines, f"{muex_path} unexpectedly empty"

        active_values = [int(val) for val in active_lines]
        muex_values = [float(val) for val in muex_lines]
        assert all(val >= 0 for val in active_values)
        assert all(abs(val - expected_muex) <= 0.02 for val in muex_values)
        active_line_counts.append(len(active_values))
        active_last_counts.append(active_values[-1])

    stats_path = Path(f"{out_prefix}_statistics.dat")
    assert stats_path.exists()
    n_total = _read_statistics_n_total(stats_path)
    assert n_total, f"{stats_path} unexpectedly empty"
    assert all(count == len(n_total) for count in active_line_counts)
    assert sum(active_last_counts) == n_total[-1]


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
    inventory = _load_inp_key_inventory(Path(test_data_dir))
    _assert_inp_key_classification(params, example="protein", inventory=inventory)


def test_opencl_protein_example_active_muex_outputs(gcmc_cpu, test_data_dir, temp_dir):
    """
    Complex deck regression: active_*/muex_* outputs must exist and stay consistent.
    """
    src = test_data_dir / "gcmc_opencl_examples" / "protein"
    assert (src / "gcmc.inp").exists()
    assert (src / "181L_apo_silcs.1.top").exists()
    assert (src / "181L_apo_silcs.1.pdb").exists()
    assert (src / "posre.itp").exists()
    assert (src / "posre_protein_ca.itp").exists()

    work = Path(temp_dir) / "opencl_protein_active_muex"
    work.mkdir(parents=True, exist_ok=True)

    shutil.copy(src / "181L_apo_silcs.1.top", work / "181L_apo_silcs.1.top")
    shutil.copy(src / "181L_apo_silcs.1.pdb", work / "181L_apo_silcs.1.pdb")
    shutil.copy(src / "posre.itp", work / "posre.itp")
    shutil.copy(src / "posre_protein_ca.itp", work / "posre_protein_ca.itp")

    _write_smoke_inp_from_opencl_example(
        src / "gcmc.inp",
        work / "run.inp",
        mcsteps=3,
        overrides={"op_top": "out.top", "op_pdb": "out.pdb"},
    )
    _symlink_forcefield_dir(work, Path(test_data_dir))

    out_prefix = work / "out" / "gcmc"
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    params_json = work / "out" / "params.json"

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
        timeout=300,
    )
    assert result.returncode == 0, result.stdout + result.stderr

    final_pdb = Path(f"{out_prefix}_final.pdb")
    assert final_pdb.exists()
    initial_pdb = work / "181L_apo_silcs.1.pdb"

    frag_names, frag_muex = _extract_frag_names_and_muex(work / "run.inp")
    assert len(frag_names) >= 2, "Expected multi-fragment deck for active/muex contract"
    assert params_json.exists()
    params = json.loads(params_json.read_text())
    assert params["bias"]["use_cavity_bias"] is True
    assert params["bias"]["use_conf_bias"] is True
    assert int(params["bias"]["num_conf_bias_trials"]) >= 2
    out_dir = out_prefix.parent
    active_line_counts: list[int] = []
    active_last_counts: list[int] = []
    resname_dupes: set[str] = set()
    resname_map: dict[str, str] = {}
    for name in frag_names:
        key = name.strip()[:3].upper()
        if key in resname_map:
            resname_dupes.add(key)
        resname_map[key] = name
    for frag_name, expected_muex in zip(frag_names, frag_muex):
        sanitized = _sanitize_fragment_name(frag_name)
        active_path = out_dir / f"active_{sanitized}.dat"
        muex_path = out_dir / f"muex_{sanitized}.dat"
        assert active_path.exists()
        assert muex_path.exists()

        active_lines = [line.strip() for line in active_path.read_text().splitlines() if line.strip()]
        muex_lines = [line.strip() for line in muex_path.read_text().splitlines() if line.strip()]
        assert active_lines, f"{active_path} unexpectedly empty"
        assert muex_lines, f"{muex_path} unexpectedly empty"

        active_values = [int(val) for val in active_lines]
        muex_values = [float(val) for val in muex_lines]
        assert all(val >= 0 for val in active_values)
        assert all(abs(val - expected_muex) <= 0.02 for val in muex_values)
        active_line_counts.append(len(active_values))
        active_last_counts.append(active_values[-1])

        resname_key = frag_name.strip()[:3].upper()
        if resname_key and resname_key not in resname_dupes:
            baseline = _count_residues_by_resname(initial_pdb, resname_key)
            if baseline == 0:
                pdb_count = _count_residues_by_resname(final_pdb, resname_key)
                assert active_values[-1] == pdb_count, f"{frag_name} count mismatch"

    stats_path = Path(f"{out_prefix}_statistics.dat")
    assert stats_path.exists()
    n_total = _read_statistics_n_total(stats_path)
    assert n_total, f"{stats_path} unexpectedly empty"
    assert all(count == len(n_total) for count in active_line_counts)
    assert sum(active_last_counts) == n_total[-1]


def test_opencl_cdk2_example_smoke_runs_and_reports_ignored_gcmc_cutoff(gcmc_cpu, test_data_dir, temp_dir):
    src = test_data_dir / "gcmc_opencl_examples" / "cdk2"
    assert (src / "test_conf.inp").exists()
    assert (src / "1h1q_silcs.1.top").exists()
    assert (src / "1h1q_silcs.1.pdb").exists()

    work = Path(temp_dir) / "opencl_cdk2"
    work.mkdir(parents=True, exist_ok=True)

    shutil.copy(src / "1h1q_silcs.1.top", work / "1h1q_silcs.1.top")
    shutil.copy(src / "1h1q_silcs.1.pdb", work / "1h1q_silcs.1.pdb")
    _write_smoke_inp_from_opencl_example(src / "test_conf.inp", work / "run.inp", mcsteps=0)
    _symlink_forcefield_dir(work, Path(test_data_dir))

    out_prefix = work / "out" / "gcmc"
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    params_json = work / "out" / "params.json"

    result = subprocess.run(
        [
            gcmc_cpu,
            "--inp",
            str(work / "run.inp"),
            "--prefix",
            str(out_prefix),
            "--seed",
            "101",
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

    op_pdb = work / "test.pdb"
    op_top = work / "test.top"
    assert op_pdb.exists()
    assert op_top.exists()
    assert op_top.read_text().strip(), "op_top unexpectedly empty"

    params = json.loads(params_json.read_text())
    assert params["basic"]["inp_units"] == "gcmc_gpu"
    assert float(params["space"]["grid_spacing_nm"]) == pytest.approx(0.1, abs=1e-6)
    assert float(params["space"]["cutoff_nm"]) == pytest.approx(1.2, abs=1e-6)
    assert params["energy"]["use_group_cutoff"] is True
    assert params["bias"]["use_cavity_bias"] is True
    assert params["bias"]["use_conf_bias"] is True
    assert int(params["bias"]["num_conf_bias_trials"]) == 10
    assert params["fragment"]["use_number_water_nbar"] is True
    assert params["basic"]["unknown_inp_keys"] == []
    ignored = set(params["basic"]["ignored_inp_keys"])
    assert "use_gcmc_cutoff" in ignored
    inventory = _load_inp_key_inventory(Path(test_data_dir))
    _assert_inp_key_classification(params, example="cdk2", inventory=inventory)


def test_opencl_lysozyme_example_smoke_runs_and_reports_ignored_map_keys(gcmc_cpu, test_data_dir, temp_dir):
    src = test_data_dir / "gcmc_opencl_examples" / "lysozyme"
    assert (src / "gcmc.0.inp").exists()
    assert (src / "181L_apo_silcs.1.top").exists()
    assert (src / "181L_apo_silcs.1.pdb").exists()
    assert (src / "181L_apo_silcs.1.equil.rec.pdb").exists()

    work = Path(temp_dir) / "opencl_lysozyme"
    work.mkdir(parents=True, exist_ok=True)

    examples_dir = work / "examples" / "lysozyme"
    examples_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy(src / "181L_apo_silcs.1.top", examples_dir / "181L_apo_silcs.1.top")
    shutil.copy(src / "181L_apo_silcs.1.pdb", examples_dir / "181L_apo_silcs.1.pdb")
    shutil.copy(
        src / "181L_apo_silcs.1.equil.rec.pdb",
        examples_dir / "181L_apo_silcs.1.equil.rec.pdb",
    )
    _write_smoke_inp_from_opencl_example(src / "gcmc.0.inp", work / "run.inp", mcsteps=0)
    _symlink_forcefield_dir_under(work, Path(test_data_dir), "data")

    out_prefix = work / "out" / "gcmc"
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    params_json = work / "out" / "params.json"

    result = subprocess.run(
        [
            gcmc_cpu,
            "--inp",
            str(work / "run.inp"),
            "--prefix",
            str(out_prefix),
            "--seed",
            "88",
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

    op_top = examples_dir / "181L_apo_silcs.1.gc.0.top"
    op_pdb = examples_dir / "181L_apo_silcs.1.gc.0.pdb"
    assert op_top.exists()
    assert op_pdb.exists()
    assert op_top.read_text().strip(), "op_top unexpectedly empty"

    params = json.loads(params_json.read_text())
    assert params["basic"]["inp_units"] == "gcmc_gpu"
    assert float(params["space"]["grid_spacing_nm"]) == pytest.approx(0.1, abs=1e-6)
    assert float(params["space"]["cutoff_nm"]) == pytest.approx(1.2, abs=1e-6)
    assert params["basic"]["unknown_inp_keys"] == []
    ignored = set(params["basic"]["ignored_inp_keys"])
    assert "initcycle" in ignored
    assert "conserve_frags" in ignored
    assert "map_generation" in ignored
    assert "map_filename_prefix" in ignored
    inventory = _load_inp_key_inventory(Path(test_data_dir))
    _assert_inp_key_classification(params, example="lysozyme", inventory=inventory)


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
    inventory = _load_inp_key_inventory(Path(test_data_dir))
    _assert_inp_key_classification(params, example="test", inventory=inventory)

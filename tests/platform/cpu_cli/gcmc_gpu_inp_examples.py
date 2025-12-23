"""
End-to-end compatibility checks using gcmc_gpu-style INP decks.

Goal: input/output/parameter semantics compatibility (not bitwise-identical trajectories).
Tests are file-driven and avoid matching stdout/stderr log text.
"""

from __future__ import annotations

import json
import subprocess
from pathlib import Path

import pytest


def _render_template(template_path: Path, *, data_dir: Path, dst_dir: Path) -> Path:
    content = template_path.read_text().replace("{DATA_DIR}", str(data_dir))
    out_path = dst_dir / template_path.name.replace("_template", "")
    out_path.write_text(content)
    return out_path


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


def _first_accept_record(path: Path, *, move: str, species: str) -> dict:
    records = [json.loads(line) for line in path.read_text().splitlines() if line.strip()]
    want_move = move.strip().lower()
    want_species = species.strip().upper()
    for rec in records:
        if str(rec.get("move", "")).strip().lower() != want_move:
            continue
        if str(rec.get("species", "")).strip().upper() != want_species:
            continue
        return rec
    raise AssertionError(f"No {move}/{species} record found in {path}; first records: {records[:3]}")


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
    inventory_path = test_data_dir / "gcmc_gpu_examples" / "inp_key_inventory.json"
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


def test_gcmc_gpu_inp_key_inventory_matches_templates(test_data_dir):
    inventory = _load_inp_key_inventory(Path(test_data_dir))
    examples = inventory.get("examples", {})
    for name, meta in examples.items():
        inp_path = Path(test_data_dir) / "gcmc_gpu_examples" / meta["inp"]
        assert inp_path.exists(), f"Missing INP for {name}: {inp_path}"
        _assert_inp_key_inventory(inp_path, meta["keys"])


def test_gcmc_gpu_water_tip3p_runs_and_honors_op_outputs(gcmc_cpu, test_data_dir, temp_dir):
    templates = test_data_dir / "gcmc_gpu_examples"
    template = templates / "water_tip3p_gcmc_gpu_template.inp"
    assert template.exists()

    work = Path(temp_dir) / "gcmc_gpu_water_tip3p"
    work.mkdir(parents=True, exist_ok=True)

    inp = _render_template(template, data_dir=test_data_dir, dst_dir=work)
    out_prefix = work / "out" / "gcmc"
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    accept_log = work / "out" / "acceptance.jsonl"
    params_json = work / "out" / "params.json"

    result = subprocess.run(
        [
            gcmc_cpu,
            "--inp",
            str(inp),
            "--prefix",
            str(out_prefix),
            "--dump-params",
            str(params_json),
            "--dump-accept",
            str(accept_log),
            "--store-probabilities",
        ],
        cwd=str(work),
        capture_output=True,
        text=True,
        timeout=60,
    )
    assert result.returncode == 0, result.stdout + result.stderr

    final_pdb = Path(f"{out_prefix}_final.pdb")
    assert final_pdb.exists()
    lx, ly, lz = _read_cryst1_box_angstrom(final_pdb)
    assert lx == pytest.approx(30.0, abs=1e-3)
    assert ly == pytest.approx(30.0, abs=1e-3)
    assert lz == pytest.approx(30.0, abs=1e-3)

    # gcmc_gpu-style "op_*" outputs should be created as additional final snapshots.
    op_pdb = work / "water_tip3p.pdb"
    op_top = work / "water_tip3p.top"
    assert op_pdb.exists()
    assert op_top.exists()
    assert op_top.read_text().strip(), "op_top unexpectedly empty"

    # Keep the comparison robust to header / formatting tweaks.
    assert _count_atom_records(op_pdb) == _count_atom_records(final_pdb)

    # Sanity: at least one SOL move is present in the structured acceptance log.
    assert accept_log.exists()
    _first_accept_record(accept_log, move="insertion", species="SOL")

    params = json.loads(params_json.read_text())
    inventory = _load_inp_key_inventory(Path(test_data_dir))
    _assert_inp_key_classification(params, example="water_tip3p", inventory=inventory)


def test_gcmc_gpu_version_gcmc_2_0_implies_angstrom_units_and_op_outputs(
    gcmc_cpu, test_data_dir, temp_dir
):
    templates = test_data_dir / "gcmc_gpu_examples"
    template = templates / "benzene_cbmc_cavity_gcmc_gpu_template.inp"
    assert template.exists()

    work = Path(temp_dir) / "gcmc_gpu_benzene"
    work.mkdir(parents=True, exist_ok=True)

    inp = _render_template(template, data_dir=test_data_dir, dst_dir=work)
    out_prefix = work / "out" / "gcmc"
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    accept_log = work / "out" / "acceptance.jsonl"
    params_json = work / "out" / "params.json"

    result = subprocess.run(
        [
            gcmc_cpu,
            "--inp",
            str(inp),
            "--prefix",
            str(out_prefix),
            "--dump-params",
            str(params_json),
            "--dump-accept",
            str(accept_log),
            "--store-probabilities",
        ],
        cwd=str(work),
        capture_output=True,
        text=True,
        timeout=120,
    )
    assert result.returncode == 0, result.stdout + result.stderr

    final_pdb = Path(f"{out_prefix}_final.pdb")
    assert final_pdb.exists()
    lx, ly, lz = _read_cryst1_box_angstrom(final_pdb)
    assert lx == pytest.approx(25.0, abs=1e-3)
    assert ly == pytest.approx(25.0, abs=1e-3)
    assert lz == pytest.approx(25.0, abs=1e-3)

    op_pdb = work / "benzene_box.pdb"
    op_top = work / "benzene_box.top"
    assert op_pdb.exists()
    assert op_top.exists()
    assert op_top.read_text().strip(), "op_top unexpectedly empty"
    assert _count_atom_records(op_pdb) == _count_atom_records(final_pdb)

    # CBMC trial count must come from "fragconf" (gcmc_gpu key).
    rec = _first_accept_record(accept_log, move="insertion", species="BENZ")
    assert int(rec["cbmcTrials"]) == 5

    params = json.loads(params_json.read_text())
    inventory = _load_inp_key_inventory(Path(test_data_dir))
    _assert_inp_key_classification(params, example="benzene_cbmc_cavity", inventory=inventory)


def test_gcmc_gpu_multi_salt_runs_and_converts_muex_list(gcmc_cpu, test_data_dir, temp_dir):
    templates = test_data_dir / "gcmc_gpu_examples"
    template = templates / "multi_salt_gcmc_gpu_template.inp"
    assert template.exists()

    work = Path(temp_dir) / "gcmc_gpu_multi_salt"
    work.mkdir(parents=True, exist_ok=True)

    inp = _render_template(template, data_dir=test_data_dir, dst_dir=work)
    out_prefix = work / "out" / "gcmc"
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    params_json = work / "out" / "params.json"

    result = subprocess.run(
        [
            gcmc_cpu,
            "--inp",
            str(inp),
            "--prefix",
            str(out_prefix),
            "--dump-params",
            str(params_json),
            "--store-probabilities",
        ],
        cwd=str(work),
        capture_output=True,
        text=True,
        timeout=120,
    )
    assert result.returncode == 0, result.stdout + result.stderr

    final_pdb = Path(f"{out_prefix}_final.pdb")
    assert final_pdb.exists()
    lx, ly, lz = _read_cryst1_box_angstrom(final_pdb)
    assert lx == pytest.approx(30.0, abs=1e-3)
    assert ly == pytest.approx(30.0, abs=1e-3)
    assert lz == pytest.approx(30.0, abs=1e-3)

    op_pdb = work / "multi_salt.pdb"
    op_top = work / "multi_salt.top"
    assert op_pdb.exists()
    assert op_top.exists()
    assert op_top.read_text().strip(), "op_top unexpectedly empty"
    assert _count_atom_records(op_pdb) == _count_atom_records(final_pdb)

    assert params_json.exists()
    params = json.loads(params_json.read_text())

    assert params["basic"]["inp_units"] == "gcmc_gpu"
    assert [float(x) for x in params["fragment"]["conc_list_M"]] == pytest.approx(
        [55.5, 0.15, 0.15], abs=1e-6
    )
    expected_muex_kj = [-5.60 * 4.184, -12.00 * 4.184, -11.50 * 4.184]
    assert [float(x) for x in params["fragment"]["muex_list_kj_mol"]] == pytest.approx(
        expected_muex_kj, abs=1e-5
    )
    inventory = _load_inp_key_inventory(Path(test_data_dir))
    _assert_inp_key_classification(params, example="multi_salt", inventory=inventory)


def test_gcmc_opencl_style_alias_keys_are_honored(gcmc_cpu, test_data_dir, temp_dir):
    templates = test_data_dir / "gcmc_gpu_examples"
    template = templates / "opencl_waterbox_style_gcmc_gpu_template.inp"
    assert template.exists()

    work = Path(temp_dir) / "gcmc_opencl_waterbox_style"
    work.mkdir(parents=True, exist_ok=True)

    inp = _render_template(template, data_dir=test_data_dir, dst_dir=work)
    out_prefix = work / "out" / "gcmc"
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    params_json = work / "out" / "params.json"

    result = subprocess.run(
        [
            gcmc_cpu,
            "--inp",
            str(inp),
            "--prefix",
            str(out_prefix),
            "--dump-params",
            str(params_json),
            "--store-probabilities",
        ],
        cwd=str(work),
        capture_output=True,
        text=True,
        timeout=120,
    )
    assert result.returncode == 0, result.stdout + result.stderr

    final_pdb = Path(f"{out_prefix}_final.pdb")
    assert final_pdb.exists()
    lx, ly, lz = _read_cryst1_box_angstrom(final_pdb)
    assert lx == pytest.approx(30.0, abs=1e-3)
    assert ly == pytest.approx(30.0, abs=1e-3)
    assert lz == pytest.approx(30.0, abs=1e-3)

    op_pdb = work / "opencl_waterbox_style.pdb"
    op_top = work / "opencl_waterbox_style.top"
    assert op_pdb.exists()
    assert op_top.exists()
    assert op_top.read_text().strip(), "op_top unexpectedly empty"
    assert _count_atom_records(op_pdb) == _count_atom_records(final_pdb)

    params = json.loads(params_json.read_text())
    inventory = _load_inp_key_inventory(Path(test_data_dir))
    _assert_inp_key_classification(params, example="opencl_waterbox_style", inventory=inventory)

    assert params_json.exists()
    params = json.loads(params_json.read_text())

    # energy_cutoff:8Å -> cutoff_nm:0.8nm under gcmc_gpu unit mode.
    assert float(params["space"]["cutoff_nm"]) == pytest.approx(0.8, abs=1e-6)
    assert float(params["energy"]["fragment_cutoff_nm"]) == pytest.approx(0.8, abs=1e-6)

    # use_vdw_radii_for_grid (plural alias) should map to SpaceInfo.use_vdw_radius_for_grid.
    assert params["space"]["use_vdw_radius_for_grid"] is True
    assert params["space"]["exclude_hydrogens_from_grid"] is False

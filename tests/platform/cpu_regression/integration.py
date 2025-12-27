"""
Integration tests for gcmc_cpu using real data files.

These tests are file-driven and validate behavior through output files / structured
dumps, not stdout/stderr log text.

External INP default semantics follow legacy gcmc_gpu/opencl conventions:
- lengths in Å
- chemical potentials in kcal/mol
"""

from __future__ import annotations

import json
import subprocess
import tempfile
from pathlib import Path

import pytest


GCMC_CPU_PATH = Path(__file__).resolve().parents[3] / "build" / "bin" / "gcmc_cpu"
TEST_DATA_DIR = Path(__file__).resolve().parents[2] / "data"


def _read_cryst1_box_angstrom(pdb_path: Path) -> tuple[float, float, float]:
    for line in pdb_path.read_text().splitlines():
        if line.startswith("CRYST1"):
            parts = line.split()
            assert len(parts) >= 4, f"Unexpected CRYST1 format: {line}"
            return float(parts[1]), float(parts[2]), float(parts[3])
    raise AssertionError(f"CRYST1 not found in {pdb_path}")


def _pdb_resnames(pdb_path: Path) -> set[str]:
    resnames: set[str] = set()
    for line in pdb_path.read_text().splitlines():
        if line.startswith(("ATOM", "HETATM")):
            res = line[17:20].strip().upper()
            if res:
                resnames.add(res)
    return resnames


def _sorted_pdb_atom_lines(pdb_path: Path) -> list[str]:
    atoms = []
    for line in pdb_path.read_text().splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        atoms.append(line[:6] + line[6:].rstrip())
    return sorted(atoms)


def _run_gcmc(inp: Path, *, cwd: Path, prefix: Path, extra_args: list[str] | None = None, timeout: int = 60):
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


class TestGCMCIntegration:
    """Integration tests demonstrating end-to-end CLI functionality."""

    @pytest.fixture
    def temp_dir(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)

    def test_fragitp_only_mode_inserts_and_writes_outputs(self, temp_dir: Path):
        sol_itp = TEST_DATA_DIR / "charmm36.ff" / "mol" / "sol.itp"
        assert sol_itp.exists()

        inp = temp_dir / "fragitp_only.inp"
        inp.write_text(
            f"""
version:gcmc_2.0
box:30.0 30.0 30.0
temperature:298.15

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

        prefix = temp_dir / "out" / "fragitp_only"
        prefix.parent.mkdir(parents=True, exist_ok=True)
        params_json = temp_dir / "out" / "params.json"
        accept_log = temp_dir / "out" / "accept.jsonl"

        _run_gcmc(
            inp,
            cwd=temp_dir,
            prefix=prefix,
            extra_args=["--seed", "42", "--dump-params", str(params_json), "--dump-accept", str(accept_log)],
            timeout=60,
        )

        final_pdb = Path(f"{prefix}_final.pdb")
        final_top = Path(f"{prefix}_final.top")
        assert final_pdb.exists()
        assert final_top.exists()
        assert _read_cryst1_box_angstrom(final_pdb) == pytest.approx((30.0, 30.0, 30.0), abs=1e-3)

        resnames = _pdb_resnames(final_pdb)
        assert "SOL" in resnames or "WAT" in resnames, f"Expected inserted water residue in output: {sorted(resnames)}"

        params = json.loads(params_json.read_text())
        assert params["basic"]["inp_units"] == "gcmc_gpu"
        assert params["basic"]["inp_units_explicit"] is False

        assert accept_log.exists()
        assert accept_log.read_text().strip(), "acceptance log unexpectedly empty"

    def test_multi_fragment_top_contains_all_species_rows(self, temp_dir: Path):
        sol_itp = TEST_DATA_DIR / "charmm36.ff" / "mol" / "sol.itp"
        na_itp = TEST_DATA_DIR / "charmm36.ff" / "mol" / "na.itp"
        cl_itp = TEST_DATA_DIR / "charmm36.ff" / "mol" / "cl.itp"
        assert sol_itp.exists()
        assert na_itp.exists()
        assert cl_itp.exists()

        inp = temp_dir / "salt_water.inp"
        inp.write_text(
            f"""
version:gcmc_2.0
box_size:40.0 40.0 40.0
temperature:298.15

fragitp:{sol_itp}
fragitp:{na_itp}
fragitp:{cl_itp}

fragname:SOL NA CL
fragconc:55.0 0.15 0.15
fragmuex:-10.5 -8.0 -7.5

mcsteps:0
nprint:1
""".strip()
            + "\n"
        )

        prefix = temp_dir / "out" / "salt"
        prefix.parent.mkdir(parents=True, exist_ok=True)
        params_json = temp_dir / "out" / "params.json"

        _run_gcmc(inp, cwd=temp_dir, prefix=prefix, extra_args=["--dump-params", str(params_json)], timeout=60)

        final_top = Path(f"{prefix}_final.top")
        assert final_top.exists()

        top_lines = final_top.read_text().splitlines()
        in_frag = False
        names: set[str] = set()
        for line in top_lines:
            if line.strip() == "[ fragments ]":
                in_frag = True
                continue
            if in_frag:
                if line.startswith("["):
                    break
                if not line.strip() or line.lstrip().startswith(";"):
                    continue
                cols = line.split()
                if cols:
                    names.add(cols[0].strip().upper())
        for want in ("SOL", "NA", "CL"):
            assert want in names, f"Missing {want} in [ fragments ]: {sorted(names)}"

        params = json.loads(params_json.read_text())
        assert [float(x) for x in params["fragment"]["conc_list_M"]] == pytest.approx([55.0, 0.15, 0.15], abs=1e-6)
        expected_muex = [-10.5 * 4.184, -8.0 * 4.184, -7.5 * 4.184]
        assert [float(x) for x in params["fragment"]["muex_list_kj_mol"]] == pytest.approx(expected_muex, abs=1e-4)

    @pytest.mark.parametrize("box_angstrom", [(20.0, 20.0, 20.0), (30.0, 40.0, 50.0)])
    def test_box_size_roundtrips_to_cryst1(self, temp_dir: Path, box_angstrom: tuple[float, float, float]):
        sol_itp = TEST_DATA_DIR / "charmm36.ff" / "mol" / "sol.itp"
        assert sol_itp.exists()

        inp = temp_dir / "box_test.inp"
        inp.write_text(
            f"""
version:gcmc_2.0
box_size:{box_angstrom[0]} {box_angstrom[1]} {box_angstrom[2]}
cutoff:10.0
temperature:298.15
fragitp:{sol_itp}
fragname:SOL
fragconc:55.0
fragmuex:1000.0
mcsteps:0
nprint:1
""".strip()
            + "\n"
        )

        prefix = temp_dir / "out" / "box"
        prefix.parent.mkdir(parents=True, exist_ok=True)

        _run_gcmc(inp, cwd=temp_dir, prefix=prefix, extra_args=["--seed", "123"], timeout=60)

        final_pdb = Path(f"{prefix}_final.pdb")
        assert final_pdb.exists()
        assert _read_cryst1_box_angstrom(final_pdb) == pytest.approx(box_angstrom, abs=1e-3)

    def test_deterministic_seed_compares_final_pdb_atom_records(self, temp_dir: Path):
        sol_itp = TEST_DATA_DIR / "charmm36.ff" / "mol" / "sol.itp"
        assert sol_itp.exists()

        inp = temp_dir / "deterministic.inp"
        inp.write_text(
            f"""
version:gcmc_2.0
box_size:25.0 25.0 25.0
temperature:298.15
fragitp:{sol_itp}
fragname:SOL
fragconc:30.0
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

        p1 = temp_dir / "out" / "det1"
        p2 = temp_dir / "out" / "det2"
        p3 = temp_dir / "out" / "det3"
        p1.parent.mkdir(parents=True, exist_ok=True)

        _run_gcmc(inp, cwd=temp_dir, prefix=p1, extra_args=["--seed", "999"], timeout=60)
        _run_gcmc(inp, cwd=temp_dir, prefix=p2, extra_args=["--seed", "999"], timeout=60)
        _run_gcmc(inp, cwd=temp_dir, prefix=p3, extra_args=["--seed", "1001"], timeout=60)

        a1 = _sorted_pdb_atom_lines(Path(f"{p1}_final.pdb"))
        a2 = _sorted_pdb_atom_lines(Path(f"{p2}_final.pdb"))
        a3 = _sorted_pdb_atom_lines(Path(f"{p3}_final.pdb"))

        assert a1, "No atom records captured"
        assert a1 == a2
        assert a1 != a3

    def test_long_run_writes_statistics_with_expected_steps(self, temp_dir: Path):
        sol_itp = TEST_DATA_DIR / "charmm36.ff" / "mol" / "sol.itp"
        assert sol_itp.exists()

        inp = temp_dir / "stability.inp"
        inp.write_text(
            f"""
version:gcmc_2.0
box_size:30.0 30.0 30.0
temperature:298.15
fragitp:{sol_itp}
fragname:SOL
fragconc:55.0
fragmuex:-10.5
mcsteps:200
nprint:50
""".strip()
            + "\n"
        )

        prefix = temp_dir / "out" / "stability"
        prefix.parent.mkdir(parents=True, exist_ok=True)

        _run_gcmc(inp, cwd=temp_dir, prefix=prefix, extra_args=["--seed", "42", "--print-freq", "50"], timeout=120)

        stats = Path(f"{prefix}_statistics.dat")
        assert stats.exists()
        steps = []
        for line in stats.read_text().splitlines():
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            steps.append(int(line.split()[0]))
        assert steps, "statistics file contained no data rows"
        assert steps[-1] == 200

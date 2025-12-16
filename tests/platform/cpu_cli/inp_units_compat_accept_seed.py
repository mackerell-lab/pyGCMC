"""
INP compatibility tests for gcmc_gpu-style units (acceptance consistency + seed behavior).
"""

from __future__ import annotations

import json
from pathlib import Path

from .inp_units_compat_helpers import _count_residues_by_resname, _run_gcmc_cpu, _write_inp


def test_dump_accept_log_consistent_with_final_pdb_counts(gcmc_cpu, test_data_dir, temp_dir):
    """
    No-cheating consistency check: the acceptance log 'accepted' flags must match the final system state.

    Run with a single-atom fragment (NA) and only insertion/deletion moves enabled.

    Verify: final molecule count (from *_final.pdb) == (#accepted insertions - #accepted deletions).
    """
    work = Path(temp_dir) / "accept_log_consistency"
    work.mkdir(parents=True, exist_ok=True)

    itp = test_data_dir / "charmm36.ff" / "mol" / "na.itp"
    assert itp.exists()

    out_prefix = work / "out" / "gcmc"
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    accept_log = work / "out" / "acceptance.jsonl"

    inp = work / "test.inp"
    _write_inp(
        inp,
        f"""
inp_units:gcmc_gpu
random_seed:999
fragitp:{itp}
fragname:NA
fragconc:55.0
fragmuex:0.0

box_size:10.0 10.0 10.0
cutoff:4.0
temperature:300.0
moves_per_step:1
mcsteps:20
nprint:1000
mc_move_prob:1 1 0 0
""",
    )

    result = _run_gcmc_cpu(
        gcmc_cpu,
        workdir=work,
        inp=inp,
        out_prefix=out_prefix,
        extra_args=["--dump-accept", str(accept_log)],
        timeout=30,
    )
    assert result.returncode == 0, result.stdout + result.stderr

    records = [json.loads(line) for line in accept_log.read_text().splitlines() if line.strip()]
    na_records = [r for r in records if str(r.get("species", "")).strip().upper() == "NA"]
    assert na_records

    assert any(r.get("move") == "insertion" for r in na_records)
    assert any(r.get("move") == "deletion" for r in na_records)
    assert any(bool(r.get("accepted")) for r in na_records if r.get("move") == "insertion")

    accepted_insert = sum(1 for r in na_records if r.get("move") == "insertion" and bool(r.get("accepted")))
    accepted_delete = sum(1 for r in na_records if r.get("move") == "deletion" and bool(r.get("accepted")))
    expected_final = accepted_insert - accepted_delete
    assert expected_final >= 0

    out_pdb = Path(f"{out_prefix}_final.pdb")
    assert out_pdb.exists()
    final_count = _count_residues_by_resname(out_pdb, "NA")

    assert final_count == expected_final


def test_inp_random_seed_used_when_cli_missing(gcmc_cpu, test_data_dir, temp_dir):
    """random_seed/seed in INP should make runs reproducible when --seed is not provided."""
    itp = test_data_dir / "charmm36.ff" / "mol" / "sol.itp"
    assert itp.exists()

    def run_once(work: Path, seed: int) -> tuple[tuple[str, str, int, float, float, float], ...]:
        work.mkdir(parents=True, exist_ok=True)
        out_prefix = work / "out" / "gcmc"
        out_prefix.parent.mkdir(parents=True, exist_ok=True)

        inp = work / "test.inp"
        _write_inp(
            inp,
            f"""
inp_units:gcmc_gpu
random_seed:{seed}
fragitp:{itp}
box_size:30.0 30.0 30.0
cutoff:12.0
temperature:300.0
moves_per_step:1
mcsteps:10
nprint:5
fragname:SOL
fragconc:55.0
fragmuex:50.0
attempt_prob_ins:1.0
attempt_prob_del:0.0
attempt_prob_trn:0.0
attempt_prob_rot:0.0
""",
        )

        result = _run_gcmc_cpu(
            gcmc_cpu,
            workdir=work,
            inp=inp,
            out_prefix=out_prefix,
        )
        assert result.returncode == 0, result.stdout + result.stderr

        out_pdb = Path(f"{out_prefix}_final.pdb")
        assert out_pdb.exists()
        atoms: list[tuple[str, str, int, float, float, float]] = []
        for line in out_pdb.read_text().splitlines():
            if not line.startswith(("ATOM", "HETATM")):
                continue
            atom_name = line[12:16].strip().upper()
            resname = line[17:20].strip().upper()
            resid = int(line[22:26])
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            atoms.append((resname, atom_name, resid, round(x, 3), round(y, 3), round(z, 3)))
        return tuple(sorted(atoms))

    run1 = Path(temp_dir) / "seed_run1"
    run2 = Path(temp_dir) / "seed_run2"
    assert run_once(run1, 123) == run_once(run2, 123)

    # And with a different seed, we should almost certainly get a different trajectory/output.
    run3 = Path(temp_dir) / "seed_run3"
    run4 = Path(temp_dir) / "seed_run4"
    assert run_once(run3, 123) != run_once(run4, 124)


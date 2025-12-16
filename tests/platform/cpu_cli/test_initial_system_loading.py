"""
Regression test: initial PDB+TOP system should populate MCState even when PAR files are GROMACS ITP.

This guards against the historical fallback-to-empty-state behavior when result.forceField is null.
"""

from __future__ import annotations

import json
import math
import subprocess
from pathlib import Path

import pytest


def test_builder_initial_system_loaded_with_itp_par(gcmc_cpu, temp_dir):
    work = Path(temp_dir) / "init_system"
    work.mkdir(parents=True, exist_ok=True)

    # Minimal PDB (Å) with one atom and CRYST1 box
    pdb = work / "sys.pdb"
    pdb.write_text(
        "\n".join(
            [
                "CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1           1",
                "ATOM      1  C   MOL A   1       0.000   0.000   0.000  1.00  0.00           C",
                "END",
                "",
            ]
        )
    )

    # Minimal GROMACS TOP with one atom matching the PDB
    top = work / "sys.top"
    top.write_text(
        "\n".join(
            [
                "[ defaults ]",
                "1 2 yes 0.5 0.8333",
                "",
                "[ atomtypes ]",
                "; name  at.num  mass     charge   ptype    sigma      epsilon",
                "C       6       12.011   0.000    A        3.00000e-01  1.00000e-01",
                "",
                "[ moleculetype ]",
                "MOL  2",
                "",
                "[ atoms ]",
                "; nr  type  resnr  residue  atom  cgnr  charge    mass",
                "1   C     1      MOL      C     1     0.000   12.011",
                "",
                "[ system ]",
                "Minimal",
                "",
                "[ molecules ]",
                "MOL  1",
                "",
            ]
        )
    )

    # Minimal nonbonded ITP for ItpNonbondedParser (sigma/epsilon are nm/kJ)
    par = work / "ffnonbonded.itp"
    par.write_text(
        "\n".join(
            [
                "[ atomtypes ]",
                "; name  at.num  mass     charge   ptype    sigma      epsilon",
                "C       6       12.011   0.000    A        3.00000e-01  1.00000e-01",
                "",
            ]
        )
    )

    out_prefix = work / "out" / "gcmc"
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    inp = work / "run.inp"
    inp.write_text(
        f"""
inp_units:gcmc_gpu
par:{par}
top:{top}
pdb:{pdb}

box_size:10.0 10.0 10.0
gc_center:15.0 15.0 15.0
cutoff:12.0
temperature:300
moves_per_step:1
mcsteps:1
nprint:1

fragname:water
fragmuex:0.0
mc_move_prob:1 0 0 0
""".strip()
        + "\n"
    )

    result = subprocess.run(
        [gcmc_cpu, "--inp", str(inp), "--prefix", str(out_prefix), "--verbose"],
        cwd=str(work),
        capture_output=True,
        text=True,
        timeout=20,
    )
    assert result.returncode == 0, result.stdout + result.stderr

    out_pdb = Path(f"{out_prefix}_final.pdb")
    assert out_pdb.exists()

    # The initial system must survive into output (i.e., we did not silently fall back to an empty state).
    resnames = set()
    for line in out_pdb.read_text().splitlines():
        if line.startswith(("ATOM", "HETATM")):
            resname = line[17:20].strip().upper()
            if resname:
                resnames.add(resname)
    assert "MOL" in resnames, f"Expected initial MOL residue in output: {sorted(resnames)}"


def test_insertion_does_not_overwrite_initial_residues(gcmc_cpu, test_data_dir, temp_dir):
    """
    When an initial PDB+TOP system is loaded, fragment insertions must not reuse residue indices
    starting at 0, otherwise the first insertion overwrites the protein/system residues.
    """
    work = Path(temp_dir) / "init_system_insert"
    work.mkdir(parents=True, exist_ok=True)

    # Initial system: 1 residue MOL with 1 atom C (Å) in a large box.
    pdb = work / "sys.pdb"
    pdb.write_text(
        "\n".join(
            [
                "CRYST1  100.000  100.000  100.000  90.00  90.00  90.00 P 1           1",
                "ATOM      1  C   MOL A   1       0.000   0.000   0.000  1.00  0.00           C",
                "END",
                "",
            ]
        )
    )

    top = work / "sys.top"
    top.write_text(
        "\n".join(
            [
                "[ defaults ]",
                "1 2 yes 0.5 0.8333",
                "",
                "[ moleculetype ]",
                "MOL  2",
                "",
                "[ atoms ]",
                "; nr  type  resnr  residue  atom  cgnr  charge    mass",
                "1   C     1      MOL      C     1     0.000   12.011",
                "",
                "[ system ]",
                "Minimal",
                "",
                "[ molecules ]",
                "MOL  1",
                "",
            ]
        )
    )

    # Minimal ITP nonbonded (sigma/epsilon are nm/kJ) covering both initial system and SOL fragment.
    par = work / "ffnonbonded.itp"
    par.write_text(
        "\n".join(
            [
                "[ atomtypes ]",
                "; name  at.num  mass     charge   ptype    sigma      epsilon",
                "C       6       12.011   0.000    A        3.00000e-01  0.00000e+00",
                "OT      8       15.9994  0.000    A        3.15000e-01  0.00000e+00",
                "HT      1       1.0080   0.000    A        0.00000e+00  0.00000e+00",
                "",
            ]
        )
    )

    sol_itp = test_data_dir / "charmm36.ff" / "mol" / "sol.itp"
    assert sol_itp.exists()

    out_prefix = work / "out" / "gcmc"
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    inp = work / "run.inp"
    inp.write_text(
        f"""
inp_units:gcmc_gpu
par:{par}
top:{top}
pdb:{pdb}
fragitp:{sol_itp}

cutoff:12.0
temperature:300
moves_per_step:1
mcsteps:1
nprint:1

fragname:SOL
fragmuex:1000.0
attempt_prob_ins:1.0
attempt_prob_del:0.0
attempt_prob_trn:0.0
attempt_prob_rot:0.0
""".strip()
        + "\n"
    )

    result = subprocess.run(
        [gcmc_cpu, "--inp", str(inp), "--prefix", str(out_prefix), "--verbose"],
        cwd=str(work),
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert result.returncode == 0, result.stdout + result.stderr

    out_pdb = Path(f"{out_prefix}_final.pdb")
    assert out_pdb.exists()

    resnames = set()
    for line in out_pdb.read_text().splitlines():
        if line.startswith("ATOM") or line.startswith("HETATM"):
            resname = line[17:20].strip().upper()
            if resname:
                resnames.add(resname)

    assert "MOL" in resnames, f"Initial residue lost/overwritten: {sorted(resnames)}"
    assert "SOL" in resnames, f"Insertion did not create SOL residue: {sorted(resnames)}"


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


def test_cavity_bias_fraction_matches_simple_grid_occupancy(gcmc_cpu, test_data_dir, temp_dir):
    """
    End-to-end cavity bias sanity check (no log matching):
    - 1 nm^3 box, grid spacing 0.5 nm -> 2x2x2 = 8 grid points
    - one atom at origin occupies exactly 1 grid point when totalRadius=0.3 nm
    => cavityFraction = 7/8 = 0.875
    """
    work = Path(temp_dir) / "cavity_fraction_simple"
    work.mkdir(parents=True, exist_ok=True)

    pdb = work / "sys.pdb"
    pdb.write_text(
        "\n".join(
            [
                "CRYST1   10.000   10.000   10.000  90.00  90.00  90.00 P 1           1",
                "ATOM      1  C   MOL A   1       0.000   0.000   0.000  1.00  0.00           C",
                "END",
                "",
            ]
        )
    )

    top = work / "sys.top"
    top.write_text(
        "\n".join(
            [
                "[ defaults ]",
                "1 2 yes 0.5 0.8333",
                "",
                "[ moleculetype ]",
                "MOL  2",
                "",
                "[ atoms ]",
                "; nr  type  resnr  residue  atom  cgnr  charge    mass",
                "1   C     1      MOL      C     1     0.000   12.011",
                "",
                "[ system ]",
                "Minimal",
                "",
                "[ molecules ]",
                "MOL  1",
                "",
            ]
        )
    )

    par = work / "ffnonbonded.itp"
    par.write_text(
        "\n".join(
            [
                "[ atomtypes ]",
                "; name  at.num  mass     charge   ptype    sigma      epsilon",
                "C       6       12.011   0.000    A        4.00000e-01  5.00000e+00",
                "",
            ]
        )
    )

    frag_itp = test_data_dir / "charmm36.ff" / "mol" / "na.itp"
    assert frag_itp.exists()

    out_prefix = work / "out" / "gcmc"
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    accept_log = work / "out" / "acceptance.jsonl"

    inp = work / "run.inp"
    inp.write_text(
        f"""
inp_units:gcmc_gpu
random_seed:123
par:{par}
top:{top}
pdb:{pdb}
fragitp:{frag_itp}

box_size:10.0 10.0 10.0
cutoff:12.0
temperature:300
moves_per_step:1
mcsteps:1
nprint:1

use_cavity_bias:yes
grid_dx:5.0
fragradius:1.0

fragname:NA
fragconc:0.01
fragmuex:0.0
mc_move_prob:1 0 0 0
""".strip()
        + "\n"
    )

    result = subprocess.run(
        [gcmc_cpu, "--inp", str(inp), "--prefix", str(out_prefix), "--dump-accept", str(accept_log)],
        cwd=str(work),
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert result.returncode == 0, result.stdout + result.stderr

    rec = _first_accept_record(accept_log, move="insertion", species="NA")
    assert float(rec["cavityFraction"]) == pytest.approx(0.875, abs=1e-12)
    assert float(rec["vEff"]) == pytest.approx(float(rec["vBox"]) * 0.875, abs=1e-12)


def test_cbmc_and_cavity_terms_reconstruct_insertion_pacc(gcmc_cpu, test_data_dir, temp_dir):
    """
    Ensure CBMC + cavity terms enter the insertion acceptance probability (pAcc) consistently.
    Reconstruct pAcc from JSONL fields (no stdout/stderr dependency).
    """
    work = Path(temp_dir) / "cbmc_cavity_pacc"
    work.mkdir(parents=True, exist_ok=True)

    pdb = work / "sys.pdb"
    pdb.write_text(
        "\n".join(
            [
                "CRYST1   10.000   10.000   10.000  90.00  90.00  90.00 P 1           1",
                "ATOM      1  C   MOL A   1       0.000   0.000   0.000  1.00  0.00           C",
                "END",
                "",
            ]
        )
    )

    top = work / "sys.top"
    top.write_text(
        "\n".join(
            [
                "[ defaults ]",
                "1 2 yes 0.5 0.8333",
                "",
                "[ moleculetype ]",
                "MOL  2",
                "",
                "[ atoms ]",
                "; nr  type  resnr  residue  atom  cgnr  charge    mass",
                "1   C     1      MOL      C     1     0.000   12.011",
                "",
                "[ system ]",
                "Minimal",
                "",
                "[ molecules ]",
                "MOL  1",
                "",
            ]
        )
    )

    par = work / "ffnonbonded.itp"
    par.write_text(
        "\n".join(
            [
                "[ atomtypes ]",
                "; name  at.num  mass     charge   ptype    sigma      epsilon",
                "C       6       12.011   0.000    A        4.00000e-01  5.00000e+00",
                "",
            ]
        )
    )

    frag_itp = test_data_dir / "charmm36.ff" / "mol" / "na.itp"
    assert frag_itp.exists()

    out_prefix = work / "out" / "gcmc"
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    accept_log = work / "out" / "acceptance.jsonl"

    inp = work / "run.inp"
    inp.write_text(
        f"""
inp_units:gcmc_gpu
random_seed:777
par:{par}
top:{top}
pdb:{pdb}
fragitp:{frag_itp}

box_size:10.0 10.0 10.0
cutoff:12.0
temperature:300
moves_per_step:1
mcsteps:1
nprint:1

use_cavity_bias:yes
grid_dx:5.0
fragradius:1.0

use_conf_bias:yes
fragconf:7

fragname:NA
fragconc:0.01
fragmuex:0.0
mc_move_prob:1 0 0 0
""".strip()
        + "\n"
    )

    result = subprocess.run(
        [gcmc_cpu, "--inp", str(inp), "--prefix", str(out_prefix), "--dump-accept", str(accept_log)],
        cwd=str(work),
        capture_output=True,
        text=True,
        timeout=60,
    )
    assert result.returncode == 0, result.stdout + result.stderr

    rec = _first_accept_record(accept_log, move="insertion", species="NA")

    cbmc_trials = int(rec["cbmcTrials"])
    assert cbmc_trials == 7

    cavity_fraction = float(rec["cavityFraction"])
    assert cavity_fraction == pytest.approx(0.875, abs=1e-12)

    rosen = float(rec["rosenbluthWeight"])
    assert 0.0 < rosen <= 1.0
    assert rosen < 0.999999

    beta = 1.0 / (8.314e-3 * 300.0)
    delta_u = float(rec["deltaU"])
    z = float(rec["z"])  # activity (1/nm^3)
    volume = float(rec["vBox"])
    n_before = int(rec["nBefore"])

    proposal_ratio = float(rec.get("proposalRatio", 1.0))
    proposal_log_ratio = -math.log(proposal_ratio) if proposal_ratio > 0.0 else 0.0

    log_ratio = (
        proposal_log_ratio
        - beta * delta_u
        + math.log(z)
        + math.log(volume)
        + math.log(cavity_fraction)
        - math.log(n_before + 1)
        + math.log(rosen)
    )
    expected_pacc = 1.0 if log_ratio >= 0.0 else math.exp(log_ratio)

    assert float(rec["pAcc"]) == pytest.approx(expected_pacc, rel=1e-12, abs=1e-12)

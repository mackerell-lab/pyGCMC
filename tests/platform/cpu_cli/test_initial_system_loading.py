"""
Regression test: initial PDB+TOP system should populate MCState even when PAR files are GROMACS ITP.

This guards against the historical fallback-to-empty-state behavior when result.forceField is null.
"""

from __future__ import annotations

import subprocess
from pathlib import Path


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

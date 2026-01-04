"""
Drude fragment GCMC contracts (CLI, file-driven).

Goal: ensure gcmc_cpu correctly handles Drude fragments provided via fragitp/monomerdir:
- Drude oscillators are included in trial energies (no silent omission).
- INSERT/DELETE ΔU closes against analytic Coulomb + Drude spring for a minimal system.

This deliberately avoids stdout/stderr text matching and uses:
- --dump-accept JSONL for ΔU
- *_final.pdb and *_statistics.dat for geometry/energy closure
"""

from __future__ import annotations

import json
import math
import subprocess
from pathlib import Path

import pytest


COULOMB_KJ_NM_PER_MOL_E2 = 138.935456  # must match DrudeConstants.ONE_4PI_EPS0


def _write_text(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content.strip() + "\n")


def _read_stats_energy_kj(stats_path: Path) -> float:
    lines = [ln for ln in stats_path.read_text().splitlines() if ln and not ln.startswith("#")]
    assert lines, f"No data lines in statistics file: {stats_path}"
    last = lines[-1].split()
    assert len(last) >= 2, f"Malformed statistics line: {lines[-1]}"
    return float(last[1])


def _read_cryst1_box_nm(pdb_path: Path) -> tuple[float, float, float]:
    for line in pdb_path.read_text().splitlines():
        if line.startswith("CRYST1"):
            # Å -> nm
            a = float(line[6:15]) * 0.1
            b = float(line[15:24]) * 0.1
            c = float(line[24:33]) * 0.1
            return a, b, c
    raise AssertionError(f"CRYST1 not found in {pdb_path}")


def _pdb_atom_xyz_angstrom(
    pdb_path: Path, *, resname: str, atom_name: str, resid: int | None = None
) -> tuple[float, float, float]:
    want_res = resname.strip().upper()
    want_atom = atom_name.strip().upper()
    for line in pdb_path.read_text().splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        if line[17:20].strip().upper() != want_res:
            continue
        if line[12:16].strip().upper() != want_atom:
            continue
        if resid is not None and int(line[22:26]) != resid:
            continue
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        return x, y, z
    suffix = f"{resname}{resid if resid is not None else ''} {atom_name}".strip()
    raise AssertionError(f"Atom not found in {pdb_path}: {suffix}")


def _minimum_image_dist_nm(
    a_a: tuple[float, float, float],
    b_a: tuple[float, float, float],
    box_nm: tuple[float, float, float],
) -> float:
    dx = (a_a[0] - b_a[0]) * 0.1
    dy = (a_a[1] - b_a[1]) * 0.1
    dz = (a_a[2] - b_a[2]) * 0.1
    bx, by, bz = box_nm
    if bx > 0:
        dx -= bx * round(dx / bx)
    if by > 0:
        dy -= by * round(dy / by)
    if bz > 0:
        dz -= bz * round(dz / bz)
    return math.sqrt(dx * dx + dy * dy + dz * dz)


def _load_accept_records(path: Path) -> list[dict]:
    return [json.loads(line) for line in path.read_text().splitlines() if line.strip()]


def _first_record(path: Path, *, move: str, species: str) -> dict:
    want_move = move.strip().lower()
    want_species = species.strip().upper()
    for rec in _load_accept_records(path):
        if str(rec.get("move", "")).strip().lower() != want_move:
            continue
        if str(rec.get("species", "")).strip().upper() != want_species:
            continue
        return rec
    raise AssertionError(f"No {move}/{species} record found in {path}")


def _setup_drude_fragment_system(work: Path) -> tuple[Path, Path, Path, Path, Path]:
    """
    Minimal system with no Drude in the initial topology:
      FIX: single +1 external point charge
      RES: insertable Drude fragment (parent + Drude)
    """
    psf = work / "system.psf"
    pdb = work / "system.pdb"
    prm = work / "ff.prm"
    frag_itp = work / "res.itp"
    monomerdir = work / "monomerdir"
    monomerdir.mkdir(parents=True, exist_ok=True)
    frag_pdb = monomerdir / "RES.pdb"

    parent_type = "P"
    drude_type = "D"
    ext_type = "Q"

    alpha_a3 = 1.0  # Å^3 -> 0.001 nm^3
    thole = 2.6
    drude_dx_a = 0.005  # Å, tiny non-zero initial offset

    _write_text(
        psf,
        f"""
PSF EXT

       1 !NTITLE
REMARKS minimal Drude fragment insert/delete regression

       1 !NATOM
       1 SYS  1 FIX  Q1   {ext_type}      1.000000  12.0110  0   0.000  0.000

       0 !NBOND: bonds
""",
    )

    cx, cy, cz = 25.0, 25.0, 25.0
    _write_text(
        pdb,
        f"""
CRYST1   50.000   50.000   50.000  90.00  90.00  90.00 P 1           1
ATOM      1  Q1  FIX A   1      {cx:8.3f}{cy:8.3f}{cz:8.3f}  1.00  0.00           Q
END
""",
    )

    _write_text(
        frag_pdb,
        f"""
ATOM      1  P1  RES A   1       0.000   0.000   0.000  1.00  0.00           C
ATOM      2 DP1  RES A   1       {drude_dx_a:6.3f}   0.000   0.000  1.00  0.00           D
END
""",
    )

    _write_text(
        prm,
        f"""
* minimal PRM for Drude fragment regression
*

NONBONDED nbxmod 5 atom cdiel -
cutnb 999.0 ctonnb 0.0 ctofnb 0.0 eps 1.0 e14fac 1.0 wmin 1.5
{parent_type:8s}  0   0.0   0.0
{drude_type:8s}   0   0.0   0.0
{ext_type:8s}     0   0.0   0.0

ALPHA
{parent_type:8s}  {alpha_a3:.3f}  {thole:.3f}
END
""",
    )

    _write_text(
        frag_itp,
        f"""
[ moleculetype ]
RES     3

[ atoms ]
; nr  type   resnr  residue  atom  cgnr  charge    mass
1   {parent_type}  1   RES  P1   1   1.000000  12.0110
2   {drude_type}   1   RES  DP1  1  -1.000000   0.4000

[ bonds ]
1  2  1
""",
    )

    return psf, pdb, prm, frag_itp, monomerdir


def test_drude_fragment_insertion_deltaU_closes_to_coulomb_plus_spring(gcmc_cpu, temp_dir):
    work = Path(temp_dir) / "drude_fragment" / "insert_energy_closure"
    work.mkdir(parents=True, exist_ok=True)

    psf, pdb, prm, frag_itp, monomerdir = _setup_drude_fragment_system(work)

    out_prefix = work / "out" / "gcmc"
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    accept_log = work / "out" / "acceptance.jsonl"

    inp = work / "run.inp"
    _write_text(
        inp,
        f"""
energy_method:direct
top:{psf}
pdb:{pdb}
par:{prm}
monomerdir:{monomerdir}
fragitp:{frag_itp}
fragname:RES
fragconc:55.0
fragmuex:0.0

box_size:50.0 50.0 50.0
cutoff:12.0
temperature:300.0
use_cavity_bias:no

moves_per_step:1
mcsteps:1
nprint:1
mc_move_prob:1 0 0 0
""",
    )

    result = subprocess.run(
        [
            gcmc_cpu,
            "--inp",
            str(inp),
            "--prefix",
            str(out_prefix),
            "--seed",
            "7",
            "--stats-interval",
            "1",
            "--dump-accept",
            str(accept_log),
        ],
        cwd=str(work),
        capture_output=True,
        text=True,
        timeout=60,
    )
    assert result.returncode == 0, result.stdout + result.stderr

    out_pdb = Path(f"{out_prefix}_final.pdb")
    stats = Path(f"{out_prefix}_statistics.dat")
    assert out_pdb.exists()
    assert stats.exists()
    assert accept_log.exists()

    rec = _first_record(accept_log, move="insertion", species="RES")
    assert bool(rec.get("accepted", False)) is True
    delta_u = float(rec["deltaU"])

    box_nm = _read_cryst1_box_nm(out_pdb)
    parent = _pdb_atom_xyz_angstrom(out_pdb, resname="RES", atom_name="P1")
    drude = _pdb_atom_xyz_angstrom(out_pdb, resname="RES", atom_name="DP1")
    ext = _pdb_atom_xyz_angstrom(out_pdb, resname="FIX", atom_name="Q1")

    r_pe = _minimum_image_dist_nm(parent, ext, box_nm)
    r_de = _minimum_image_dist_nm(drude, ext, box_nm)
    disp_pd = _minimum_image_dist_nm(parent, drude, box_nm)

    q_parent, q_drude, q_ext = 1.0, -1.0, 1.0
    alpha_nm3 = 0.001
    cutoff_nm = 12.0 * 0.1
    coulomb = 0.0
    if r_pe <= cutoff_nm:
        coulomb += COULOMB_KJ_NM_PER_MOL_E2 * q_parent * q_ext / r_pe
    if r_de <= cutoff_nm:
        coulomb += COULOMB_KJ_NM_PER_MOL_E2 * q_drude * q_ext / r_de
    k_spring = COULOMB_KJ_NM_PER_MOL_E2 * (q_drude * q_drude) / alpha_nm3
    spring = 0.5 * k_spring * disp_pd * disp_pd
    expected = coulomb + spring

    energy_stats = _read_stats_energy_kj(stats)
    assert energy_stats == pytest.approx(expected, abs=5e-3, rel=1e-6)
    assert delta_u == pytest.approx(expected, abs=5e-3, rel=1e-6)


def test_drude_fragment_delete_deltaU_is_negative_of_insert(gcmc_cpu, temp_dir):
    work = Path(temp_dir) / "drude_fragment" / "insert_delete_symmetry"
    work.mkdir(parents=True, exist_ok=True)

    psf, pdb, prm, frag_itp, monomerdir = _setup_drude_fragment_system(work)

    out_prefix = work / "out" / "gcmc"
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    ins_accept = work / "out" / "ins_accept.jsonl"

    inp_ins = work / "insert.inp"
    _write_text(
        inp_ins,
        f"""
energy_method:direct
top:{psf}
pdb:{pdb}
par:{prm}
monomerdir:{monomerdir}
fragitp:{frag_itp}
fragname:RES
fragconc:55.0
fragmuex:0.0

box_size:50.0 50.0 50.0
cutoff:12.0
temperature:300.0
use_cavity_bias:no

moves_per_step:1
mcsteps:2
nprint:1
mc_move_prob:1 0 0 0
""",
    )

    res_ins = subprocess.run(
        [
            gcmc_cpu,
            "--inp",
            str(inp_ins),
            "--prefix",
            str(out_prefix),
            "--seed",
            "11",
            "--dump-accept",
            str(ins_accept),
            "--checkpoint-freq",
            "1",
            "--max-molecules-per-type",
            "1",
        ],
        cwd=str(work),
        capture_output=True,
        text=True,
        timeout=60,
    )
    assert res_ins.returncode == 0, res_ins.stdout + res_ins.stderr

    ins_rec = _first_record(ins_accept, move="insertion", species="RES")
    assert bool(ins_rec.get("accepted", False)) is True
    delta_u_ins = float(ins_rec["deltaU"])
    assert abs(delta_u_ins) > 1e-3

    checkpoint = Path(f"{out_prefix}_checkpoint_1.dat")
    assert checkpoint.exists()

    del_accept = work / "out" / "del_accept.jsonl"
    inp_del = work / "delete.inp"
    _write_text(
        inp_del,
        f"""
energy_method:direct
top:{psf}
pdb:{pdb}
par:{prm}
monomerdir:{monomerdir}
fragitp:{frag_itp}
fragname:RES
fragconc:55.0
fragmuex:0.0

box_size:50.0 50.0 50.0
cutoff:12.0
temperature:300.0
use_cavity_bias:no

moves_per_step:1
mcsteps:1
nprint:1
mc_move_prob:0 1 0 0
""",
    )

    res_del = subprocess.run(
        [
            gcmc_cpu,
            "--inp",
            str(inp_del),
            "--prefix",
            str(out_prefix),
            "--seed",
            "13",
            "--resume",
            str(checkpoint),
            "--dump-accept",
            str(del_accept),
        ],
        cwd=str(work),
        capture_output=True,
        text=True,
        timeout=60,
    )
    assert res_del.returncode == 0, res_del.stdout + res_del.stderr
    del_rec = _first_record(del_accept, move="deletion", species="RES")
    delta_u_del = float(del_rec["deltaU"])

    assert delta_u_del == pytest.approx(-delta_u_ins, rel=1e-4, abs=5e-3)

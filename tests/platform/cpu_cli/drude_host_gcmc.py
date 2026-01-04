"""
Drude fragment GCMC in a Drude host (CLI, file-driven).

Goal: ensure gcmc_cpu can insert/delete a Drude fragment when the initial host system
already contains Drude particles (e.g., Drude water / Drude protein environments).

Strictness:
- No stdout/stderr log-string assertions.
- Use --dump-accept JSONL for per-move ΔU.
- Use *_final.pdb for geometry and analytic closure (Coulomb + Drude spring).

Design:
- Build a minimal Drude "host" (one polarizable parent + one Drude) near a fixed charge QH.
- Force insertions to occur far away (gcmc_region box), near a separate fixed charge QI.
- With direct+cutoff, host↔insert interactions are beyond cutoff, so ΔU closes to the same
  analytic expression as the vacuum-host case, while still exercising "hostDrudeParticles != 0".
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


def _pdb_contains_resname(pdb_path: Path, resname: str) -> bool:
    want = resname.strip().upper()
    for line in pdb_path.read_text().splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        if line[17:20].strip().upper() == want:
            return True
    return False


def _setup_drude_host_and_insert_region(work: Path) -> tuple[Path, Path, Path, Path, Path]:
    """
    Host system (Drude-enabled):
      HST: one polarizable parent + one Drude
      FIX: QH (fixed +1) near host to induce host polarization
      FAR: QI (fixed +1) near insertion region to induce inserted-fragment polarization

    Insertable fragment:
      RES: parent + Drude (from fragitp/monomerdir)
    """
    psf = work / "system.psf"
    pdb = work / "system.pdb"
    prm = work / "ff.prm"
    frag_itp = work / "res.itp"
    monomerdir = work / "monomerdir"
    monomerdir.mkdir(parents=True, exist_ok=True)
    frag_pdb = monomerdir / "RES.pdb"

    host_parent_type = "O"
    frag_parent_type = "P"
    drude_type = "D"
    q_type = "Q"

    alpha_a3 = 1.0  # Å^3 -> 0.001 nm^3
    thole = 2.6
    drude_dx_a = 0.005  # Å

    _write_text(
        psf,
        f"""
PSF EXT DRUDE

       1 !NTITLE
REMARKS drude host + drude fragment insertion/deletion regression

       4 !NATOM
       1 SYS  1 HST  O1   {host_parent_type}   1.000000  12.0110  0  {thole:.3f}  {alpha_a3:.3f}
       2 SYS  1 HST  D1   {drude_type}        -1.000000   0.4000  0  0.000  0.000
       3 SYS  2 FIX  QH   {q_type}             1.000000  12.0110  0  0.000  0.000
       4 SYS  3 FAR  QI   {q_type}             1.000000  12.0110  0  0.000  0.000

       1 !NBOND: bonds
       1       2
""",
    )

    # Host: centered around (10,10,10). Insertion region is far away around ~ (38.5,38.5,38.5).
    host_c = (10.0, 10.0, 10.0)
    qhost = (13.0, 10.0, 10.0)  # Å (near host)
    # Place QI at least 3 Å away from the insertion box, while keeping all insertion points within cutoff.
    qins = (45.0, 38.5, 38.5)  # Å (near insertion region)

    _write_text(
        pdb,
        f"""
CRYST1   50.000   50.000   50.000  90.00  90.00  90.00 P 1           1
ATOM      1  O1  HST A   1    {host_c[0]:8.3f}{host_c[1]:8.3f}{host_c[2]:8.3f}  1.00  0.00           O
ATOM      2  D1  HST A   1    {host_c[0] + drude_dx_a:8.3f}{host_c[1]:8.3f}{host_c[2]:8.3f}  1.00  0.00           D
ATOM      3  QH  FIX A   2    {qhost[0]:8.3f}{qhost[1]:8.3f}{qhost[2]:8.3f}  1.00  0.00           Q
ATOM      4  QI  FAR A   3    {qins[0]:8.3f}{qins[1]:8.3f}{qins[2]:8.3f}  1.00  0.00           Q
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

    # Minimal PRM: LJ is zero for all types; alpha/thole provided for fragment parent type P.
    _write_text(
        prm,
        f"""
* minimal PRM for Drude host GCMC regression
*

NONBONDED nbxmod 5 atom cdiel -
cutnb 999.0 ctonnb 0.0 ctofnb 0.0 eps 1.0 e14fac 1.0 wmin 1.5
{host_parent_type:8s}  0   0.0   0.0
{frag_parent_type:8s}  0   0.0   0.0
{drude_type:8s}        0   0.0   0.0
{q_type:8s}            0   0.0   0.0

ALPHA
{frag_parent_type:8s}  {alpha_a3:.3f}  {thole:.3f}
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
1   {frag_parent_type}  1   RES  P1   1   1.000000  12.0110
2   {drude_type}        1   RES  DP1  1  -1.000000   0.4000

[ bonds ]
1  2  1
""",
    )

    return psf, pdb, prm, frag_itp, monomerdir


def test_drude_fragment_insertion_deltaU_closes_with_drude_host_present(gcmc_cpu, temp_dir):
    work = Path(temp_dir) / "drude_host_gcmc" / "insert_closure"
    work.mkdir(parents=True, exist_ok=True)

    psf, pdb, prm, frag_itp, monomerdir = _setup_drude_host_and_insert_region(work)

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

gcmc_region:box 35.0 35.0 35.0 42.0 42.0 42.0

moves_per_step:1
mcsteps:1
nprint:1
mc_move_prob:1 0 0 0
""",
    )

    res = subprocess.run(
        [
            gcmc_cpu,
            "--inp",
            str(inp),
            "--prefix",
            str(out_prefix),
            "--seed",
            "5",
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
    assert res.returncode == 0, res.stdout + res.stderr

    out_pdb = Path(f"{out_prefix}_final.pdb")
    assert out_pdb.exists()
    assert accept_log.exists()

    rec = _first_record(accept_log, move="insertion", species="RES")
    assert bool(rec.get("accepted", False)) is True
    delta_u = float(rec["deltaU"])

    box_nm = _read_cryst1_box_nm(out_pdb)
    parent = _pdb_atom_xyz_angstrom(out_pdb, resname="RES", atom_name="P1")
    drude = _pdb_atom_xyz_angstrom(out_pdb, resname="RES", atom_name="DP1")
    qins = _pdb_atom_xyz_angstrom(out_pdb, resname="FAR", atom_name="QI")

    r_pe = _minimum_image_dist_nm(parent, qins, box_nm)
    r_de = _minimum_image_dist_nm(drude, qins, box_nm)
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

    assert delta_u == pytest.approx(expected, abs=5e-3, rel=1e-6)

    # Host Drude must also be present and relaxed (proves hostDrudeParticles path is active).
    assert _pdb_contains_resname(out_pdb, "HST")
    host_p1 = _pdb_atom_xyz_angstrom(out_pdb, resname="HST", atom_name="O1")
    host_d1 = _pdb_atom_xyz_angstrom(out_pdb, resname="HST", atom_name="D1")
    host_disp_a = _minimum_image_dist_nm(host_p1, host_d1, box_nm) * 10.0  # nm -> Å
    assert host_disp_a > 0.02
    assert host_disp_a < 5.0


def test_drude_fragment_delete_deltaU_is_negative_of_insert_with_drude_host_present(gcmc_cpu, temp_dir):
    work = Path(temp_dir) / "drude_host_gcmc" / "insert_delete_symmetry"
    work.mkdir(parents=True, exist_ok=True)

    psf, pdb, prm, frag_itp, monomerdir = _setup_drude_host_and_insert_region(work)

    ins_prefix = work / "out_ins" / "gcmc"
    ins_prefix.parent.mkdir(parents=True, exist_ok=True)
    ins_accept = work / "out_ins" / "accept.jsonl"

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

    gcmc_region:box 35.0 35.0 35.0 42.0 42.0 42.0

moves_per_step:1
mcsteps:2
nprint:1
mc_move_prob:1 0 0 0
max_translation:0.0
max_rotation:0.0
""",
    )

    res_ins = subprocess.run(
        [
            gcmc_cpu,
            "--inp",
            str(inp_ins),
            "--prefix",
            str(ins_prefix),
            "--seed",
            "17",
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
    assert abs(delta_u_ins) > 1e-4

    checkpoint = Path(f"{ins_prefix}_checkpoint_1.dat")
    assert checkpoint.exists()

    del_prefix = work / "out_del" / "gcmc"
    del_prefix.parent.mkdir(parents=True, exist_ok=True)
    del_accept = work / "out_del" / "accept.jsonl"

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
# Make deletion acceptance deterministic by driving activity (z) very low.
# This does not affect ΔU symmetry (energy-only), but ensures the molecule is actually removed.
fragconc:0.001
fragmuex:0.0

box_size:50.0 50.0 50.0
cutoff:12.0
temperature:300.0
use_cavity_bias:no

    gcmc_region:box 35.0 35.0 35.0 42.0 42.0 42.0

moves_per_step:1
mcsteps:1
nprint:1
mc_move_prob:0 1 0 0
max_translation:0.0
max_rotation:0.0
""",
    )

    res_del = subprocess.run(
        [
            gcmc_cpu,
            "--inp",
            str(inp_del),
            "--prefix",
            str(del_prefix),
            "--seed",
            "19",
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
    assert bool(del_rec.get("accepted", False)) is True
    delta_u_del = float(del_rec["deltaU"])

    assert delta_u_del == pytest.approx(-delta_u_ins, rel=1e-4, abs=5e-3)

    out_pdb = Path(f"{del_prefix}_final.pdb")
    assert out_pdb.exists()
    assert not _pdb_contains_resname(out_pdb, "RES")

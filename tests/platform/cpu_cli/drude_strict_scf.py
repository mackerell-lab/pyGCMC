"""
gcmc_cpu Drude (polarizable) strict SCF integration tests.

These tests are file-driven and avoid stdout/stderr log-string assertions.
"""

from __future__ import annotations

import math
import subprocess
from pathlib import Path

import pytest


def _write_text(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content.strip() + "\n")


def _pdb_atom_xyz(pdb_path: Path, *, resname: str, atom_name: str, resid: int) -> tuple[float, float, float]:
    for line in pdb_path.read_text().splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        parts = line.split()
        # Expected minimal PDB tokenization:
        # ATOM serial name resname chain resid x y z ...
        if len(parts) < 9:
            continue
        if parts[2].upper() != atom_name.upper():
            continue
        if parts[3].upper() != resname.upper():
            continue
        if int(parts[5]) != resid:
            continue
        x = float(parts[6])
        y = float(parts[7])
        z = float(parts[8])
        return (x, y, z)
    raise AssertionError(f"Atom not found in {pdb_path}: {resname}{resid} {atom_name}")


def _dist(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    dx = a[0] - b[0]
    dy = a[1] - b[1]
    dz = a[2] - b[2]
    return math.sqrt(dx * dx + dy * dy + dz * dz)


def _setup_minimal_drude_system(work: Path) -> tuple[Path, Path, Path, Path, Path]:
    """
    Build a minimal 2-residue system:
      RES: parent + Drude (bonded)
      FIX: external charge (fixed framework)

    Topology uses PSF EXT DRUDE format so alpha/thole can be read without extra files.
    """
    psf = work / "system.psf"
    pdb = work / "system.pdb"
    prm = work / "ff.prm"
    frag_itp = work / "res.itp"
    monomerdir = work / "monomerdir"
    monomerdir.mkdir(parents=True, exist_ok=True)
    frag_pdb = monomerdir / "RES.pdb"

    # Topology atom types must have the same leading element letter as the PDB atom names,
    # otherwise MolecularCombiner will reject the structure/topology match.
    parent_type = "P"
    drude_type = "D"
    ext_type = "Q"

    # Polarizable atom parameters (CHARMM units in PSF/PRM):
    # alpha is in Å^3 and will be converted to nm^3 internally (×1e-3).
    alpha_a3 = 1.0
    thole = 2.6

    # Small initial parent-drude separation to avoid degenerate fits/poses; must be << final induced displacement.
    drude_dx_a = 0.005

    _write_text(
        psf,
        f"""
PSF EXT DRUDE

       1 !NTITLE
REMARKS minimal Drude SCF integration test

       3 !NATOM
       1 SYS  1 RES  P1   {parent_type}   1.000000  12.0110  0  {thole:.3f}  {alpha_a3:.3f}
       2 SYS  1 RES  DP1  {drude_type}   -1.000000   0.4000  0  0.000  0.000
       3 SYS  2 FIX  Q1   {ext_type}      1.000000  12.0110  0  0.000  0.000

       1 !NBOND: bonds
       1       2
""",
    )

    # PDB coordinates are in Å (gcmc_gpu/CHARMM-style external units).
    # Place the external charge 3 Å away from the Drude parent along +x.
    cx, cy, cz = 25.0, 25.0, 25.0
    ext_dx_a = 3.0
    _write_text(
        pdb,
        f"""
CRYST1   50.000   50.000   50.000  90.00  90.00  90.00 P 1           1
ATOM      1  P1  RES A   1      {cx:8.3f}{cy:8.3f}{cz:8.3f}  1.00  0.00           C
ATOM      2 DP1  RES A   1      {cx + drude_dx_a:8.3f}{cy:8.3f}{cz:8.3f}  1.00  0.00           D
ATOM      3  Q1  FIX A   2      {cx + ext_dx_a:8.3f}{cy:8.3f}{cz:8.3f}  1.00  0.00           Q
END
""",
    )

    # Template geometry for RES (relative positions; FragmentLibrary centers to COM).
    _write_text(
        frag_pdb,
        f"""
ATOM      1  P1  RES A   1       0.000   0.000   0.000  1.00  0.00           C
ATOM      2 DP1  RES A   1       {drude_dx_a:6.3f}   0.000   0.000  1.00  0.00           D
END
""",
    )

    # Minimal PRM: LJ is set to zero to isolate electrostatics + Drude spring energy.
    # Alpha/Thole is included as a fallback; PSF should already provide it.
    _write_text(
        prm,
        f"""
* minimal PRM for Drude SCF CLI integration test
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
; name  nrexcl
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


def test_gcmc_cpu_drude_strict_scf_relaxes_and_limits_displacement(gcmc_cpu, temp_dir):
    """
    Phase-I Drude acceptance criteria for gcmc_cpu (direct cutoff):
      - Drude SCF must run (Drude-parent distance changes from the initial near-zero value).
      - With strict SCF (hard wall), displacement must stay within maxDrudeDistance (0.02 nm = 0.2 Å).
    """
    work = Path(temp_dir) / "drude_strict_scf"
    work.mkdir(parents=True, exist_ok=True)

    psf, pdb, prm, frag_itp, monomerdir = _setup_minimal_drude_system(work)

    out_prefix = work / "out" / "gcmc"
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    inp = work / "run.inp"
    _write_text(
        inp,
        f"""
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

# Translation-only with zero step: ensures no geometry/logical changes except Drude relaxation.
moves_per_step:1
mcsteps:1
nprint:1
mc_move_prob:0 0 1 0
max_translation:0.0
max_rotation:0.0
""",
    )

    result = subprocess.run(
        [gcmc_cpu, "--inp", str(inp), "--prefix", str(out_prefix), "--seed", "123"],
        cwd=str(work),
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert result.returncode == 0, result.stdout + result.stderr

    out_pdb = Path(f"{out_prefix}_final.pdb")
    assert out_pdb.exists()

    init_parent = _pdb_atom_xyz(pdb, resname="RES", atom_name="P1", resid=1)
    init_drude = _pdb_atom_xyz(pdb, resname="RES", atom_name="DP1", resid=1)
    final_parent = _pdb_atom_xyz(out_pdb, resname="RES", atom_name="P1", resid=1)
    final_drude = _pdb_atom_xyz(out_pdb, resname="RES", atom_name="DP1", resid=1)

    initial_dist_a = _dist(init_parent, init_drude)
    final_dist_a = _dist(final_parent, final_drude)

    # Must move a visible amount (avoid false positives from template geometry).
    assert final_dist_a > initial_dist_a + 0.02

    # Strict SCF hard wall: 0.02 nm = 0.2 Å.
    assert final_dist_a <= 0.21
